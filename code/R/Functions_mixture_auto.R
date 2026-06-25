suppressPackageStartupMessages({
  library(monocle)
  library(scran)
  library(Seurat)
  library(diptest)
  library(SC3)
  library(FEAST)
  library(dplyr)
  library(glmGamPoi)
  library(SingleCellExperiment)})

prepDir <- function(path){
  if (!dir.exists(path)) {
    dir.create(path, showWarnings = F, recursive = T)
  }
  return(path)
}

whichPC <- function(seurat){
  pct <- seurat[["pca"]]@stdev / sum(seurat[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  pcs <- min(co1, co2)
  return(pcs)
}

######################################################################################################################################
###################################################Preprocessing####################################################################
######################################################################################################################################

preprocessingData <- function(count.mat, min.cells = 3, min.features = 200){
  #using the default setting
  seurat <-  CreateSeuratObject(count.mat, min.cells = min.cells, min.features=min.features)
  #End at FindVariableFeatures to make sure the seurat object has the Variable Features slot
  seurat <- seurat %>% SCTransform(verbose = F)
  return(seurat)
}

######################################################################################################################################
###################################################Doing Dip tests####################################################################
######################################################################################################################################


dip.test.v1 <- function (x, simulate.p.value = T, B = 2000){
  x <- sort(x[complete.cases(x)])
  stopifnot(is.numeric(x))
  n <- length(x)
  D <- dip(x)
  P <- mean(D <= replicate(B, dip(rnorm(n, mean(x), sd(x)))))
  return(P)
}


## In the dip tests, the binary bound is set to filter out the low expression counts.
dip_test_pip.v1 <- function(data.log.mat, gene, min.expression = 0.5, min.expression.proportion = 0.05){
  thresh1 = min.expression.proportion * dim(data.log.mat)[2]
  test_gene <- data.log.mat[which(rownames(data.log.mat) == gene), ]
  if(suppressWarnings(max(test_gene) <= min.expression)){
    #If the conditions are not met then the "p-value" will be 1.1 which means We don't do the test in the case.
    return(1.1)
  }else{
    gene.log <- test_gene[test_gene>min.expression]
    #checking the expressed genes' expression
    if (length(gene.log) < max(thresh1, 35) ){
      return(1.1)
    } else{
      return(dip.test.v1(gene.log))
    }
  }
}

preprocessing_dip_results <- function(data.log.mat, min.expression = 0.5, min.expression.proportion ){
  genenames <- rownames(data.log.mat)
  p.value.dip = suppressWarnings(sapply(1:length(genenames), FUN = function(x) dip_test_pip.v1(data.log.mat, genenames[x] , min.expression= min.expression, min.expression.proportion )))
  names(p.value.dip) <- genenames
  sorted_p <- sort(p.value.dip)
  return(sorted_p)
}



######################################################################################################################################
###################################################Fisher exact test method###########################################################
######################################################################################################################################
binary.bound.decider <- function(data.log.mat){
  #binary.bound<- quantile(data.log.mat[data.log.mat>0], 0.5)
  ##the binary bound is used to do Fisher's exact test. It is a indicator of high expression
  binary.bound = log(5)
  return(binary.bound)
}


transform.binary <- function(data.log.mat, binary.bound) {
  result <- ifelse(data.log.mat > binary.bound, 1, 0)
  return(result)
}


splited.binary.mat <- function(data.seu, cluster.ident, binary.bound){
  splited_seu <- SplitObject(data.seu, split.by = cluster.ident)
  num_clusters <- length(splited_seu)
  splited_mat <- lapply(1: num_clusters, function(x) as.matrix(GetAssayData(splited_seu[[x]], slot = "data")))
  splited_binary_mat <- lapply(1:num_clusters, function(x) transform.binary(splited_mat[[x]], binary.bound))
  return(splited_binary_mat)
}


excute.fisher.test <- function(vector1, vector2) {
  levels <- 0:1
  vector1 <- factor(vector1, levels = levels)
  vector2 <- factor(vector2, levels = levels)
  # Count occurrences of 1s and 0s in each vector
  count_vector1 <- table(vector1)
  count_vector2 <- table(vector2)
  
  # Create a 2x2 contingency table
  contingency_table <- matrix(c(count_vector1[1], count_vector2[1],
                                count_vector1[2], count_vector2[2]), nrow = 2)
  result <- fisher.test(contingency_table)
  
  return(result$p.value)
}

fisher_test <- function(data.seu, cluster.ident, binary.bound, min.mean.expr = 0.01){
  splited_binary_mat <- splited.binary.mat(data.seu, cluster.ident, binary.bound)
  # total_expressed_counts <- lapply(1: length(splited_binary_mat), function(x) apply(splited_binary_mat[[x]], 1, sum))
  # total_expressed_mat <- t(do.call(rbind, total_expressed_counts))
  mean_expressed_counts <- lapply(1: length(splited_binary_mat), function(x) apply(splited_binary_mat[[x]], 1, mean))
  mean_expressed_mat <- t(do.call(rbind, mean_expressed_counts))
  
  ## the genes will be tested as the stage II markers only in the clusters where its mean expression is the highest.
  cluster_testing <- data.frame(clusters = apply(mean_expressed_mat, 1, which.max))
  
  
  p_values <- sapply(1:length(cluster_testing$clusters), function(i) {
    if(max(mean_expressed_mat[i, ]) < min.mean.expr){
      return(p = 1.1)
    }
    vector2 = c()
    p = excute.fisher.test(
      vector1 = splited_binary_mat[[cluster_testing$clusters[i]]][rownames(cluster_testing)[i], ],
      vector2 = unlist(sapply((1:length(splited_binary_mat))[-cluster_testing$clusters[i]], function(x) {
        c(vector2, splited_binary_mat[[x]][rownames(cluster_testing)[i], ])
      }))
    )
    return(p)
  })
  names(p_values) = rownames(data.seu)
  return(sort(p_values))
}

############################################################mix two stages markers #############################################
### Mix stage 1 and stage 2 markers one by one.
interleave <- function(vec1, vec2, ratio) {
  
  loops <- sum(ratio)
  n <- max(length(vec1), length(vec2))
  vec1_new <- rep(vec1, length.out = n*ratio[1])
  vec2_new <- rep(vec2, length.out = n*ratio[2])
  mixed_vec <- numeric(n*loops)
  
  for(i in 1: (n*loops)){
    if(i %% loops == 1){
      mixed_vec[i: (i-1+ratio[1])] <- vec1_new[(i%/%loops * ratio[1]+1): (i%/%loops * ratio[1] + ratio[1])]
      mixed_vec[(i+ratio[1]): (i+loops-1)] <- vec2_new[(i%/%loops * ratio[2] +1): (i%/%loops * ratio[2] + ratio[2])]
    }
  }
  return(unique(mixed_vec))
}





######################################################################################################################################
###################################################clustering methods##################################################################
######################################################################################################################################

sc3_clustering <- function (data.mat, cluster_count, input_markers = NULL) {
  sce = SingleCellExperiment(assays = list(counts = as.matrix(data.mat), 
                                           logcounts = log2(as.matrix(data.mat) + 1)))
  rowData(sce)$feature_symbol <- rownames(sce)
  sce = sce[!duplicated(rowData(sce)$feature_symbol), ]
  if (!is.null(input_markers)) {
    rowD = rowData(sce)
    sce = sc3_prepare(sce, n_cores = 2)
    sc3_gene_filter = rep(FALSE, nrow(sce))
    ix = match(input_markers, rownames(sce))
    col_sds = colVars(data.mat[ix, ],useNames = T)
    col_ix = which(col_sds == 0)
    if (length(col_ix) > 0) {
      for (col_id in col_ix) {
        id_add = which.max(data.mat[, col_id])
        ix = c(ix, id_add)
      }
    }
    markers = rownames(data.mat)[ix]
    sc3_gene_filter[ix] = TRUE
    rowD$sc3_gene_filter = sc3_gene_filter
    markers = rownames(rowD)[rowD$sc3_gene_filter == TRUE]
    rowData(sce) = rowD
  }
  else {
    sce = sc3_prepare(sce, n_cores = 1)
    rowD = rowData(sce)
    markers = rownames(rowD)[rowD$sc3_gene_filter == TRUE]
  }
  if (is.null(cluster_count)) {
    sce = sc3_estimate_k(sce)
    cat("estimated k clusters: ", metadata(sce)$sc3$k_estimation, 
        "\n")
    cluster_count = metadata(sce)$sc3$k_estimation
  }
  sce = sc3_calc_dists(sce)
  sce = sc3_calc_transfs(sce)
  sce = sc3_kmeans(sce, ks = cluster_count)
  sce = sc3_calc_consens(sce)
  colTb = data.frame(colData(sce), stringsAsFactors = FALSE)
  cluster = colTb[, 1]
  names(cluster) = rownames(colTb)
  return(list(cluster = cluster, markers = markers))
}


monocle_clustering <- function(cds, cluster_count, input_markers, nPC) {
  cds = cds[input_markers, ]
  cds <- reduceDimension(cds, max_components=2, norm_method = 'log', reduction_method = 'tSNE', num_dim = nPC,  verbose = F)
  cds <- monocle::clusterCells(cds, verbose = F, num_clusters = cluster_count)
  cluster <- cds$Cluster
  names(cluster) <- colnames(cds)
  return(cluster)
}

seurat_k_cluster = function(data.seu, cluster_count, res = c(seq(0.01, 0.1, 0.01), seq(0.1, 1.5, 0.05), seq(1.5, 3, 0.1))){
  cluster_lists <- paste("RNA_snn_res.",res, sep = "")
  num.clusters = numeric()
  k = 1
  for(i in cluster_lists){
    num.clusters[k] <- length(levels(data.seu@meta.data[, i]))
    k= k+1
  }
  data.seu <- AddMetaData(data.seu, data.seu@meta.data[, cluster_lists[which.min(abs(num.clusters - cluster_count))]], "seurat_clusters")
  return(data.seu)
}

seudo_clustering <- function(data.seu, markers, cluster_count, clustering.method, res = c(seq(0.01, 0.1, 0.01), seq(0.1, 1.5, 0.05), seq(1.5, 3, 0.1)), cds = NULL){
  VariableFeatures(data.seu) = markers
  data.seu <- data.seu %>% ScaleData(verbose = F) %>%RunPCA(verbose = F)
  nPC= whichPC(data.seu)
  if(clustering.method %in% c("Refined Louvain","SLM")){
    data.seu <- FindNeighbors(data.seu, verbose =F, dims = 1:nPC)
    if(clustering.method == "Refined Louvain"){
      data.seu <- FindClusters(data.seu, resolution = res, verbose = F, algorithm = 2)
    }else if(clustering.method == "SLM"){
      data.seu <- FindClusters(data.seu, resolution = res, verbose = F, algorithm = 3)
    }
    data.seu = seurat_k_cluster(data.seu, cluster_count = cluster_count)
  }else if(clustering.method == "SC3"){
    data.mat <- as.matrix(GetAssayData(data.seu, slot = "counts"))
    sc3_cluster <- sc3_clustering(data.mat, cluster_count = cluster_count, markers)
    data.seu <- AddMetaData(data.seu, sc3_cluster$cluster , "seurat_clusters")
  }else if(clustering.method == "TSCAN"){
    data.mat <- as.matrix(GetAssayData(data.seu, slot = "data"))
    tscan_cluster <- TSCAN_Clust(data.mat, k = cluster_count, input_markers = markers)$cluster
    data.seu <- AddMetaData(data.seu, tscan_cluster , "seurat_clusters")
  }else if(clustering.method =="monocle"){
    monocle_cluster <-suppressWarnings( monocle_clustering(cds, cluster_count  = cluster_count, input_markers = markers, nPC))
    data.seu <- AddMetaData(data.seu, monocle_cluster , "seurat_clusters")
  }
    return(data.seu)
}


###########################################Optimal feature set##################################

num_features_decider <- function (data.seu, cluster_count, tops = seq(500, 2500, 100), top_genes){
  data.log.mat = as.matrix(GetAssayData(data.seu, slot = "data"))
  ntop = length(tops)
  mse = sc3_res = NULL
  for (i in seq_len(ntop)) {
    top = tops[i]
    markers = top_genes[1:top]
    data.seu = seudo_clustering(data.seu, markers = markers, cluster_count = cluster_count, clustering.method = "Refined Louvain")
    sc3_rslt <- list(cluster = data.seu$seurat_clusters, markers = markers)
    sc3_res[[toString(top)]] = sc3_rslt
    mse = c(mse, cal_MSE(as.matrix(data.log.mat), sc3_rslt$cluster))
  }
  names(mse) = tops
  res = list(mse = mse, sc3_res = sc3_res)
  which.num.features <- which.min(mse)
  optimal_markers <- res$sc3_res[[which.num.features]]$markers
  return(list(markers = optimal_markers, mse = min(mse)))
}

###########################Feature selection methods ###############################################
monocle_gene_filter <- function(cds, top_gene_count = 2000) {
  disp_table <- dispersionTable(cds)
  ordered_disp_table <- disp_table[order(disp_table$dispersion_fit, decreasing = F), ]
  top_genes <- head(ordered_disp_table, top_gene_count)$gene_id
  return(top_genes)
}


SC3_filter_genes <- function(data.log.mat,  lower_threshold = 0.06, upper_threshold =0.94) {
  dropout_rates <- rowSums(data.log.mat == 0) / ncol(data.log.mat)
  filtered_genes <- names(dropout_rates)[dropout_rates > lower_threshold & dropout_rates < upper_threshold]
  return(filtered_genes)
}


gm_markers_finder <- function(data.log.mat, cluster_count, binary.bound= binary.bound,
                              min.expression.proportion = 0.05, data.seu){
  genelist <- rownames(data.log.mat)
  gene_ranking <- preprocessing_dip_results(data.log.mat = data.log.mat, binary.bound,
                                            min.expression.proportion)
  stage1_markers <- num_features_decider(data.seu, cluster_count = cluster_count, top_genes = names(gene_ranking))$markers
  
  #stage1_markers = names(gene_ranking)[1:2000]
  data.seu <- seudo_clustering(data.seu, markers = stage1_markers, cluster_count, res = c(seq(0.01, 0.1, 0.01), seq(0.1, 1.5, 0.05), seq(1.5, 3, 0.1)),
                               clustering.method = "SC3") #1_30
  if( any(is.na(data.seu$seurat_clusters))){ 
    filtered.seu <- data.seu[, -which(is.na(data.seu$seurat_clusters))] 
  } else{
    filtered.seu <- data.seu
  }
  sorted_p_values <- fisher_test(filtered.seu, "seurat_clusters", binary.bound = binary.bound, min.mean.expr = 0.1)
  stage2_markers <- names(sorted_p_values[sorted_p_values<0.05])
  if(length(stage2_markers) ==0){
    gm_markers = stage1_markers
  }else{
    ratio <- rbind(c(1, 1), c(2, 3), c(1, 2), c(1, 3), c(3, 2), c(2, 1), c(3,1))
    mixture_lists <- apply(ratio, 1, function(x){
      num_features_decider(filtered.seu, cluster_count = cluster_count, tops = 1000, top_genes =interleave(stage2_markers, stage1_markers, ratio = x)) 
    })
    
    gm_markers <- mixture_lists[[which.min(lapply(mixture_lists, function(x){
      return(x$mse)
    }))]]$markers
  }
  return(gm_markers)
}




########################################## Validation #############################################################


# gene_hist <- function(mat, genename, binary.bound.gene){
#   slcted_gene = mat[rownames(mat) == genename]
#   slcted_gene = slcted_gene[slcted_gene >binary.bound.gene]
#   return(slcted_gene)
# }




Run.method.comparison <- function(data.path, clustering.method, feature.selection.method){
  data.name = gsub(".*/([^.]*).*", "\\1", data.path)
  data.seu <- readRDS(file = data.path)
  trueclass <- data.seu$trueclass
  data.log.mat <-  as.matrix(GetAssayData(data.seu, slot = "data"))
  genelist <- rownames(data.seu)
  data.mat = as.matrix(GetAssayData(data.seu, slot = "count"))
  sce = SingleCellExperiment(assays = list(counts = as.matrix(data.mat),
                                           logcounts = log2(as.matrix(data.mat) + 1)))
  sce = sc3_estimate_k(sce)
  cluster_count = metadata(sce)$sc3$k_estimation
  
 
  cds <- suppressWarnings(convertTo(sce, type = "monocle"))
  cds = suppressWarnings(estimateSizeFactors(cds))
  cds <- suppressWarnings(estimateDispersions(cds))
  
  if(feature.selection.method == "gm"){
    binary.bound <- binary.bound.decider(data.log.mat = data.log.mat)
    markers <-gm_markers_finder(data.log.mat = data.log.mat, cluster_count = length(unique(trueclass)), binary.bound = binary.bound, data.seu = data.seu)
  }else if(feature.selection.method == "Seurat"){
    data.seu <- FindVariableFeatures(data.seu, nfeatures = 1000, verbose = F)
    markers <- VariableFeatures(data.seu)[1:1000]
  } else if(feature.selection.method == "Feast"){
    if(dim(data.log.mat )[2] >= 1000){
      ixs = FEAST_fast(data.log.mat , k=cluster_count, batch_size = 1000)
    }else{
      ixs = FEAST(data.log.mat , k=cluster_count)
    }
    feast_markers <- rownames(data.log.mat)[ixs]
    #markers = num_features_decider(data.seu,  cluster_count =cluster_count, top_genes = feast_markers)$markers
    markers = feast_markers[1:1000]
  }else if( feature.selection.method == "SC3"){
    sc3_markers <- SC3_filter_genes(data.log.mat,  lower_threshold = 0.06, upper_threshold =0.94)
    markers <- sample(sc3_markers, 1000, replace = F)
  }
  data.seu <- seudo_clustering(data.seu, markers, cluster_count = length(unique(trueclass)), clustering.method = clustering.method,  cds = cds)
  model_perform <- eval_Cluster(data.seu$seurat_clusters, trueclass)
  
  return(list(data.name = data.name, clustering.method =clustering.method,
              feature.selection.method = feature.selection.method, ARI = model_perform[1], Purity = model_perform[2], Jaccard = model_perform[3], FM =model_perform[4]))
}

