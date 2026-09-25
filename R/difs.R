# Public API. Scientific functions are read from the shared benchmark sources.
.difs_api_file <- tryCatch(sys.frame(1)$ofile,error=function(e) NULL)
if (is.null(.difs_api_file)) stop("Load this API with source('R/difs.R')")
.difs_repo <- normalizePath(file.path(dirname(.difs_api_file),".."),mustWork=TRUE)

difs_require <- function(packages) {
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop("Required packages unavailable: ",paste(missing,collapse=", "),
    ". See environment/README.md; no pipeline step was skipped.")
}
difs_load_definitions <- function(path,wanted,env) {
  found <- character()
  for (expr in parse(path,keep.source=FALSE)) {
    if (is.call(expr) && as.character(expr[[1]]) %in% c("<-","=") &&
        is.symbol(expr[[2]]) && as.character(expr[[2]]) %in% wanted) {
      name <- as.character(expr[[2]])
      if (name %in% found) stop("Duplicate definition: ",name)
      eval(expr,envir=env); found <- c(found,name)
    }
  }
  if (!setequal(found,wanted)) stop("Missing definitions in ",basename(path),": ",paste(setdiff(wanted,found),collapse=", "))
  invisible(env)
}
difs_load_engine <- function(root=.difs_repo,full=TRUE) {
  difs_require("diptest"); e <- new.env(parent=globalenv())
  difs_load_definitions(file.path(root,"benchmark/Functions_mixture_auto.R"),
    c("whichPC","dip.test.v1","binary.bound.decider","transform.binary","splited.binary.mat",
      "excute.fisher.test","fisher_test","interleave","seurat_k_cluster","seudo_clustering"),e)
  difs_load_definitions(file.path(root,"benchmark/Functions_controlled.R"),
    c("num_features_decider","DIFS_RES_GRID","DIFS_TOPS_GRID","DIFS_MIX_RATIO","DIFS_GATE_RULES",
      "difs_gate_threshold","difs_gate_genes","difs_null_path",".difs_null_cache","difs_null_tables",
      "dip_pvalue_lookup_vec","DIFS_STAGE1_KEYS","difs_stage1_ranking","difs_features"),e)
  difs_load_definitions(file.path(root,"benchmark/difs_sc3_complete.R"),"sc3_clustering",e)
  e$difs_null_path <- local({table <- file.path(root,"inst/extdata/dip_null_tables.rds")
    function() {if (!file.exists(table)) stop("Missing bundled null table: ",table);table}})
  if (full) {
    difs_require(c("Seurat","SeuratObject","SC3","SingleCellExperiment","SummarizedExperiment",
      "S4Vectors","matrixStats","FEAST","magrittr"))
    for (name in c("VariableFeatures<-","ScaleData","RunPCA","FindNeighbors","FindClusters","AddMetaData","SplitObject"))
      e[[name]] <- getExportedValue("Seurat",name)
    e$`%>%` <- getExportedValue("magrittr","%>%")
    e$cal_MSE <- getExportedValue("FEAST","cal_MSE")
    # API adaptation only: read the same RNA assay under SeuratObject v4/v5.
    e$GetAssayData <- function(object,slot="data",...) {
      if (utils::packageVersion("SeuratObject") >= "5.0.0") SeuratObject::GetAssayData(object,layer=slot,...)
      else SeuratObject::GetAssayData(object,slot=slot,...)
    }
  }
  e
}
difs_validate_counts <- function(counts) {
  if (length(dim(counts))!=2L || nrow(counts)<2L || ncol(counts)<3L) stop("counts must be a genes-by-cells matrix")
  if (is.null(rownames(counts)) || is.null(colnames(counts)) || anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
      any(!nzchar(rownames(counts))) || any(!nzchar(colnames(counts))) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) stop("Unique, nonempty gene and cell identifiers are required")
  vals <- if (inherits(counts,"sparseMatrix")) counts@x else as.vector(counts)
  if (!is.numeric(vals)) stop("counts must contain numeric values")
  if (any(!is.finite(vals))) stop("counts contain NA, NaN or infinite values")
  if (any(vals<0)) stop("counts contain negative values")
  # Preserve fractional count estimates, as the benchmark preparation does.
  # Do not round, rescale or replace the input to satisfy validation.
  libs <- if (inherits(counts,"sparseMatrix")) Matrix::colSums(counts) else colSums(counts)
  if (any(libs<=0)) stop("Every cell must have a positive library size")
  invisible(TRUE)
}
difs_check_integer <- function(x,name,lower=1L,upper=.Machine$integer.max) {
  if (length(x)!=1L || !is.numeric(x) || !is.finite(x) || x!=floor(x) || x<lower || x>upper)
    stop(name," must be one integer in [",lower,", ",upper,"]")
  as.integer(x)
}
difs_fit <- function(counts,k,n_features=100L,seed=1L,root=.difs_repo) {
  difs_validate_counts(counts)
  k <- difs_check_integer(k,"k",2L,ncol(counts)-1L)
  n_features <- difs_check_integer(n_features,"n_features"); seed <- difs_check_integer(seed,"seed",0L)
  e <- difs_load_engine(root,full=TRUE)
  previous <- options(Seurat.object.assay.version="v3"); on.exit(options(previous),add=TRUE)
  set.seed(seed)
  seu <- Seurat::CreateSeuratObject(counts,min.cells=0,min.features=0)
  if (!identical(colnames(seu),colnames(counts)) || !identical(rownames(seu),rownames(counts)))
    stop("Seurat changed identifiers; prepare a stable input first")
  seu <- Seurat::NormalizeData(seu,normalization.method="LogNormalize",scale.factor=1e4,verbose=FALSE)
  logmat <- as.matrix(e$GetAssayData(seu,slot="data"))
  ranking <- e$difs_stage1_ranking(logmat,gate_rule="submitted",gate_k=k,key="normal_lookup",seed=seed)
  trace <- new.env(parent=emptyenv()); repaired <- e$sc3_clustering
  # Preserve the historical repaired-wrapper defaults: seed=1, n_cores=2.
  e$sc3_clustering <- function(...) {z <- repaired(...);trace$preliminary <- z;z}
  started <- proc.time()[["elapsed"]]
  selected <- e$difs_features(logmat,seu,cluster_count=k,gate_rule="submitted",gate_k=k,
    tune_stage1=TRUE,tune_ratio=TRUE,stage1_key="normal_lookup",use_stage2=TRUE,seed=seed,n_fixed=n_features)
  if (is.null(trace$preliminary) || is.null(selected$stage2_n) || is.na(selected$stage2_n)) stop("The stage-II path was not executed")
  if (!length(selected$markers) || anyNA(selected$markers) || anyDuplicated(selected$markers) ||
      !all(selected$markers %in% rownames(counts))) stop("Invalid returned gene set")
  clustered <- e$seudo_clustering(seu,markers=selected$markers,cluster_count=k,
    clustering.method="Refined Louvain",res=e$DIFS_RES_GRID)
  labels <- setNames(as.character(clustered$seurat_clusters),colnames(clustered))
  if (!identical(names(labels),colnames(counts)) || anyNA(labels) || any(!nzchar(labels))) stop("Final labels do not cover all input cells in order")
  list(features=selected$markers,stage1_ranking=ranking,preliminary_labels=trace$preliminary$cluster,labels=labels,
    stage1_n=selected$stage1_n,stage2_n=selected$stage2_n,ratio=selected$ratio,
    n_features_realised=length(selected$markers),n_features_requested=n_features,gate_n=length(ranking),
    k_input=k,n_clusters_found=length(unique(labels)),seed=seed,sc3_audit=trace$preliminary$audit,
    stage2_executed=TRUE,elapsed_selection_and_final_clustering_sec=proc.time()[["elapsed"]]-started,
    gate_rule="submitted",session=utils::sessionInfo())
}
difs_metrics <- function(truth,predicted) {
  if (is.null(names(truth)) || is.null(names(predicted)) || anyDuplicated(names(truth)) ||
      anyDuplicated(names(predicted)) || !setequal(names(truth),names(predicted))) stop("Labels require matching unique cell IDs")
  truth <- as.character(truth[names(predicted)]); predicted <- as.character(predicted)
  if (anyNA(truth) || anyNA(predicted) || any(!nzchar(truth)) || any(!nzchar(predicted))) stop("Missing labels cannot be silently excluded")
  tab <- table(truth,predicted); n <- sum(tab); if (n<2) stop("At least two labeled cells are required")
  tp <- sum(choose(tab,2)); a <- sum(choose(rowSums(tab),2)); b <- sum(choose(colSums(tab),2))
  expectation <- a*b/choose(n,2); denominator <- (a+b)/2-expectation
  c(ARI=if (denominator==0) 1 else (tp-expectation)/denominator,
    FM=if (a*b==0) NA_real_ else tp/sqrt(a*b),Jaccard=if (a+b-tp==0) NA_real_ else tp/(a+b-tp),
    Purity=sum(apply(tab,2,max))/n)
}
