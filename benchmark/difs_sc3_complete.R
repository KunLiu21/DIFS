sc3_clustering <- function(data.mat, cluster_count, input_markers = NULL,
                           seed = 1L, n_cores = 2L) {
  data.mat <- as.matrix(data.mat)
  if (is.null(rownames(data.mat)) || is.null(colnames(data.mat)) ||
      anyDuplicated(colnames(data.mat))) stop("Gene/cell identifiers are required; cell IDs must be unique")
  data.mat <- data.mat[!duplicated(rownames(data.mat)), , drop = FALSE]
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = data.mat, logcounts = log2(data.mat + 1)))
  SummarizedExperiment::rowData(sce)$feature_symbol <- rownames(sce)
  sce <- SC3::sc3_prepare(sce, n_cores = as.integer(n_cores),
                         rand_seed = as.integer(seed))
  if (!inherits(sce, "SingleCellExperiment")) stop("sc3_prepare failed")
  if (!is.null(input_markers)) {
    ix <- match(unique(input_markers), rownames(sce))
    if (!length(ix) || anyNA(ix)) stop("Empty or unmatched input markers")
    col_sds <- matrixStats::colVars(data.mat[ix, , drop = FALSE], useNames = TRUE)
    if (anyNA(col_sds)) stop("Undefined cell variance in the selected feature set")
    for (col_id in which(col_sds == 0)) ix <- c(ix, which.max(data.mat[, col_id]))
    ix <- unique(ix)
    filt <- rep(FALSE, nrow(sce)); filt[ix] <- TRUE
    SummarizedExperiment::rowData(sce)$sc3_gene_filter <- filt
  }
  markers <- rownames(sce)[SummarizedExperiment::rowData(sce)$sc3_gene_filter]
  if (is.null(cluster_count)) {
    sce <- SC3::sc3_estimate_k(sce)
    cluster_count <- S4Vectors::metadata(sce)$sc3$k_estimation
  }
  if (length(cluster_count) != 1L || !is.finite(cluster_count) ||
      cluster_count < 2 || cluster_count != as.integer(cluster_count)) stop("Invalid cluster_count")
  cluster_count <- as.integer(cluster_count)
  sce <- SC3::sc3_calc_dists(sce)
  sce <- SC3::sc3_calc_transfs(sce)
  sce <- SC3::sc3_kmeans(sce, ks = cluster_count)
  sce <- SC3::sc3_calc_consens(sce)
  label_col <- paste0("sc3_", cluster_count, "_clusters")
  if (!label_col %in% names(SummarizedExperiment::colData(sce))) stop("SC3 label column missing")
  before <- SummarizedExperiment::colData(sce)[[label_col]]
  names(before) <- colnames(sce)
  meta <- S4Vectors::metadata(sce)$sc3
  train_ix <- meta$svm_train_inds
  study_ix <- meta$svm_study_inds
  if (length(train_ix) && length(study_ix)) sce <- SC3::sc3_run_svm(sce, ks = cluster_count)
  cluster <- SummarizedExperiment::colData(sce)[[label_col]]
  names(cluster) <- colnames(sce)
  if (!identical(names(cluster), colnames(data.mat))) stop("Cell order changed")
  if (anyNA(cluster)) stop("SC3 still has unassigned cells after hybrid completion")
  if (length(train_ix) && !identical(as.character(before[train_ix]), as.character(cluster[train_ix])))
    stop("SVM completion changed training-cell labels")
  list(cluster = cluster, markers = markers, cluster_before_svm = before,
       svm_train_cells = colnames(sce)[train_ix],
       svm_study_cells = colnames(sce)[study_ix],
       audit = list(version = "sc3-complete-v1", SC3_version = as.character(utils::packageVersion("SC3")),
                    seed = as.integer(seed), n_cores = as.integer(n_cores),
                    k = cluster_count, n_total = ncol(sce),
                    n_training = length(train_ix), n_prediction = length(study_ix),
                    n_missing_before = sum(is.na(before)), n_missing_after = sum(is.na(cluster)),
                    n_input_markers = length(unique(input_markers)), n_effective_markers = length(markers)))
}
