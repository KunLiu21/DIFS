# difs_rev2_overlay.R: revision-2 analysis source.
# Source: code/as_run/02_core/difs_rev2_overlay.R
# Usage, scope and limitations: benchmark/REVISION2.md.
if (!exists("difs_stage1_ranking") || !exists("interleave") ||
    !exists("seudo_clustering") || !exists("sc3_clustering"))
  stop("difs_rev2_overlay.R must be sourced after Functions_controlled.R and difs_sc3_complete.R")
if (!"seed" %in% names(formals(sc3_clustering)))
  stop("sc3_clustering has no seed argument -- source difs_sc3_complete.R before this overlay")

.difs_orig_stage1      <- difs_stage1_ranking
.difs_orig_interleave  <- interleave
.difs_orig_seudo       <- seudo_clustering

difs_stage1_ranking <- function(data.log.mat, min.expression = log(5),
                                min.expression.proportion = 0.05,
                                gate_rule = "submitted", gate_k = NULL,
                                key = "normal_lookup", B = 2000, seed = 1) {
  src <- getOption("difs.stage1_source", "dip")
  if (identical(src, "dip"))
    return(.difs_orig_stage1(data.log.mat, min.expression, min.expression.proportion,
                             gate_rule, gate_k, key, B, seed))
  if (!identical(src, "hvg")) stop("unknown difs.stage1_source: ", src)

  g <- difs_gate_genes(data.log.mat, min.expression, min.expression.proportion,
                       gate_rule = gate_rule, gate_k = gate_k)
  if (!length(g)) stop("no gene passed the stage I expression gate")
  seu <- get0("DIFS_CURRENT_SEU", envir = globalenv(), inherits = FALSE)
  if (is.null(seu)) stop("HVG stage I needs DIFS_CURRENT_SEU (set by Method_comparison_rev2.R)")
  if (ncol(seu) != ncol(data.log.mat) || !all(g %in% rownames(seu)))
    stop("DIFS_CURRENT_SEU does not match the matrix being ranked")
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst",
                                      nfeatures = nrow(seu), verbose = FALSE)
  hv  <- Seurat::HVFInfo(seu)
  col <- grep("variance\\.standardized$", colnames(hv), value = TRUE)[1]
  if (is.na(col)) stop("HVFInfo has no standardised-variance column; found: ",
                       paste(colnames(hv), collapse = ", "))
  v <- stats::setNames(hv[[col]], rownames(hv))
  if (anyNA(v[g])) stop(sum(is.na(v[g])), " gated gene(s) have no vst value")
  g[order(v[g], decreasing = TRUE)]
}

interleave <- function(vec1, vec2, ratio) {
  if (identical(getOption("difs.final_set", "mixed"), "stage2only")) return(vec1)
  .difs_orig_interleave(vec1, vec2, ratio)
}

seudo_clustering <- function(data.seu, markers, cluster_count, clustering.method,
                             res = c(seq(0.01, 0.1, 0.01), seq(0.1, 1.5, 0.05), seq(1.5, 3, 0.1)),
                             cds = NULL) {
  s <- getOption("difs.cluster_seed")
  if (is.null(s) || !clustering.method %in% c("Refined Louvain", "SLM", "SC3"))
    return(.difs_orig_seudo(data.seu, markers, cluster_count, clustering.method, res, cds))
  s <- as.integer(s)
  VariableFeatures(data.seu) <- markers
  data.seu <- ScaleData(data.seu, verbose = FALSE)
  data.seu <- RunPCA(data.seu, verbose = FALSE, seed.use = s)
  nPC <- whichPC(data.seu)
  if (clustering.method %in% c("Refined Louvain", "SLM")) {
    data.seu <- FindNeighbors(data.seu, verbose = FALSE, dims = 1:nPC)
    data.seu <- FindClusters(data.seu, resolution = res, verbose = FALSE,
                             algorithm = if (clustering.method == "SLM") 3 else 2,
                             random.seed = s)
    data.seu <- seurat_k_cluster(data.seu, cluster_count = cluster_count)
  } else {
    data.mat <- as.matrix(GetAssayData(data.seu, slot = "counts"))
    sc3_cluster <- sc3_clustering(data.mat, cluster_count = cluster_count, markers, seed = s)
    data.seu <- AddMetaData(data.seu, sc3_cluster$cluster, "seurat_clusters")
  }
  data.seu
}

sc3_estimate_k <- function(object) {
  k <- getOption("difs.k_override")
  if (is.null(k)) return(SC3::sc3_estimate_k(object))
  md <- S4Vectors::metadata(object)
  if (is.null(md$sc3)) md$sc3 <- list()
  md$sc3$k_estimation <- as.integer(k)
  S4Vectors::metadata(object) <- md
  object
}

difs_rev2_options_reset <- function()
  options(difs.cluster_seed = NULL, difs.stage1_source = "dip",
          difs.final_set = "mixed", difs.k_override = NULL)
difs_rev2_options_reset()
