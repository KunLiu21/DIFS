# difs_export_features_v2.R: revision-2 analysis source.
# Source: code/as_run/07_feature_sets_gate_ranking/difs_export_features_v2.R
# Usage, scope and limitations: benchmark/REVISION2.md.
difs_args <- function(args, known, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    if (!kv[2] %in% known)
      stop("unknown argument '", kv[2], "'. This script accepts: ",
           paste(known, collapse = ", "))
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
  invisible(NULL)
}
task <- ""; dataset <- "Baron"; method <- ""
base_dir <- "../results/feature_sets"; out_dir <- "../results/feature_sets_v2"
n_fixed <- 1000; seed <- 1
difs_args(commandArgs(TRUE),
          c("task", "dataset", "method", "base_dir", "out_dir", "n_fixed", "seed"))
if (!task %in% c("select", "universe")) stop("task= must be 'select' or 'universe'")
if (normalizePath(out_dir, mustWork = FALSE) == normalizePath(base_dir, mustWork = FALSE))
  stop("out_dir must differ from base_dir: this script never writes over the v1 export")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_info_path <- file.path(base_dir, "feature_sets_info.csv")
if (!file.exists(base_info_path)) stop("v1 export info not found: ", base_info_path)
base_info <- utils::read.csv(base_info_path, stringsAsFactors = FALSE)

norm_id <- function(x) {
  x <- sub("--.*$", "", x)            # Muraro's "GENE--chrN"
  x <- sub("\\.[0-9]+$", "", x)       # Ensembl version suffix
  toupper(x)
}
id_space <- function(x) {
  f <- function(p) mean(grepl(p, x))
  if (f("^ENS[A-Z]*G[0-9]{6,}") > 0.5) return("ensembl")
  if (f("^(NM_|NR_|XM_|XR_)[0-9]+")  > 0.5) return("refseq")
  if (f("^[0-9]+$")                  > 0.5) return("entrez")
  if (f("^[A-Za-z][A-Za-z0-9._-]*$") > 0.5) return("symbol")
  "other"
}
source_path <- function(d) {
  fs <- list.files(file.path("..", "source"), pattern = "\\.rds$", full.names = TRUE)
  hit <- fs[gsub(".*/([^.]*).*", "\\1", fs) == d]
  if (length(hit) != 1L)
    stop("expected exactly one source file for ", d, " under ../source, found ",
         length(hit), if (length(hit)) paste0(": ", paste(basename(hit), collapse = ", ")))
  hit
}

if (task == "universe") {
  suppressPackageStartupMessages(library(Seurat))
  U <- list()
  for (d in unique(base_info$dataset)) {
    message("=== universe: ", d, " ===")
    s  <- readRDS(source_path(d))
    rn <- rownames(s)
    U[[d]] <- list(ids = unique(norm_id(rn)), idspace = id_space(rn),
                   n_genes_raw = length(rn), n_cells = ncol(s))
    if (length(U[[d]]$ids) != length(rn))
      message("  note: ", length(rn) - length(U[[d]]$ids),
              " identifier(s) collapse under normalisation")
    v1 <- unique(base_info$n_genes_total[base_info$dataset == d])
    if (length(v1) == 1L && v1 != length(rn))
      stop(d, ": ", length(rn), " genes now but the v1 export recorded ", v1,
           ". The source object has changed; the v1 sets cannot be compared.")
    cat(sprintf("  %-12s %-8s %6d genes  %5d cells\n", d, U[[d]]$idspace,
                length(rn), ncol(s)))
  }
  saveRDS(U, file.path(out_dir, "universes.rds"))
  cat("written:", file.path(out_dir, "universes.rds"), "\n")
  quit(save = "no", status = 0)
}

RERUN <- c(DIFS = "normal_lookup", DIFS_hartigan = "hartigan")
if (!method %in% names(RERUN))
  stop("method= must be one of ", paste(names(RERUN), collapse = ", "),
       ". The other methods' Baron sets do not pass through SC3 and are reused.")

source("Functions_controlled.R")
if (!file.exists("difs_sc3_complete.R"))
  stop("difs_sc3_complete.R must be in the code directory -- without it this task ",
       "reproduces the defect it exists to repair")
source("difs_sc3_complete.R")
cat("SC3 missing-label repair: ON\n")

t0 <- proc.time()[["elapsed"]]
data.seu     <- readRDS(source_path(dataset))
data.log.mat <- as.matrix(GetAssayData(data.seu, slot = "data"))
data.mat     <- as.matrix(GetAssayData(data.seu, slot = "counts"))

sce <- SingleCellExperiment(assays = list(counts = data.mat,
                                          logcounts = log2(data.mat + 1)))
sce <- sc3_estimate_k(sce)
k_est <- metadata(sce)$sc3$k_estimation
v1_k  <- unique(base_info$k_estimated[base_info$dataset == dataset])
if (length(v1_k) != 1L) stop("v1 export has no single k_estimated for ", dataset)
if (k_est != v1_k)
  stop("sc3_estimate_k gives ", k_est, " on ", dataset, " but the v1 export used ",
       v1_k, ". The two sets would differ by k as well as by the repair.")
cat(sprintf("%s: %d cells, %d genes, k_est = %d (matches v1)\n",
            dataset, ncol(data.log.mat), nrow(data.log.mat), k_est))

res <- difs_features(data.log.mat = data.log.mat, data.seu = data.seu,
                     cluster_count = k_est, binary.bound = log(5),
                     min.expression = log(5),
                     tune_stage1 = TRUE, tune_ratio = TRUE,
                     stage1_key = RERUN[[method]], seed = seed, n_fixed = n_fixed)
genes <- norm_id(res$markers)
if (length(genes) != n_fixed)
  stop(method, " on ", dataset, " returned ", length(genes), " genes, not ", n_fixed,
       ". Recording that as a successful export is what v1 did; this stops instead.")
if (anyDuplicated(genes)) stop("duplicate genes after identifier normalisation")

out <- list(dataset = dataset, method = method, genes = genes,
            k_estimated = k_est, n_fixed = n_fixed, seed = seed,
            stage1_n = res$stage1_n, stage2_n = res$stage2_n, ratio = res$ratio,
            sc3_repair = TRUE, elapsed_sec = proc.time()[["elapsed"]] - t0,
            written = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
            R = R.version.string,
            SC3 = as.character(utils::packageVersion("SC3")))
f <- file.path(out_dir, sprintf("select__%s__%s.rds", dataset, method))
saveRDS(out, f)
cat(sprintf("written: %s  (%d genes, stage I %d, stage II %d, ratio %s, %.0f s)\n",
            f, length(genes), res$stage1_n, res$stage2_n, res$ratio, out$elapsed_sec))
