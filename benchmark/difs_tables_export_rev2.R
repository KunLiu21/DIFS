# difs_tables_export_rev2.R: revision-2 analysis source.
# Source: code/as_run/10_tables_figures/difs_tables_export_rev2.R
# Usage, scope and limitations: benchmark/REVISION2.md.
SETS <- c("abl2", "reps", "krange")
out <- "../results/supp"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

g <- function(x, n, def = NA) { v <- x[[n]]; if (is.null(v) || !length(v)) def else v[1] }
read_set <- function(set) {
  fs <- list.files(file.path("..", "results", set), pattern = "\\.rds$", full.names = TRUE)
  if (!length(fs)) stop("no run files in ../results/", set)
  do.call(rbind, lapply(fs, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x)) stop("unreadable run file: ", f)
    data.frame(
      result_set = set,
      status = as.character(g(x, "status", "ok")),
      dataset = as.character(g(x, "data.name", "")),
      arm = as.character(g(x, "arm", "")),
      runner_method = as.character(g(x, "runner_fmethod", NA)),
      stage1_source = as.character(g(x, "stage1_source", NA)),
      final_set = as.character(g(x, "final_set", NA)),
      clustering_method = as.character(g(x, "clustering.method", "")),
      k_policy = as.character(g(x, "k_policy", "")),
      k_offset = as.integer(g(x, "k_offset", NA)),
      cluster_seed = as.integer(g(x, "cluster_seed", NA)),
      subsample_rep = as.integer(g(x, "subsample_rep", 0L)),
      subsample_frac = as.numeric(g(x, "subsample_frac", 1)),
      seed = as.integer(g(x, "seed", NA)),
      n_requested = as.integer(g(x, "n_requested", NA)),
      n_features_realised = as.integer(g(x, "n_features", NA)),
      gate_n = as.integer(g(x, "gate_n", NA)),
      k_true = as.integer(g(x, "k_true", g(x, "k_true_labels", NA))),
      k_used_clustering = as.integer(g(x, "k_used_clustering", NA)),
      n_clusters_found = as.integer(g(x, "n_clusters_found", NA)),
      n_cells_used = as.integer(g(x, "n_cells_used", NA)),
      n_cells_dropped = as.integer(g(x, "n_cells_dropped", NA)),
      difs_stage1_n = as.integer(g(x, "difs_stage1_n", NA)),
      difs_stage2_n = as.integer(g(x, "difs_stage2_n", NA)),
      difs_ratio_nominal = as.character(g(x, "difs_ratio", NA)),
      stage2only_fell_back = as.logical(g(x, "stage2only_fell_back", NA)),
      ARI = as.numeric(g(x, "ARI", NA)), FM = as.numeric(g(x, "FM", NA)),
      Jaccard = as.numeric(g(x, "Jaccard", NA)), Purity = as.numeric(g(x, "Purity", NA)),
      elapsed_sec = as.numeric(g(x, "elapsed_sec", NA)),
      source_file = basename(f), stringsAsFactors = FALSE)
  }))
}

R <- do.call(rbind, lapply(SETS, read_set))
key <- function(d) paste(d$result_set, d$dataset, d$arm, d$clustering_method, d$n_requested,
                         d$k_offset, d$cluster_seed, d$subsample_rep, d$seed, sep = "|")
if (anyDuplicated(key(R))) stop(sum(duplicated(key(R))), " duplicated run(s); resolve before export")

E <- do.call(rbind, lapply(SETS, function(set) {
  e <- new.env(); load(sprintf("../scenarios/scenarios_%s.Rdata", set), envir = e); sc <- e$scenarios
  data.frame(result_set = set, task = seq_len(nrow(sc)),
             dataset = gsub(".*/([^.]*).*", "\\1", sc$data.path), arm = sc$arm,
             clustering_method = sc$method, n_requested = as.integer(sc$n_features),
             k_offset = as.integer(sc$k_offset), cluster_seed = as.integer(sc$cluster_seed),
             subsample_rep = as.integer(sc$subsample_rep), seed = as.integer(sc$seed),
             stringsAsFactors = FALSE)
}))
miss <- E[!key(E) %in% key(R), ]
extra <- R[!key(R) %in% key(E), ]
if (nrow(extra)) stop(nrow(extra), " run file(s) match no scenario row -- stale files in a results directory?")

utils::write.csv(R[order(R$result_set, R$dataset, R$arm, R$clustering_method, R$n_requested,
                         R$k_offset, R$cluster_seed, R$subsample_rep), ],
                 file.path(out, "S0b_rev2_per_run.csv"), row.names = FALSE)
utils::write.csv(miss, file.path(out, "S5b_rev2_missing.csv"), row.names = FALSE)

cat("second-round runs written: ", nrow(R), "\n", sep = "")
print(table(R$result_set, R$status))
cat("\nexpected tasks: ", nrow(E), "   present: ", nrow(E) - nrow(miss),
    "   missing: ", nrow(miss), "\n", sep = "")
if (nrow(miss)) print(miss, row.names = FALSE)
cat("\nwritten to ", normalizePath(out), ": S0b_rev2_per_run.csv, S5b_rev2_missing.csv\n", sep = "")
