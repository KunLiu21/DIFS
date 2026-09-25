
difs_args <- function(args, known, env = parent.frame()) {
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
args <- commandArgs(TRUE)
grid <- "../results/diag/grid_corrected.csv"
abl  <- "../results/diag/abl_corrected.csv"
out  <- "../results/supp"
difs_args(args, c("grid", "abl", "out"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("difs_validate.R")

read_runs <- function(set) {
  d <- file.path("..", "results", set)
  fs <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  if (!length(fs)) return(NULL)
  g <- function(x, n, def = NA) { v <- x[[n]]; if (is.null(v) || !length(v)) def else v[1] }
  do.call(rbind, lapply(fs, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(x)) return(NULL)
    data.frame(
      result_set = set,
      dataset = as.character(g(x, "data.name", "")),
      feature_method = as.character(g(x, "feature.selection.method", "")),
      clustering_method = as.character(g(x, "clustering.method", "")),
      k_policy = as.character(g(x, "k_policy", "")),
      gate_rule = as.character(g(x, "gate_rule", "submitted")),
      gate_rule_recorded = !is.null(x[["gate_rule"]]) && length(x[["gate_rule"]]) > 0,
      seed = as.integer(g(x, "seed", NA)),
      n_requested = as.integer(g(x, "n_requested", NA)),
      n_features_realised = as.integer(g(x, "n_features", NA)),
      gate_n = as.integer(g(x, "gate_n", NA)),
      k_true = as.integer(g(x, "k_true", NA)),
      k_estimated = as.integer(g(x, "k_estimated", NA)),
      k_used_clustering = as.integer(g(x, "k_used_clustering", NA)),
      n_clusters_found = as.integer(g(x, "n_clusters_found", NA)),
      n_cells_total = as.integer(g(x, "n_cells_total", NA)),
      n_cells_scored = as.integer(g(x, "n_cells_total", NA)) -
                       as.integer(g(x, "n_cells_dropped", 0)),
      n_cells_dropped = as.integer(g(x, "n_cells_dropped", NA)),
      n_cells_unlabelled = as.integer(g(x, "n_cells_unlabelled", NA)),
      difs_stage1_n = as.integer(g(x, "difs_stage1_n", NA)),
      difs_stage2_n = as.integer(g(x, "difs_stage2_n", NA)),
      difs_ratio_nominal = as.character(g(x, "difs_ratio", NA)),
      ARI = as.numeric(g(x, "ARI", NA)),
      FM = as.numeric(g(x, "FM", NA)),
      Jaccard = as.numeric(g(x, "Jaccard", NA)),
      Purity = as.numeric(g(x, "Purity", NA)),
      elapsed_sec = as.numeric(g(x, "elapsed_sec", NA)),
      source_file = basename(f), stringsAsFactors = FALSE)
  }))
}

R <- do.call(rbind, Filter(Negate(is.null),
       lapply(c("ncurve", "abl", "hart", "baronfix"), read_runs)))
cat("per-run rows read: ", nrow(R), "\n", sep = "")
print(table(R$result_set))

BARON_KEY <- c("dataset","feature_method","clustering_method","n_requested",
               "k_policy","gate_rule","seed")
METRICS   <- c("ARI","FM","Jaccard","Purity")
CARRY     <- c(METRICS, "n_features_realised", "n_clusters_found",
               "n_cells_total", "n_cells_scored", "n_cells_dropped",
               "difs_stage1_n", "difs_stage2_n", "difs_ratio_nominal", "elapsed_sec",
               "n_cells_unlabelled", "k_true", "k_estimated", "k_used_clustering", "gate_n")

key <- function(d, cols) do.call(paste, c(lapply(cols, function(k) as.character(d[[k]])), sep = "|"))

fix <- R[R$result_set == "baronfix", ]
expected_fix <- expected_from_scenarios("baronfix")
checked_fix <- data.frame(data.name=fix$dataset, fmethod=fix$feature_method,
  clust=fix$clustering_method, n=fix$n_requested, k_policy=fix$k_policy,
  gate_rule=fix$gate_rule, seed=fix$seed, ARI=fix$ARI, FM=fix$FM,
  Jaccard=fix$Jaccard, Purity=fix$Purity, nfeat=fix$n_features_realised,
  total=fix$n_cells_total, dropped=fix$n_cells_dropped,
  unlabelled=fix$n_cells_unlabelled, stringsAsFactors=FALSE)
difs_gate(checked_fix, expected_fix, require_full_labels=TRUE,
  required_new=c(METRICS,"nfeat","total","unlabelled"), expect_total=8569L)
if (any(checked_fix$unlabelled != 0L)) stop("repaired Baron results contain unlabelled cells")
if (anyDuplicated(paste(R$result_set, key(R,BARON_KEY))))
  stop("duplicate configurations within a result set")
tgt <- R$result_set %in% c("ncurve", "abl", "hart") & R$dataset == "Baron"
kf  <- key(fix, BARON_KEY)
if (anyDuplicated(kf)) stop("baronfix runs are not unique on ", paste(BARON_KEY, collapse="+"))
j   <- match(key(R, BARON_KEY), kf)
hit <- tgt & !is.na(j)
required_hit <- R$dataset == "Baron" & (
  (R$result_set == "ncurve" & (R$clustering_method == "SC3" | R$feature_method == "DIFS")) |
  (R$result_set == "abl" & R$feature_method == "DIFS") |
  R$result_set == "hart")
if (!identical(hit, required_hit)) stop("the replacement keys do not cover exactly the affected analysis rows")

old_ARI <- R$ARI
for (cc in CARRY) R[[cc]][hit] <- fix[[cc]][j[hit]]
R$baron_repaired <- hit
R$ARI_before_repair <- ifelse(hit, old_ARI, NA_real_)
R$source_result_set <- ifelse(hit, "baronfix", R$result_set)

cat("Baron repair substituted on ", sum(hit), " rows (",
    sum(hit & R$result_set == "ncurve"), " main grid, ",
    sum(hit & R$result_set == "abl"), " ablation, ",
    sum(hit & R$result_set == "hart"), " reference-distribution); all four metrics carried.\n",
    sep = "")
moved <- hit & abs(R$ARI - old_ARI) > 1e-12
cat("  of these, ", sum(moved), " actually changed ARI; range ",
    sprintf("%+.4f .. %+.4f", min((R$ARI - old_ARI)[moved]),
                              max((R$ARI - old_ARI)[moved])), "\n", sep = "")
if (!nrow(fix))
  stop("no baronfix runs under ../results/baronfix. Without them this table would\n",
       "  carry the DEFECTIVE Baron values, so it is not written at all.")
if (sum(hit) == 0L)
  stop("the baronfix runs matched no ncurve/abl row. Check ", paste(BARON_KEY, collapse="+"), ".")
if (!all(is.finite(as.matrix(R[hit, METRICS, drop=FALSE]))))
  stop("a substituted row has a non-finite metric")
if (any(hit & (R$n_cells_dropped != 0 | R$n_cells_total != 8569)))
  stop("a repaired row does not cover all 8569 cells. The repaired runs all do,\n",
       "  so this means the substitution did not carry the coverage fields.")

xcheck <- function(path, set) {
  if (!file.exists(path)) {
    stop("required corrected table not found for ", set, ": ", path)
  }
  cg <- utils::read.csv(path, stringsAsFactors = FALSE)
  cg <- cg[cg$data.name == "Baron", ]
  names(cg)[match(c("data.name","fmethod","clust","n"), names(cg))] <-
    c("dataset","feature_method","clustering_method","n_requested")
  a <- R[R$result_set == set & R$dataset == "Baron", ]
  m <- match(key(a, BARON_KEY), key(cg, BARON_KEY))
  if (anyNA(m)) stop(sum(is.na(m)), " Baron ", set, " row(s) absent from ", path)
  pairs <- list(ARI = "ARI", FM = "FM", Jaccard = "Jaccard", Purity = "Purity",
                n_cells_dropped = "dropped", n_cells_total = "total",
                n_features_realised = "nfeat", difs_stage2_n = "stage2_n")
  have <- pairs[vapply(pairs, function(z) z %in% names(cg), logical(1))]
  worst <- 0; checked <- character()
  for (nm in names(have)) {
    x <- suppressWarnings(as.numeric(a[[nm]])); y <- suppressWarnings(as.numeric(cg[[have[[nm]]]][m]))
    if (any(is.na(x) != is.na(y)) || any(!is.na(x) & !is.finite(x)) ||
        any(!is.na(y) & !is.finite(y))) stop("invalid or mismatched values in cross-check field ", nm)
    ok <- is.finite(x) & is.finite(y)
    if (!any(ok)) next
    worst <- max(worst, max(abs(x[ok] - y[ok]))); checked <- c(checked, nm)
  }
  cat(sprintf("  cross-check %-6s vs %-22s %3d rows, %d field(s) [%s], max |diff| = %.2e %s\n",
              set, basename(path), nrow(a), length(checked),
              paste(checked, collapse = ","), worst, if (worst < 1e-9) "OK" else "*** MISMATCH"))
  if (length(checked) < 4L)
    cat("    NOTE: that corrected table carries only ", length(checked),
        " comparable field(s). Re-run difs_baronfix_merge.R so it writes the whole\n",
        "    record, then this check covers all four metrics and the coverage fields.\n", sep = "")
  if (worst >= 1e-9) stop("this script and the merge script disagree on Baron -- stop")
}
xcheck(grid, "ncurve")
xcheck(abl,  "abl")

R$is_analysis_row <- R$result_set != "baronfix"
cat("S0: ", sum(R$is_analysis_row), " analysis row(s) + ", sum(!R$is_analysis_row),
    " provenance row(s) (the repaired runs themselves). Filter on is_analysis_row\n",
    "    before summarising, or Baron is counted twice.\n", sep = "")
utils::write.csv(R, file.path(out, "S0_per_run.csv"), row.names = FALSE)

main <- R[R$result_set == "ncurve", ]
for (i in seq_along(c("ARI","FM","Jaccard","Purity"))) {
  m <- c("ARI","FM","Jaccard","Purity")[i]
  s <- main[, c("dataset","feature_method","clustering_method","k_policy",
                "n_requested","n_features_realised", m, "baron_repaired")]
  names(s)[ncol(s) - 1L] <- "value"
  s <- s[, c(setdiff(names(s), "baron_repaired"), "baron_repaired")]
  utils::write.csv(s[order(s$dataset, s$k_policy, s$n_requested,
                           s$clustering_method, s$feature_method), ],
                   file.path(out, sprintf("S%d_%s_per_run.csv", i, m)), row.names = FALSE)
}

exp_grid <- expand.grid(dataset = unique(main$dataset),
                        feature_method = unique(main$feature_method),
                        clustering_method = unique(main$clustering_method),
                        k_policy = unique(main$k_policy),
                        n_requested = unique(main$n_requested),
                        stringsAsFactors = FALSE)
kk <- function(d) do.call(paste, c(d[, c("dataset","feature_method","clustering_method",
                                         "k_policy","n_requested")], sep="|"))
missing <- exp_grid[!kk(exp_grid) %in% kk(main), ]
utils::write.csv(missing, file.path(out, "S5_missing_configurations.csv"), row.names = FALSE)
cat("\nexpected main-grid configurations: ", nrow(exp_grid),
    "   present: ", nrow(main), "   missing: ", nrow(missing), "\n", sep = "")
if (nrow(missing)) print(table(missing$dataset, missing$feature_method))

cat("\nwritten to ", normalizePath(out), ":\n",
    "  S0_per_run.csv            every run of the ncurve/abl/hart/baronfix sets,\n",
    "                            every recorded field. It does NOT include the\n",
    "                            Array A controlled comparison, which is a separate\n",
    "                            experiment under ../results/controlled_comparison.\n",
    "  S1..S4_<metric>_per_run.csv  one file per metric, with baron_repaired\n",
    "  S5_missing_configurations.csv  which configurations are absent. It lists the\n",
    "                            coordinates only; the REASON each failed is in the\n",
    "                            manuscript text and is not recoverable from the\n",
    "                            stored objects, because a run that failed wrote none.\n", sep = "")
cat("\nNOT produced, because the pipeline never saved it: the returned gene set,\n",
    "the realised stage I / stage II split (only the two list sizes and the\n",
    "nominal ratio were stored), and per-cell-type recall.\n", sep = "")
