args <- commandArgs(TRUE)
for (a in args) {
  if (!startsWith(a, "out_dir=")) stop("Only out_dir=... is accepted")
  out_dir <- substring(a, 9L)
}
if (!exists("out_dir")) out_dir <- "../results/simulation"
rep_dir <- file.path(out_dir, "report")
dir.create(rep_dir, recursive = TRUE, showWarnings = FALSE)
rd <- function(f) { p <- file.path(out_dir, f)
  if (file.exists(p)) readRDS(p) else { message("missing: ", p); NULL } }
wr <- function(x, f) if (!is.null(x)) {
  utils::write.csv(x, file.path(rep_dir, f), row.names = FALSE)
  cat("\n----", f, "----\n"); print(utils::head(x, 40), row.names = FALSE) }

inv <- rd("E1_invariance.rds"); mono <- rd("E1_monotonicity.rds")
if (!is.null(inv)) {
  t1 <- data.frame(quantity = c("max |dip(x) - dip(a+bx)| over all replicates"),
                   value = max(inv$max_abs_diff))
  if (!is.null(mono)) t1 <- rbind(t1, data.frame(
    quantity = c(sprintf("Spearman(D, p_normal) at n=%d", mono$n),
                 sprintf("strictly discordant pairs at n=%d", mono$n),
                 sprintf("tied Monte-Carlo p-values (of 200) at n=%d", mono$n)),
    value = c(mono$spearman_D_vs_pnormal_lookup,
              mono$n_strictly_discordant_pairs, mono$n_tied_pvalues_mc)))
  wr(t1, "T1_invariance_monotonicity.csv")
}

e2 <- rd("E2_null_calibration.rds")
if (!is.null(e2)) {
  a <- subset(e2, alpha == 0.05)
  t2 <- a[, c("n_cells", "mu", "null_type", "n_tested", "median_n_expressing",
              "median_distinct_frac", "fpr_modified", "fpr_original",
              "ks_modified", "ks_original")]
  t2 <- t2[order(t2$n_cells, t2$mu, t2$null_type), ]
  wr(t2, "T2_type1_error.csv")
}

e3 <- rd("E3_power.rds")
if (!is.null(e3)) wr(e3[order(e3$n_cells, e3$mu, e3$pi, e3$log2fc),
  c("n_cells", "mu", "pi", "log2fc", "frac_testable",
    "power_nominal_modified", "power_nominal_original",
    "power_matched_modified", "power_matched_original")], "T3_power.csv")

e4 <- rd("E4_ranking.rds")
if (!is.null(e4)) {
  ag <- stats::aggregate(cbind(auc, precision_500, n_for_recall50, n_for_recall80) ~
                           method + log2fc + disp + n_cells, e4, mean)
  wr(ag[order(ag$n_cells, ag$disp, ag$log2fc, -ag$auc), ], "T4_ranking.csv")
  w <- stats::reshape(stats::aggregate(auc ~ method + log2fc, e4, mean),
                      idvar = "method", timevar = "log2fc", direction = "wide")
  utils::write.csv(w, file.path(rep_dir, "T4_ranking_auc_wide.csv"), row.names = FALSE)
  cat("\n---- AUC by effect size (rows = method, cols = log2 fold change) ----\n")
  print(w, row.names = FALSE, digits = 3)
}

wr(rd("E5_ties_runtime.rds"),      "T5_ties_runtime.csv")
e6 <- rd("E6_gate_reference.rds")
if (!is.null(e6)) {
  e6$gate_label <- ifelse(abs(e6$min_expression - log(5)) < 1e-8,
                          "log(5)  <- as implemented",
                          ifelse(abs(e6$min_expression - 0.5) < 1e-8,
                                 "0.5  <- as documented", ""))
  wr(e6, "T6_gate_sensitivity.csv")
}
e7 <- rd("E7_stage2_robustness.rds")
if (!is.null(e7)) wr(stats::aggregate(
  cbind(ari_vs_truth, n_selected, precision, recall, enrichment) ~ scenario, e7, mean),
  "T7_stage2_robustness.csv")

cat("\nreport written to ", normalizePath(rep_dir), "\n")
