# difs_fig_dip_null.R: revision-2 analysis source.
# Source: code/as_run/10_tables_figures/difs_fig_dip_null.R
# Usage, scope and limitations: benchmark/REVISION2.md.
suppressPackageStartupMessages(library(diptest))
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1])))
source(file.path(here, "difs_sim_core.R"))                     # dip_pvalue_lookup_vec()
args <- commandArgs(TRUE)
tf <- if (length(args)) args[1] else file.path(here, "../inst/extdata/dip_null_tables.rds")
TABS <- readRDS(tf)
set.seed(20261005)
n_samples <- 1000; n <- 1000
res <- t(vapply(seq_len(n_samples), function(i) {
  x <- stats::rnorm(n, mean = stats::runif(1, 0, 5), sd = stats::runif(1, 0.5, 2))
  D <- diptest::dip(sort(x))
  c(D = D, p_original = suppressWarnings(diptest::dip.test(x)$p.value),
    p_modified = dip_pvalue_lookup_vec(D, n, TABS$normal))
}, numeric(3)))
utils::write.csv(data.frame(sample = seq_len(n_samples), res), "fig_dip_null_pvalues.csv", row.names = FALSE)
grDevices::pdf("modified_dip.pdf", width = 9, height = 3.6)
op <- graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.5, 2.6, 1), cex.main = 1.05)
graphics::hist(res[, "p_original"], breaks = seq(0, 1, 0.1), col = "skyblue", main = "Original dip test (uniform reference)",
               xlab = "p-value", ylab = "Number of samples")
graphics::hist(res[, "p_modified"], breaks = seq(0, 1, 0.1), col = "lightcoral", main = "Modified dip test (Gaussian reference)",
               xlab = "p-value", ylab = "Number of samples")
graphics::par(op); grDevices::dev.off()
cat(sprintf("original: share of p > 0.9 = %.3f, min p = %.3f\nmodified: KS test against U(0,1) p = %.3f; share p < 0.05 = %.3f\n",
            mean(res[, "p_original"] > 0.9), min(res[, "p_original"]),
            suppressWarnings(stats::ks.test(res[, "p_modified"], "punif")$p.value), mean(res[, "p_modified"] < 0.05)))
