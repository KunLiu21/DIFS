# difs_sc3fix_compare.R: revision-2 analysis source.
# Source: code/as_run/05_repairs_baron_sc3_hartigan/difs_sc3fix_compare.R
# Usage, scope and limitations: benchmark/REVISION2.md.
difs_args <- function(args, env = parent.frame()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
  invisible(NULL)
}
args <- commandArgs(TRUE)
out <- "../results/diag"; new_set <- "sc3fix"; old_set <- "ncurve"
difs_args(args)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("difs_validate.R")

read_set <- function(set) {
  d  <- file.path("..", "results", set)
  fs <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  if (!length(fs)) stop("no results under ", d)
  g <- function(x, n, def = NA) { v <- x[[n]]; if (is.null(v) || !length(v)) def else v[1] }
  do.call(rbind, lapply(fs, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x)) { cat("*** unreadable, counted as missing: ", f, "\n", sep = ""); return(NULL) }
    data.frame(
      data.name = as.character(g(x, "data.name", "")),
      fmethod   = as.character(g(x, "feature.selection.method", "")),
      clust     = as.character(g(x, "clustering.method", "")),
      n         = as.integer(g(x, "n_requested", NA)),
      k_policy  = as.character(g(x, "k_policy", "")),
      gate_rule = as.character(g(x, "gate_rule", "submitted")),
      gate_rule_src = if (is.null(x[["gate_rule"]]) || !length(x[["gate_rule"]]))
                        "default" else "object",
      gate_n    = as.integer(g(x, "gate_n", NA)),
      seed      = as.integer(g(x, "seed", NA)),
      ARI       = as.numeric(g(x, "ARI", NA)),
      nfeat     = as.integer(g(x, "n_features", NA)),
      total     = as.integer(g(x, "n_cells_total", NA)),
      dropped   = as.integer(g(x, "n_cells_dropped", NA)),
      stringsAsFactors = FALSE)
  }))
}

new <- read_set(new_set); old <- read_set(old_set)
expected <- expected_from_scenarios(new_set)
cat("repaired runs read: ", nrow(new), "   expected by the scenario table: ",
    nrow(expected), "   stored grid runs: ", nrow(old), "\n\n", sep = "")

old_pair <- merge(old, expected, by = DIFS_KEY)
show_gr <- function(df, what) {
  tb <- table(df$gate_rule, df$gate_rule_src)
  cat(what, " gate_rule provenance:\n", sep = ""); print(tb)
  chk <- unique(df[, c("data.name", "gate_n", "gate_rule_src")])
  bad <- do.call(rbind, lapply(split(chk, chk$data.name), function(z)
           if (length(unique(z$gate_n)) > 1L) z else NULL))
  if (!is.null(bad)) {
    cat("*** a dataset reports more than one gate_n across provenances --\n",
        "    the default fill cannot be trusted here:\n", sep = "")
    print(bad, row.names = FALSE)
    quit(save = "no", status = 2)
  }
}
show_gr(new, "new")
if (exists("old_pair")) show_gr(old_pair, "stored")
cat("\n")

difs_gate(new, expected, old = old_pair, require_full_labels = TRUE,
          required_new = c("ARI", "nfeat", "total"),
          required_old = c("ARI", "total"),
          expect_total = 8569L)

m <- merge(new, old_pair, by = DIFS_KEY, suffixes = c("_new", "_old"))
stopifnot(nrow(m) == nrow(new))
m$dARI <- m$ARI_new - m$ARI_old
m <- m[order(-abs(m$dARI)), ]

cat(strrep("=", 76), "\nHOW MUCH DID THE ARI MOVE?\n", strrep("=", 76), "\n", sep = "")
pr <- m[, c("fmethod", "clust", "n", "ARI_old", "ARI_new", "dARI",
            "dropped_old", "dropped_new")]
names(pr) <- c("feature", "clusterer", "n", "ARI_stored", "ARI_repaired",
               "dARI", "drop_stored", "drop_repaired")
print(pr, row.names = FALSE, digits = 4)
md <- stats::median(abs(m$dARI)); mx <- max(abs(m$dARI))
cat(sprintf("\nmedian |dARI| = %.4f   max |dARI| = %.4f  (%s, %s, n = %d)\n",
            md, mx, m$fmethod[1], m$clust[1], m$n[1]))
cat(sprintf("mean dARI = %+.4f   %d of %d repaired runs scored higher\n",
            mean(m$dARI), sum(m$dARI > 0), nrow(m)))

cat("\n", strrep("=", 76), "\nDOES THE ORDERING OF METHODS CHANGE?\n",
    strrep("=", 76), "\n", sep = "")
reordered <- 0L; skipped <- character(0)
for (cl in unique(m$clust)) for (bud in sort(unique(m$n[m$clust == cl]))) {
  s <- m[m$clust == cl & m$n == bud, ]
  if (nrow(s) < 2L) {
    skipped <- c(skipped, sprintf("%s n=%d (%d method)", cl, bud, nrow(s))); next
  }
  r_old <- rank(-s$ARI_old, ties.method = "min")
  r_new <- rank(-s$ARI_new, ties.method = "min")
  same  <- identical(r_old, r_new)
  if (!same) reordered <- reordered + 1L
  ord <- function(r) paste(s$fmethod[order(r)], collapse = " > ")
  cat(sprintf("  %-16s n = %4d  %s\n", cl, bud, if (same) "unchanged" else "CHANGED"))
  cat("      stored  : ", ord(r_old), "\n", sep = "")
  if (!same) cat("      repaired: ", ord(r_new), "\n", sep = "")
}
if (length(skipped))
  cat("\n  not checked (only one feature method, no ordering exists): ",
      paste(skipped, collapse = "; "), "\n",
      "  This check therefore covers the SC3 groupings only.\n", sep = "")

cat("\n", strrep("=", 76), "\nDECISION\n", strrep("=", 76), "\n", sep = "")
cat("Rule fixed in scenario_generator_sc3fix.R before these runs existed.\n\n")
if (mx < 0.02 && reordered == 0L) {
  cat("max |dARI| = ", sprintf("%.4f", mx), " < 0.02, and no ordering changed.\n\n", sep = "")
  cat("-> In these twelve prespecified configurations the repair moves the ARI\n")
  cat("   by less than 0.02 and changes no method ordering. Report the table\n")
  cat("   above as the sensitivity analysis and do not rerun the grid.\n\n")
  cat("   WORDING, and this is not optional: these twelve cells cover k = true,\n")
  cat("   two budgets and one seed on ONE dataset. They do NOT show that the\n")
  cat("   other stored Baron results are correct, and the repaired runs remain\n")
  cat("   the only ones that scored all 8569 cells. Write \"the difference is\n")
  cat("   small in these prespecified configurations\", never \"below the\n")
  cat("   resolution of anything the paper claims\".\n")
} else if (mx < 0.02 && reordered > 0L) {
  cat("max |dARI| = ", sprintf("%.4f", mx), " < 0.02 BUT ", reordered,
      " grouping(s) reordered.\n\n", sep = "")
  cat("-> The magnitude threshold passes and the thing it was a proxy for does\n")
  cat("   not. A shift too small to notice that changes which method wins is\n")
  cat("   exactly the case the threshold was supposed to catch. Treat this as\n")
  cat("   the 0.02-0.05 branch: rerun the full Baron x SC3 set at k = true.\n")
} else if (mx < 0.05) {
  cat("max |dARI| = ", sprintf("%.4f", mx), " in [0.02, 0.05)\n\n", sep = "")
  cat("-> Report it, AND rerun the full Baron x SC3 set at k = true (45 cells).\n")
} else {
  cat("max |dARI| = ", sprintf("%.4f", mx), " >= 0.05\n\n", sep = "")
  cat("-> Every stored Baron number is provisional. Rerun before using any.\n")
}

cat("\n", strrep("=", 76), "\nIF A RERUN IS CALLED FOR, THE AFFECTED SET IS\n",
    strrep("=", 76), "\n", sep = "")
cat("  90  main grid  Baron x SC3-as-final-clusterer (both k policies)\n")
cat("  18  main grid  Baron x DIFS x Refined Louvain -- no cells dropped in the\n")
cat("      final scoring, but stage II's preliminary clustering is the same SC3\n")
cat("      call, so the FEATURE SET was chosen from 5000 of 8569 cells\n")
cat("   3  ablation   Baron x DIFS x {100, 300, 1000} -- same reason as above.\n")
cat("      These are easy to miss because the ablation is Louvain-only.\n")
cat("  --\n")
cat(" 111  cells total. Which of them actually need rerunning depends on which\n")
cat("      results the final manuscript keeps; decide that first.\n")

cat("\n", strrep("=", 76), "\nFOR THE WRITE-UP\n", strrep("=", 76), "\n", sep = "")
cat("Dropping Baron from the 13-dataset mean rank moves DIFS from 1.692 to\n")
cat("1.615 (k = true). That single sensitivity does not by itself license the\n")
cat("sentence \"every conclusion is unchanged\" -- recompute each claim you want\n")
cat("to make that statement about, including the component analysis, which has\n")
cat("not been recomputed without Baron.\n")

verdict <- if (mx >= 0.05) {
  "rerun_required_large_shift"
} else if (reordered > 0L) {
  "rerun_required_reordering"
} else if (mx >= 0.02) {
  "rerun_required_moderate_shift"
} else {
  "small_in_checked_configs"
}
cat("\nVERDICT (machine-readable): ", verdict, "\n", sep = "")
utils::write.csv(data.frame(verdict = verdict, max_abs_dARI = mx,
                            median_abs_dARI = md, mean_dARI = mean(m$dARI),
                            groupings_reordered = reordered, n_cells = nrow(m),
                            stringsAsFactors = FALSE),
                 file.path(out, "sc3fix_verdict.csv"), row.names = FALSE)
utils::write.csv(pr, file.path(out, "sc3fix_comparison.csv"), row.names = FALSE)
cat("\nwritten: ", normalizePath(out), "/sc3fix_comparison.csv\n", sep = "")
