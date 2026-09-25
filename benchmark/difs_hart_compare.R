
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
out <- "../results/diag"; new_set <- "hart"; old_set <- "ncurve"
corrected <- ""
difs_args(args, c("out", "new_set", "old_set", "corrected"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("difs_validate.R")
source("difs_signrank.R")

NOT_INDEPENDENT <- c("SimKumar4easy", "SimKumar4hard", "Zhengmix4eq", "Zhengmix8eq")

read_set <- function(set, keep = NULL) {
  d  <- file.path("..", "results", set)
  fs <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  if (!length(fs)) stop("no results under ", d)
  g <- function(x, n, def = NA) { v <- x[[n]]; if (is.null(v) || !length(v)) def else v[1] }
  z <- do.call(rbind, lapply(fs, function(f) {
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
  if (!is.null(keep)) z <- z[z$fmethod %in% keep, , drop = FALSE]
  z
}

new <- read_set(new_set)
old <- read_set(old_set, keep = "DIFS")
expected <- expected_from_scenarios(new_set)
pair_scope <- setdiff(DIFS_KEY, "fmethod")
old <- old[difs_key_string(old, pair_scope) %in%
             difs_key_string(expected, pair_scope), , drop = FALSE]

SUBKEY <- c("data.name", "clust", "n", "k_policy", "gate_rule", "seed")
BOTH   <- c("ARI", "nfeat", "total", "dropped")
if (nzchar(corrected)) {
  if (!file.exists(corrected)) stop("corrected table not found: ", corrected)
  cg <- utils::read.csv(corrected, stringsAsFactors = FALSE)

  wide_cols <- c("ARI_d", "ARI_h")
  if (!all(wide_cols %in% names(cg)))
    stop("corrected= expects the PAIRED table written by difs_baronfix_merge.R\n",
         "  (baron_hart_repaired.csv: one row per Baron cell, columns ARI_d and ARI_h).\n",
         "  The file given has columns: ", paste(names(cg), collapse = ", "), "\n",
         "  grid_corrected.csv carries only the main grid and has no DIFS_hartigan arm,\n",
         "  so passing it would repair one arm and leave the other defective.")
  if (!all(SUBKEY %in% names(cg))) stop("corrected table lacks configuration keys")
  expected_baron <- expected[expected$data.name == "Baron", SUBKEY, drop = FALSE]
  if (anyDuplicated(difs_key_string(cg, SUBKEY)) ||
      !setequal(difs_key_string(cg, SUBKEY), difs_key_string(expected_baron, SUBKEY)))
    stop("corrected table must contain each expected Baron pair exactly once, without extras")

  sub_in <- function(df, src, suffix, label) {
    j <- match(difs_key_string(df, SUBKEY), difs_key_string(src, SUBKEY))
    hit <- !is.na(j)
    n_baron <- sum(df$data.name == "Baron")
    if (sum(hit) != n_baron)
      stop(label, ": the corrected table covers ", sum(hit), " of this arm's ",
           n_baron, " Baron row(s). A partial substitution is not a substitution.")
    if (anyDuplicated(j[hit])) stop(label, ": a corrected row matches more than one run")
    for (cc in BOTH) {
      sc <- paste0(cc, suffix)
      if (!sc %in% names(src)) stop(label, ": corrected table lacks column ", sc)
      df[[cc]][hit] <- src[[sc]][j[hit]]
    }
    if (!all(is.finite(as.matrix(df[hit, BOTH, drop = FALSE]))))
      stop(label, ": a substituted row has non-finite metrics or coverage")
    if (any(df$dropped[hit] != 0L) || any(df$total[hit] != 8569L))
      stop(label, ": a substituted Baron row does not cover all 8569 cells")
    cat(sprintf("  %-28s %d row(s) replaced, whole record, all covering 8569 cells\n",
                label, sum(hit)))
    df
  }
  cat("\nSubstituting repaired Baron runs on BOTH arms from ", basename(corrected), ":\n", sep = "")
  old <- sub_in(old, cg, "_d", "DIFS (stored arm)")
  new <- sub_in(new, cg, "_h", "DIFS_hartigan arm")
  cat("\n")
} else {
  cat("\n*** corrected= NOT GIVEN. Both arms of this comparison carry the SC3\n",
      "    missing-label defect on Baron. The comparison remains internally valid,\n",
      "    but its Baron cells swung by -0.092 and +0.138 between the two clustering\n",
      "    methods and the 100-feature result depended on them. DO NOT QUOTE THIS RUN.\n\n",
      sep = "")
}
expected <- expected_from_scenarios(new_set)
cat("hartigan runs read: ", nrow(new), "   expected by the scenario table: ",
    nrow(expected), "   stored DIFS runs available: ", nrow(old), "\n\n", sep = "")

old_pair <- merge(old, expected[, setdiff(DIFS_KEY, "fmethod")],
                  by = setdiff(DIFS_KEY, "fmethod"))
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

difs_gate(new, expected, old = old_pair, require_full_labels = FALSE,
          pair_key     = setdiff(DIFS_KEY, "fmethod"),
          required_new = c("ARI", "nfeat", "total"),
          required_old = c("ARI", "total"))

m <- merge(new, old_pair, by = setdiff(DIFS_KEY, "fmethod"), suffixes = c("_h", "_d"))
stopifnot(nrow(m) == nrow(new))
m$dARI <- m$ARI_d - m$ARI_h          # positive = DIFS (Gaussian) ahead

coverage <- m[m$data.name == "Baron" | m$total_d > 5000 | m$total_h > 5000, ]
if (nrow(coverage)) {
  cc <- c("total_d", "total_h", "dropped_d", "dropped_h")
  if (!all(is.finite(as.matrix(coverage[, cc, drop=FALSE]))) ||
      any(coverage$dropped_d != 0L | coverage$dropped_h != 0L) ||
      any(coverage$data.name == "Baron" &
          (coverage$total_d != 8569L | coverage$total_h != 8569L)))
    stop("paired label coverage is incomplete; pass corrected=baron_hart_repaired.csv")
}

CLUSTERERS <- sort(unique(expected$clust))
cnt <- as.data.frame(table(data.name = m$data.name, n = m$n), stringsAsFactors = FALSE)
cnt$n <- as.integer(as.character(cnt$n))
short <- cnt[cnt$Freq != length(CLUSTERERS), ]
if (nrow(short)) {
  cat("\n*** these dataset x budget cells are not backed by all ",
      length(CLUSTERERS), " clustering methods:\n", sep = "")
  print(short, row.names = FALSE)
  cat("Averaging them would compare a one-clusterer mean against two-clusterer\n",
      "means elsewhere. Stopping.\n", sep = "")
  quit(save = "no", status = 2)
}
agg <- aggregate(m[, c("ARI_d", "ARI_h")],
                 by = list(data.name = m$data.name, n = m$n), FUN = mean)
stopifnot(!any(is.na(agg$ARI_d)), !any(is.na(agg$ARI_h)))
agg$dARI <- agg$ARI_d - agg$ARI_h

report <- function(a, label) {
  d <- a$dARI
  w  <- difs_signrank(a$ARI_d, a$ARI_h)
  tt <- stats::t.test(a$ARI_d, a$ARI_h, paired = TRUE)
  cat(sprintf("%-34s %2d datasets   mean %+.4f   DIFS ahead %d/%d   exact p = %.4f\n",
              label, nrow(a), mean(d), sum(d > 0), nrow(a), w$p.value))
  cat(sprintf("%-34s paired t p = %.4f   95%% CI [%+.4f, %+.4f]   (%s)\n",
              "", tt$p.value, tt$conf.int[1], tt$conf.int[2], w$method))
  list(p = w$p.value, mean = mean(d), n = nrow(a),
       lo = tt$conf.int[1], hi = tt$conf.int[2])
}

cat(strrep("=", 76), "\nRESULTS  (9 independent = primary inference; 13 = full panel alongside; 300 = secondary)\n",
    strrep("=", 76), "\n", sep = "")
res <- list()
for (bud in sort(unique(agg$n))) {
  cat("\n-- ", bud, " features --\n", sep = "")
  a13 <- agg[agg$n == bud, ]
  a9  <- a13[!a13$data.name %in% NOT_INDEPENDENT, ]
  a8  <- a9[a9$data.name != "Baron", ]
  res[[paste0("all13_", bud)]] <- report(a13, "all 13 (full panel)")
  res[[paste0("ind9_",  bud)]] <- report(a9,  "9 independent real (PRIMARY)")
  res[[paste0("ind8_",  bud)]] <- report(a8,  "  same, Baron removed")
}

cat("\n", strrep("=", 76), "\nPER DATASET\n", strrep("=", 76), "\n", sep = "")
for (bud in sort(unique(agg$n))) {
  a <- agg[agg$n == bud, ]
  a$independent <- ifelse(a$data.name %in% NOT_INDEPENDENT, "", "*")
  cat("\n-- ", bud, " features  (* = mutually independent real dataset) --\n", sep = "")
  print(a[order(-a$dARI), c("data.name", "independent", "ARI_d", "ARI_h", "dARI")],
        row.names = FALSE, digits = 4)
}

cat("\n", strrep("=", 76), "\nBARON, REPORTED SEPARATELY\n", strrep("=", 76), "\n", sep = "")
if (nzchar(corrected)) {
  cat("Both Baron arms use the repaired results and cover 8569 cells.\n",
      "The Baron-excluded sensitivity is reported separately above.\n\n", sep = "")
} else {
  cat("Both Baron arms retain the historical missing-label defect.\n",
      "This is a historical comparison, not a corrected result.\n\n", sep = "")
}
print(m[m$data.name == "Baron", c("clust", "n", "ARI_d", "ARI_h", "dARI")],
      row.names = FALSE, digits = 4)

cat("\n", strrep("=", 76), "\nDECISION\n", strrep("=", 76), "\n", sep = "")
cat("Rule fixed in scenario_generator_hart.R before these runs existed and\n",
    "amended in PREREG-hart-amendment-20260919.md before any result was read.\n",
    "Branch is read from the NINE independent real datasets at 100 features.\n\n", sep = "")
r <- res[["ind9_100"]]
cat(sprintf("  9 independent, n = 100:  mean %+.4f   95%% CI [%+.4f, %+.4f]   p = %.4f\n\n",
            r$mean, r$lo, r$hi, r$p))
if (r$p >= 0.05 && abs(r$mean) < 0.02) {
  cat("-> No stable difference observed under the conditions tested.\n")
  cat("   Section 2.2 stands as written. Report this together with the Array A\n")
  cat("   1000-feature result as the answer to Reviewer 1's first major comment.\n")
  cat("   WORDING: \"no stable difference was observed under the conditions\n")
  cat("   tested\" -- NOT \"the two are equivalent\". A null result at n = 9 does\n")
  cat("   not establish equivalence; quote the interval above.\n")
} else if (r$p < 0.05 && r$mean > 0) {
  cat("-> The Gaussian reference is ahead on these data.\n")
  cat("   Section 2.2 may say so, citing THIS experiment, and must also report\n")
  cat("   the 1000-feature null result.\n")
  cat("   WORDING: the branch was read at 100 features ONLY. It licenses a\n")
  cat("   claim at 100 features, true k, seed 1, on these datasets. The\n")
  cat("   300-feature result stands on its own terms above and does NOT\n")
  cat("   support the same statement -- do not write \"100-300\". Nothing here\n")
  cat("   licenses a general superiority of the reference distribution.\n")
} else if (r$p < 0.05 && r$mean < 0) {
  cat("-> HARTIGAN is ahead where the paper's claim lives.\n")
  cat("   This belongs in the Results and the Discussion, not only the response.\n")
  cat("   Report it, narrow the contribution to the gate and the two stages, and\n")
  cat("   do NOT switch the reference this close to the deadline: the whole grid\n")
  cat("   was run with the Gaussian one, and a new stage I cannot inherit it.\n")
} else {
  cat("-> Underpowered: not significant, but the point estimate is not small.\n")
  cat("   Report the estimate with its interval and say that nine datasets\n")
  cat("   cannot resolve a difference of this size.\n")
}
cat("\nRegardless of branch, the response and supplement must carry: both budgets,\n",
    "effect sizes with intervals, the per-dataset table, Baron on its own, the\n",
    "Baron-removed sensitivity, and the Array A 1000-feature result.\n", sep = "")

utils::write.csv(m[order(m$n, m$data.name, m$clust),
                   c("data.name","clust","n","ARI_d","ARI_h","dARI")],
                 file.path(out, "hart_cells.csv"), row.names = FALSE)
utils::write.csv(agg, file.path(out, "hart_bydataset.csv"), row.names = FALSE)
cat("\nwritten: ", normalizePath(out), "/hart_cells.csv and hart_bydataset.csv\n", sep = "")
