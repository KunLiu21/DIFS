# difs_rev2_analysis.R: revision-2 analysis source.
# Source: code/as_run/06_revision2_reps_krange_abl2_runtime/difs_rev2_analysis.R
# Usage, scope and limitations: benchmark/REVISION2.md.
difs_args <- function(args, known, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    if (!kv[2] %in% known) stop("unknown argument '", kv[2], "'. Accepts: ", paste(known, collapse = ", "))
    assign(kv[2], kv[3], envir = env)
  }
}
only <- "abl2,reps,krange"; out <- "../results/diag"
difs_args(commandArgs(TRUE), c("only", "out"))
SETS <- trimws(strsplit(only, ",")[[1]])
source("difs_signrank.R")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
NOT_IND <- c("SimKumar4easy", "SimKumar4hard", "Zhengmix4eq", "Zhengmix8eq")

read_set <- function(set) {
  load(sprintf("../scenarios/scenarios_%s.Rdata", set))
  fs <- list.files(file.path("..", "results", set), pattern = "\\.rds$", full.names = TRUE)
  g <- function(x, n, d = NA) { v <- x[[n]]; if (is.null(v) || !length(v)) d else v[1] }
  z <- do.call(rbind, lapply(fs, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(x)) return(NULL)
    data.frame(data.name = g(x, "data.name", ""), arm = g(x, "arm", ""),
               clust = g(x, "clustering.method", ""), n = as.integer(g(x, "n_requested")),
               k_offset = as.integer(g(x, "k_offset")), cseed = as.integer(g(x, "cluster_seed")),
               sub = as.integer(g(x, "subsample_rep", 0L)), seed = as.integer(g(x, "seed")),
               status = g(x, "status", "ok"), ARI = as.numeric(g(x, "ARI")),
               nfeat = as.integer(g(x, "n_features")), stage2_n = as.integer(g(x, "difs_stage2_n")),
               fell_back = isTRUE(g(x, "stage2only_fell_back", FALSE)),
               k_used = as.integer(g(x, "k_used_clustering")), k_found = as.integer(g(x, "n_clusters_found")),
               dropped = as.integer(g(x, "n_cells_dropped")), cells = as.integer(g(x, "n_cells_used")),
               stringsAsFactors = FALSE)
  }))
  ex <- data.frame(data.name = gsub(".*/([^.]*).*", "\\1", scenarios$data.path), arm = scenarios$arm,
                   clust = scenarios$method, n = as.integer(scenarios$n_features),
                   k_offset = as.integer(scenarios$k_offset), cseed = as.integer(scenarios$cluster_seed),
                   sub = as.integer(scenarios$subsample_rep), seed = as.integer(scenarios$seed),
                   stringsAsFactors = FALSE)
  key <- function(d) do.call(paste, c(d[, c("data.name","arm","clust","n","k_offset","cseed","sub","seed")], sep = "|"))
  have <- if (is.null(z)) character() else key(z)
  if (anyDuplicated(have)) stop(set, ": duplicate result files for ", sum(duplicated(have)), " task(s)")
  miss <- ex[!key(ex) %in% have, ]
  complete <- nrow(miss) == 0L
  cat(sprintf("\n%s\n%s: %d of %d tasks present%s\n", strrep("=", 78), set, length(have), nrow(ex),
              if (complete) "" else "  *** INCOMPLETE -- tables below exclude the missing tasks"))
  if (!complete) print(utils::head(miss[, c("data.name","arm","clust","n","k_offset","cseed","sub")], 30), row.names = FALSE)
  if (!is.null(z)) {
    sk <- z[z$status != "ok", ]
    if (nrow(sk)) cat(sprintf("  %d task(s) skipped by design (%s)\n", nrow(sk),
                              paste(unique(sk$status), collapse = ", ")))
    attr_skipped <- nrow(sk)
    bad <- z[z$status == "ok" & (!is.finite(z$ARI) | z$dropped != 0), ]
    if (nrow(bad)) { cat("*** runs with no usable ARI or dropped cells:\n"); print(bad[, 1:10], row.names = FALSE) }
    z <- z[z$status == "ok" & is.finite(z$ARI) & z$dropped == 0, ]
  }
  attr(z, "complete") <- complete; z
}
tag <- function(d) if (isTRUE(attr(d, "complete"))) "" else "  [INCOMPLETE]"
paired <- function(a, b, lab) {
  m <- merge(a, b, by = "data.name", suffixes = c(".a", ".b"))
  if (!nrow(m)) return(invisible())
  r <- difs_signrank(m$ARI.a, m$ARI.b); ci <- difs_ci(m$ARI.a, m$ARI.b)
  cat(sprintf("  %-42s %2d ds  mean %+.4f  [%+.4f, %+.4f]  ahead %2d/%-2d  exact p = %.4f\n",
              lab, nrow(m), mean(m$ARI.a - m$ARI.b), ci[1], ci[2], r$n_pos, r$n, r$p.value))
  invisible(data.frame(comparison = lab, n_ds = nrow(m), mean = mean(m$ARI.a - m$ARI.b),
                       lo = ci[1], hi = ci[2], ahead = r$n_pos, of = r$n, p = r$p.value))
}
per_ds <- function(d) aggregate(ARI ~ data.name, data = d, FUN = mean)

if ("abl2" %in% SETS) {
  z <- read_set("abl2")
  old <- utils::read.csv(file.path(out, "abl_corrected.csv"), stringsAsFactors = FALSE)
  old <- data.frame(data.name = old$data.name, arm = old$fmethod, n = old$n, ARI = old$ARI,
                    stage2_n = old$stage2_n, fell_back = FALSE, stringsAsFactors = FALSE)
  A <- rbind(old, z[, names(old)])
  fb <- z[z$arm == "DIFS_stage2only" & z$fell_back, ]
  if (nrow(fb)) { cat("\nDIFS_stage2only runs where stage II found NO gene (the set is stage I, excluded):\n")
    print(fb[, c("data.name","n","stage2_n")], row.names = FALSE) }
  A <- A[!(A$arm == "DIFS_stage2only" & A$fell_back), ]
  LADDER <- c("GateOnly","DIFS_stage1only","DIFS_fixed","DIFS_fixratio","DIFS",
              "HVG_stage1only","HVG_fixed","HVG_full","DIFS_stage2only")
  res <- list()
  for (N in c(100, 300, 1000)) for (sub in c("13", "9")) {
    s <- A[A$n == N & (sub == "13" | !A$data.name %in% NOT_IND), ]
    cat(sprintf("\n-- abl2, %d features, %s datasets%s --\n", N, sub, tag(z)))
    mm <- tapply(s$ARI, s$arm, mean)[intersect(LADDER, unique(s$arm))]
    cat("  mean ARI: ", paste(sprintf("%s %.3f", names(mm), mm), collapse = " | "), "\n", sep = "")
    arm <- function(a) s[s$arm == a, c("data.name","ARI")]
    C <- list(
      c("DIFS_stage1only","GateOnly",  "dip ranking over the gate"),
      c("HVG_stage1only", "GateOnly",  "HVG ranking over the gate"),
      c("DIFS_stage1only","HVG_stage1only", "stage I: dip minus HVG"),
      c("DIFS_fixed","DIFS_stage1only","+ stage II (no MSE), dip"),
      c("DIFS_fixratio","DIFS_fixed",  "+ MSE ratio search, dip"),
      c("DIFS","DIFS_fixratio",        "+ MSE size sweep, dip"),
      c("HVG_fixed","HVG_stage1only",  "+ stage II (no MSE), HVG"),
      c("HVG_full","HVG_fixed",        "+ MSE (size and ratio), HVG"),
      c("DIFS","HVG_full",             "full pipeline: dip minus HVG"),
      c("DIFS_fixed","HVG_fixed",      "no-MSE pipeline: dip minus HVG"),
      c("DIFS_stage2only","DIFS_stage1only", "stage II only minus stage I only"),
      c("DIFS","DIFS_stage2only",      "combined minus stage II only"))
    for (cc in C) {
      r <- paired(arm(cc[1]), arm(cc[2]), cc[3])
      if (!is.null(r)) res[[length(res) + 1L]] <- cbind(set = "abl2", n = N, datasets = sub, r)
    }
  }
  utils::write.csv(do.call(rbind, res), file.path(out, "rev2_abl2_comparisons.csv"), row.names = FALSE)
}

if ("reps" %in% SETS) {
  z <- read_set("reps")
  g <- utils::read.csv(file.path(out, "grid_corrected.csv"), stringsAsFactors = FALSE)
  g <- g[g$k_policy == "true" & g$n %in% c(100, 300) & g$fmethod %in% c("DIFS","Seurat","FEAST"), ]
  base <- data.frame(data.name = g$data.name, arm = g$fmethod, clust = g$clust, n = g$n,
                     family = "stored", ARI = g$ARI, stringsAsFactors = FALSE)
  z$family <- ifelse(z$sub > 0, "subsample", "seed")
  R <- rbind(base, z[, names(base)])
  cat(sprintf("\n-- reps: spread across replicates (SD of ARI over the runs of one cell)%s --\n", tag(z)))
  sp <- aggregate(ARI ~ data.name + arm + clust + n + family, data = R[R$family != "stored", ],
                  FUN = function(x) c(mean = mean(x), sd = stats::sd(x), n = length(x)))
  sp <- do.call(data.frame, sp)
  names(sp) <- sub("^ARI\\.", "", names(sp))
  print(aggregate(sd ~ arm + clust + family, data = sp, FUN = function(x) round(stats::median(x), 4)), row.names = FALSE)
  utils::write.csv(sp, file.path(out, "rev2_reps_spread.csv"), row.names = FALSE)

  cat("\n-- reps: is the DIFS advantage larger than run-to-run variation? --\n")
  res <- list()
  for (fam in c("seed", "subsample")) for (cl in c("Refined Louvain", "SC3")) for (N in c(100, 300)) {
    s <- R[R$family %in% c(fam) & R$clust == cl & R$n == N, ]
    for (sub in c("13", "9")) {
      s2 <- s[sub == "13" | !s$data.name %in% NOT_IND, ]
      for (cmp in c("Seurat", "FEAST")) {
        r <- paired(per_ds(s2[s2$arm == "DIFS", ]), per_ds(s2[s2$arm == cmp, ]),
                    sprintf("%s %s n=%d %s: DIFS - %s (replicate means)", fam, substr(cl, 1, 7), N, sub, cmp))
        if (!is.null(r)) res[[length(res) + 1L]] <- cbind(family = fam, clust = cl, n = N, datasets = sub, r)
      }
    }
    per <- do.call(rbind, lapply(unique(s$data.name), function(d) {
      a <- s$ARI[s$data.name == d & s$arm == "DIFS"]
      do.call(rbind, lapply(c("Seurat", "FEAST"), function(cmp) {
        b <- s$ARI[s$data.name == d & s$arm == cmp]
        if (length(a) < 3 || length(b) < 3) return(NULL)
        p <- suppressWarnings(stats::wilcox.test(a, b, exact = FALSE)$p.value)
        data.frame(data.name = d, cmp = cmp, diff = mean(a) - mean(b),
                   verdict = if (is.finite(p) && p < 0.05) (if (mean(a) > mean(b)) "DIFS higher" else "DIFS lower") else "not separated")
      }))
    }))
    if (!is.null(per) && nrow(per)) {
      cat(sprintf("  %s %s n=%d, per dataset (two-sample rank test on the 5 replicates, 0.05):\n", fam, cl, N))
      print(table(per$cmp, per$verdict))
    }
  }
  utils::write.csv(do.call(rbind, res), file.path(out, "rev2_reps_comparisons.csv"), row.names = FALSE)
}

if ("krange" %in% SETS) {
  z <- read_set("krange")
  g <- utils::read.csv(file.path(out, "grid_corrected.csv"), stringsAsFactors = FALSE)
  g <- g[g$k_policy == "true" & g$n == 300 & g$fmethod %in% c("DIFS","Seurat","FEAST"), ]
  K <- rbind(data.frame(data.name = g$data.name, arm = g$fmethod, clust = g$clust, off = 0L, ARI = g$ARI),
             data.frame(data.name = z$data.name, arm = z$arm, clust = z$clust, off = z$k_offset, ARI = z$ARI))
  cat(sprintf("\n-- krange: 300 features, k = k_true + offset, every method alike%s --\n", tag(z)))
  res <- list()
  for (cl in c("Refined Louvain", "SC3")) {
    s <- K[K$clust == cl, ]
    mm <- aggregate(ARI ~ arm + off, data = s, FUN = mean)
    w <- stats::reshape(mm, idvar = "arm", timevar = "off", direction = "wide")
    names(w) <- sub("^ARI\\.", "k", names(w)); cat(" ", cl, "mean ARI by offset:\n"); print(w, row.names = FALSE, digits = 3)
    for (o in sort(unique(s$off))) for (cmp in c("Seurat", "FEAST")) {
      r <- paired(s[s$arm == "DIFS" & s$off == o, c("data.name","ARI")],
                  s[s$arm == cmp & s$off == o, c("data.name","ARI")],
                  sprintf("%s offset %+d: DIFS - %s", substr(cl, 1, 7), o, cmp))
      if (!is.null(r)) res[[length(res) + 1L]] <- cbind(clust = cl, offset = o, r)
    }
  }
  cat("  (offsets giving k < 2 are recorded as skipped by the runner and excluded; count above)\n")
  utils::write.csv(do.call(rbind, res), file.path(out, "rev2_krange_comparisons.csv"), row.names = FALSE)
}
cat("\nwritten to", normalizePath(out), ": rev2_*_comparisons.csv, rev2_reps_spread.csv\n")
