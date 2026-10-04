# difs_key_check.R: revision-2 analysis source.
# Source: code/as_run/07_feature_sets_gate_ranking/difs_key_check.R
# Usage, scope and limitations: benchmark/REVISION2.md.
source("Functions_controlled.R")

args <- commandArgs(TRUE)
for (a in args) {
  kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
  if (length(kv) == 3L) assign(kv[2],
    if (grepl("^-?[0-9.]+$", kv[3])) as.numeric(kv[3]) else kv[3], envir = globalenv())
}
if (!exists("dataset", inherits = FALSE)) dataset <- "Koh"
if (!exists("out_dir", inherits = FALSE)) out_dir <- "../results/key_check"
if (!exists("min.expression", inherits = FALSE)) min.expression <- log(5)
if (!exists("topN",    inherits = FALSE)) topN <- 1000

KEYS <- c("normal_lookup", "tnormal_lookup", "uniform_lookup",
          "hartigan", "mc_normal", "D_raw")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))

for (dn in trimws(strsplit(dataset, ",")[[1]])) {
  pth <- file.path("../source", paste0(dn, ".rds"))
  if (!file.exists(pth)) { message("missing: ", pth); next }
  message("\n########## ", dn, " ##########")
  seu <- readRDS(pth)
  m   <- as.matrix(GetAssayData(seu, slot = "data"))
  g   <- difs_gate_genes(m, min.expression)
  P   <- length(g)
  cat(sprintf("%d cells, %d genes, %d classes; %d genes pass the gate (%.1f%%)\n",
              ncol(m), nrow(m), length(unique(seu$trueclass)), P, 100 * P / nrow(m)))

  Dn <- vapply(g, function(gn) { x <- m[gn, ]; x <- x[x > min.expression]
                                 c(diptest::dip(sort(x)), length(x)) }, numeric(2))
  cat(sprintf("expressing cells among tested genes: median %.0f, range %d-%d\n",
              stats::median(Dn[2, ]), min(Dn[2, ]), max(Dn[2, ])))

  ord <- list(); val <- list()
  for (k in KEYS) {
    t0 <- proc.time()[["elapsed"]]
    r  <- difs_stage1_ranking(m, min.expression, key = k, seed = 1)
    el <- proc.time()[["elapsed"]] - t0
    ord[[k]] <- r
    v <- switch(k,
      mc_normal = NA, hartigan = NA,
      D_raw = -Dn[1, ],
      dip_pvalue_lookup_vec(Dn[1, ], Dn[2, ], difs_null_tables(sub("_lookup$", "", k))))
    val[[k]] <- v
    cat(sprintf("  %-15s %6.1f s\n", k, el))
  }
  rep2 <- list(mc_normal     = difs_stage1_ranking(m, min.expression, key = "mc_normal", seed = 2),
               normal_lookup = difs_stage1_ranking(m, min.expression, key = "normal_lookup", seed = 2))

  cat("\n---------- ties ----------\n")
  tie <- do.call(rbind, lapply(KEYS, function(k) {
    v <- val[[k]]
    if (all(is.na(v))) {
      v <- if (k == "hartigan")
             vapply(g, function(gn){x <- m[gn,]; suppressWarnings(diptest::dip.test(x[x>min.expression])$p.value)}, 0)
           else NA
    }
    if (all(is.na(v))) return(data.frame(key = k, distinct = NA, at_min = NA, at_max = NA))
    data.frame(key = k, distinct = length(unique(v)),
               at_min = sum(v == min(v)), at_max = sum(v == max(v)))
  }))
  tie$n_genes <- P; tie$distinct_frac <- round(tie$distinct / P, 4)
  print(tie, row.names = FALSE)
  cat("  at_min is the number of genes tied at the END THAT DECIDES SELECTION;\n")
  cat("  sort() breaks those by matrix row order, not by evidence.\n")

  cat("\n---------- reproducibility (same key, different seed) ----------\n")
  for (k in names(rep2))
    cat(sprintf("  %-15s top-%d Jaccard vs rerun: %.4f\n", k, topN,
                jacc(head(ord[[k]], topN), head(rep2[[k]], topN))))

  cat("\n---------- chance-corrected top-", topN, " agreement ----------\n", sep = "")
  ch <- topN / (2 * P - topN)
  cat(sprintf("  (chance Jaccard for N=%d from a pool of %d is %.3f)\n", topN, P, ch))
  ks <- KEYS
  A <- matrix(NA_real_, length(ks), length(ks), dimnames = list(ks, ks))
  for (i in seq_along(ks)) for (j in seq_along(ks)) {
    jj <- jacc(head(ord[[ks[i]]], topN), head(ord[[ks[j]]], topN))
    A[i, j] <- round((jj - ch) / (1 - ch), 3)
  }
  print(A)
  saveRDS(list(order = ord, values = val, tie = tie, agreement = A, gate_n = P),
          file.path(out_dir, paste0("key_check_", dn, ".rds")))
}
cat("\nwritten to: ", normalizePath(out_dir), "\n", sep = "")
