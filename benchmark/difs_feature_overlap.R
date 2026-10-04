# difs_feature_overlap.R: revision-2 analysis source.
# Source: code/as_run/07_feature_sets_gate_ranking/difs_feature_overlap.R
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
base_dir <- "../results/feature_sets"; v2_dir <- "../results/feature_sets_v2"
n_fixed <- 1000
difs_args(commandArgs(TRUE), c("base_dir", "v2_dir", "n_fixed"))

overlap_stats <- function(a, b, U) {
  a2 <- intersect(a, U); b2 <- intersect(b, U)
  N <- length(U); n1 <- length(a2); n2 <- length(b2)
  x <- length(intersect(a2, b2))
  if (n1 == 0L || n2 == 0L) stop("an empty selection after restriction to the common universe")
  k  <- seq.int(max(0L, n1 + n2 - N), min(n1, n2))
  pk <- stats::dhyper(k, n1, N - n1, n2)
  EJ <- sum(pk * k / (n1 + n2 - k))
  J  <- x / (n1 + n2 - x)
  data.frame(universe = N, n1_raw = length(a), n2_raw = length(b), n1 = n1, n2 = n2,
             shared = x, expected = n1 * n2 / N, fold = x / (n1 * n2 / N),
             p_enrich = stats::phyper(x - 1L, n1, N - n1, n2, lower.tail = FALSE),
             jaccard = J, jaccard_chance = EJ, jaccard_adj = (J - EJ) / (1 - EJ))
}

local({
  set.seed(20260924)
  U <- sprintf("g%05d", 1:15000)
  sims <- replicate(400, {
    a <- sample(U, 1000); b <- sample(U, 700)
    overlap_stats(a, b, U)$jaccard_adj
  })
  st <- overlap_stats(sample(U, 1000), sample(U, 700), U)
  mc <- mean(replicate(4000, { a <- sample(U, 1000); b <- sample(U, 700)
                               x <- length(intersect(a, b)); x / (1700 - x) }))
  cat(sprintf("self-test: mean J_adj over random pairs = %+.5f (must be ~0)\n", mean(sims)))
  cat(sprintf("           exact E[J] = %.6f   simulated = %.6f\n", st$jaccard_chance, mc))
  if (abs(mean(sims)) > 0.003 || abs(st$jaccard_chance - mc) > 0.001)
    stop("self-test failed -- the chance correction is wrong; nothing written")
})

base_sets <- readRDS(file.path(base_dir, "feature_sets.rds"))
base_info <- utils::read.csv(file.path(base_dir, "feature_sets_info.csv"), stringsAsFactors = FALSE)
U_path <- file.path(v2_dir, "universes.rds")
if (!file.exists(U_path)) stop("gene universes not found: ", U_path,
                               " -- run difs_export_features_v2.R task=universe")
UNI <- readRDS(U_path)

REPAIRED <- c("Baron|DIFS", "Baron|DIFS_hartigan")
sets <- base_sets; prov <- list()
for (key in REPAIRED) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  f <- file.path(v2_dir, sprintf("select__%s__%s.rds", parts[1], parts[2]))
  if (!file.exists(f)) stop("repaired selection not found: ", f,
                            ". The v1 ", key, " set was chosen on 5000 of 8569 cells ",
                            "and may not be used.")
  r <- readRDS(f)
  if (!isTRUE(r$sc3_repair)) stop(f, " was not produced with the SC3 repair")
  if (length(r$genes) != n_fixed) stop(f, " holds ", length(r$genes), " genes, not ", n_fixed)
  old <- base_sets[[key]]
  prov[[key]] <- data.frame(entry = key, source = "v2_repaired",
                            n_old = length(old), n_new = length(r$genes),
                            shared_old_new = length(intersect(old, r$genes)),
                            stage1_n = r$stage1_n, stage2_n = r$stage2_n, ratio = r$ratio,
                            elapsed_sec = round(r$elapsed_sec), stringsAsFactors = FALSE)
  sets[[key]] <- r$genes
}
for (key in setdiff(names(sets), REPAIRED))
  prov[[key]] <- data.frame(entry = key, source = "v1", n_old = length(sets[[key]]),
                            n_new = length(sets[[key]]), shared_old_new = length(sets[[key]]),
                            stage1_n = NA, stage2_n = NA, ratio = NA, elapsed_sec = NA,
                            stringsAsFactors = FALSE)
prov <- do.call(rbind, prov)

bad <- names(sets)[vapply(sets, length, 1L) != n_fixed]
if (length(bad)) {
  cat("\nentries EXCLUDED because they do not hold ", n_fixed, " genes:\n", sep = "")
  for (b in bad) cat(sprintf("  %-26s %d genes\n", b, length(sets[[b]])))
}

METHODS <- c("DIFS", "DIFS_hartigan", "Seurat", "FEAST", "GateOnly", "Variance", "Random")
dsets <- names(UNI)

rows <- list()
for (m in METHODS) for (i in seq_along(dsets)) for (j in seq_along(dsets)) {
  if (j <= i) next
  d1 <- dsets[i]; d2 <- dsets[j]
  k1 <- paste(d1, m, sep = "|"); k2 <- paste(d2, m, sep = "|")
  if (!k1 %in% names(sets) || !k2 %in% names(sets)) next
  base <- data.frame(method = m, d1 = d1, d2 = d2,
                     id1 = UNI[[d1]]$idspace, id2 = UNI[[d2]]$idspace,
                     stringsAsFactors = FALSE)
  if (!identical(UNI[[d1]]$idspace, UNI[[d2]]$idspace) || k1 %in% bad || k2 %in% bad) {
    rows[[length(rows) + 1L]] <- cbind(base, comparable = FALSE,
      reason = if (k1 %in% bad || k2 %in% bad) "incomplete selection" else "identifier spaces differ")
    next
  }
  U <- intersect(UNI[[d1]]$ids, UNI[[d2]]$ids)
  if (length(U) < 1000) stop(d1, " and ", d2, " are both '", UNI[[d1]]$idspace,
                             "' but share only ", length(U), " genes -- identifier detection is wrong")
  rows[[length(rows) + 1L]] <- cbind(base, comparable = TRUE, reason = "",
                                     overlap_stats(sets[[k1]], sets[[k2]], U))
}
ov <- do.call(rbind, lapply(rows, function(r) {
  full <- c("method","d1","d2","id1","id2","comparable","reason","universe","n1_raw","n2_raw",
            "n1","n2","shared","expected","fold","p_enrich","jaccard","jaccard_chance","jaccard_adj")
  for (c in setdiff(full, names(r))) r[[c]] <- NA
  r[, full]
}))

rnd <- ov[ov$method == "Random" & ov$comparable, "jaccard_adj"]
cat(sprintf("\nRandom arm, chance-corrected Jaccard over %d comparable pair(s): %s\n",
            length(rnd), paste(sprintf("%+.4f", rnd), collapse = "  ")))
if (length(rnd) && any(abs(rnd) > 0.05))
  cat("  *** a Random pair is far from zero. The correction or the universe is wrong.\n")

dir.create(v2_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(sets, file.path(v2_dir, "feature_sets_merged.rds"))
utils::write.csv(prov, file.path(v2_dir, "feature_sets_provenance.csv"), row.names = FALSE)
utils::write.csv(ov,   file.path(v2_dir, "feature_set_overlap_v2.csv"), row.names = FALSE)

cat("\n========== repaired Baron entries: how much the repair changed them ==========\n")
print(prov[prov$source == "v2_repaired", ], row.names = FALSE)

cat("\n========== comparable pairs, chance-corrected on the common universe ==========\n")
show <- ov[ov$comparable, c("method","d1","d2","universe","n1","n2","shared","expected",
                            "fold","jaccard","jaccard_adj","p_enrich")]
show$expected <- round(show$expected, 1); show$fold <- round(show$fold, 2)
show$jaccard <- round(show$jaccard, 4); show$jaccard_adj <- round(show$jaccard_adj, 4)
show$p_enrich <- signif(show$p_enrich, 3)
print(show[order(show$d1, show$d2, -show$jaccard_adj), ], row.names = FALSE)

cat("\nNot comparable (reported, not scored as zero):\n")
nc <- unique(ov[!ov$comparable, c("d1","d2","id1","id2","reason")])
print(nc, row.names = FALSE)

cat("\nwritten to ", normalizePath(v2_dir), ":\n",
    "  feature_sets_merged.rds       v1 sets with the two Baron entries replaced\n",
    "  feature_sets_provenance.csv   which entry came from where\n",
    "  feature_set_overlap_v2.csv    every pair, every method\n", sep = "")
