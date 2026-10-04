# Regenerate corrected summaries from public per-run/per-dataset CSVs; no clustering rerun.
script <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))
here <- dirname(dirname(script))
args <- commandArgs(TRUE)
if (length(args)>1) stop("usage: Rscript scripts/regen_tables_20261004.R [new-output-directory]")
out <- if(length(args)) args[1] else file.path(here,"work/recomputed-tables")
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop("Output directory must be empty")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(here, "benchmark/difs_signrank.R"))
NOT_IND <- c("SimKumar4easy", "SimKumar4hard", "Zhengmix4eq", "Zhengmix8eq")

old <- read.csv(file.path(here, "results/revision2/02_ablation/abl_corrected.csv"))
old <- data.frame(data.name = old$data.name, arm = old$fmethod, n = old$n, ARI = old$ARI)
s0b <- read.csv(file.path(here,"results/revision2/supplementary_tables_S0-S5/S0b_rev2_per_run.csv"))
s0b <- s0b[s0b$result_set=="abl2" & s0b$status=="ok", ]
if(any(s0b$stage2only_fell_back %in% TRUE)) stop("Stage-II fallback rows require explicit exclusion")
new <- data.frame(data.name = s0b$dataset, arm = s0b$arm, n = s0b$n_requested, ARI = s0b$ARI)
A <- rbind(old, new)
paired <- function(a, b, lab) {
  m <- merge(a, b, by = "data.name", suffixes = c(".a", ".b"))
  r <- difs_signrank(m$ARI.a, m$ARI.b); ci <- difs_ci(m$ARI.a, m$ARI.b)
  data.frame(comparison = lab, n_ds = nrow(m), mean = mean(m$ARI.a - m$ARI.b),
             lo = ci[1], hi = ci[2], ahead = r$n_pos, of = r$n, p = r$p.value)
}
CMP <- list(
  c("DIFS_stage1only","GateOnly","dip ranking over the gate"), c("HVG_stage1only","GateOnly","HVG ranking over the gate"),
  c("DIFS_stage1only","HVG_stage1only","stage I: dip minus HVG"), c("DIFS_fixed","DIFS_stage1only","+ stage II (no MSE), dip"),
  c("DIFS_fixratio","DIFS_fixed","+ MSE ratio search, dip"), c("DIFS","DIFS_fixratio","+ MSE size sweep, dip"),
  c("HVG_fixed","HVG_stage1only","+ stage II (no MSE), HVG"), c("HVG_full","HVG_fixed","+ MSE (size and ratio), HVG"),
  c("DIFS","HVG_full","full pipeline: dip minus HVG"), c("DIFS_fixed","HVG_fixed","no-MSE pipeline: dip minus HVG"),
  c("DIFS_stage2only","DIFS_stage1only","stage II only minus stage I only"), c("DIFS","DIFS_stage2only","combined minus stage II only"))
res <- list()
for (N in c(100, 300, 1000)) for (sub in c("13", "9")) {
  s <- A[A$n == N & (sub == "13" | !A$data.name %in% NOT_IND), ]
  arm <- function(a) s[s$arm == a, c("data.name", "ARI")]
  for (cc in CMP) res[[length(res) + 1]] <- cbind(set = "abl2", n = N, datasets = sub, paired(arm(cc[1]), arm(cc[2]), cc[3]))
}
ab <- do.call(rbind, res)
write.csv(ab, file.path(out, "rev2_abl2_comparisons.csv"), row.names = FALSE)
cat("abl2 comparisons:", nrow(ab), "rows written\n")

signrank_exact <- function(d) {
  d <- round(d[is.finite(d)], 12); d <- d[d != 0]; m <- length(d)
  if (m < 2) return(NA_real_)
  r <- rank(abs(d)); obs <- sum(r[d > 0])
  S <- as.matrix(expand.grid(rep(list(c(0, 1)), m))); null <- S %*% r
  mean(abs(null - sum(r) / 2) >= abs(obs - sum(r) / 2) - 1e-9)
}
C <- read.csv(file.path(here, "results/revision2/08_markers/marker_coverage_per_dataset.csv"))
INDEPENDENT <- c("Kumar","Trapnell","Romanov","Lawlor","Koh","Darmanis","Muraro","Fletcher","Baron")
BUDGETS <- sort(unique(C$budget)); ind <- C[C$dataset %in% INDEPENDENT, ]
tests <- list()
for (L in unique(C$list)) for (lab in c("independent9", "all"))
  for (cmp in c("Seurat", "FEAST", "Seurat_gated", "DIFS_stage1only"))
  for (b in BUDGETS) for (v in c("cov10", "types1", "precision")) {
    sub <- if (lab == "all") C else ind; sub <- sub[sub$list == L, ]
    a <- sub[sub$method == "DIFS" & sub$budget == b, c("dataset", v)]
    z <- sub[sub$method == cmp & sub$budget == b, c("dataset", v)]
    m <- merge(a, z, by = "dataset"); if (nrow(m) < 2) next
    dd <- round(m[[2]] - m[[3]], 12)
    tests[[length(tests) + 1]] <- data.frame(list = L, datasets = lab, vs = cmp, budget = b, metric = v,
      n = nrow(m), mean_DIFS = mean(m[[2]]), mean_other = mean(m[[3]]), mean_diff = mean(dd),
      wins = sum(dd > 0), losses = sum(dd < 0), p_exact = signrank_exact(dd),
      primary = L == "hurdle" && lab == "independent9" && cmp %in% c("Seurat", "FEAST") && b <= 300)
  }
Tt <- do.call(rbind, tests)
write.csv(Tt, file.path(out, "marker_coverage_tests.csv"), row.names = FALSE)
cat("marker tests:", nrow(Tt), "rows written\n")
