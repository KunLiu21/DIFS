# difs_rev2_smoke.R: revision-2 analysis source.
# Source: code/as_run/06_revision2_reps_krange_abl2_runtime/difs_rev2_smoke.R
# Usage, scope and limitations: benchmark/REVISION2.md.
source("Functions_controlled.R")
source("difs_sc3_complete.R")
source("difs_rev2_overlay.R")
ok_all <- TRUE
check <- function(label, cond, detail = "") {
  cat(sprintf("  [%s] %s  %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label, detail))
  if (!isTRUE(cond)) ok_all <<- FALSE
}
koh <- list.files("../source", pattern = "^Koh\\.rds$", full.names = TRUE)
if (length(koh) != 1L) stop("../source/Koh.rds not found")

row <- function(arm, fm, n, method = "Refined Louvain", n_policy = "original",
                stage1 = "dip", final = "mixed", cseed = NA_integer_, sub = 0L,
                frac = 1, koff = NA_integer_, seed = 1)
  data.frame(arm = arm, feature.selection.method = fm, n_policy = n_policy,
             stage1_source = stage1, final_set = final, n_features = n, data.path = koh,
             method = method, k_policy = "true", seed = seed, min.expression = log(5),
             gate_rule = "submitted", cluster_seed = cseed, subsample_rep = sub,
             subsample_frac = frac, k_offset = koff, stringsAsFactors = FALSE)
scenarios <- rbind(
  row("DIFS", "DIFS", 100), row("Seurat", "Seurat", 100),                          # 1, 2
  row("HVG_stage1only", "DIFS_stage1only", 100, stage1 = "hvg"),                   # 3
  row("DIFS_stage2only", "DIFS", 100, final = "stage2only"),                       # 4
  row("DIFS_fixed", "DIFS", 100, n_policy = "fixed_all"),                          # 5
  row("DIFS_fixratio", "DIFS", 100, n_policy = "fixed_stage1"),                    # 6
  row("Seurat", "Seurat", 300, koff = 1L),                                         # 7
  row("Seurat", "Seurat", 100, sub = 1L, frac = 0.9),                              # 8
  row("Seurat", "Seurat", 100, method = "SC3", cseed = 101L, seed = 101),          # 9
  row("Seurat", "Seurat", 100, cseed = 101L, seed = 101),                          # 10
  row("Seurat", "Seurat", 100, cseed = 102L, seed = 102),                          # 11
  row("HVG_full", "DIFS", 100, stage1 = "hvg"))                                    # 12
prepDir("../scenarios"); save(scenarios, file = "../scenarios/scenarios_rev2smoke.Rdata")
unlink("../results/rev2smoke", recursive = TRUE)
for (i in seq_len(nrow(scenarios))) {
  st <- system2("Rscript", c("Method_comparison_rev2.R", "set=rev2smoke", paste0("k=", i)))
  if (st != 0) { cat("  [FAIL] task", i, "exited with status", st, "\n"); ok_all <- FALSE }
}
R <- lapply(list.files("../results/rev2smoke", full.names = TRUE), readRDS)
get_run <- function(arm, n, method = "Refined Louvain", cs = NA, sub = 0L, koff = NA) {
  hit <- Filter(function(r) identical(r$arm, arm) && r$n_requested == n &&
                  identical(r$clustering.method, method) &&
                  identical(is.na(r$cluster_seed), is.na(cs)) &&
                  (is.na(cs) || r$cluster_seed == cs) && r$subsample_rep == sub &&
                  identical(is.na(r$k_offset), is.na(koff)) && (is.na(koff) || r$k_offset == koff), R)
  if (length(hit) != 1L) stop("expected one run for ", arm, " n=", n, ", found ", length(hit))
  hit[[1]]
}
stored <- function(fm, n) {
  f <- Sys.glob(sprintf("../results/ncurve/Koh_%s_Refined Louvain_n%04d_k-true*seed-1.rds", fm, n))
  if (length(f) != 1L) stop("stored Koh ", fm, " n=", n, ": found ", length(f), " files")
  readRDS(f)
}

cat("\n1  defaults reproduce the stored runs\n")
for (fm in c("DIFS", "Seurat")) {
  a <- get_run(fm, 100)$ARI; b <- stored(fm, 100)$ARI
  check(sprintf("%s ARI identical to results/ncurve", fm), isTRUE(all.equal(a, b, tolerance = 1e-12)),
        sprintf("new %.10f  stored %.10f", a, b))
}

cat("\n2  the clustering seed arrives\n")
seu <- readRDS(koh)
top <- VariableFeatures(FindVariableFeatures(seu, nfeatures = 100, verbose = FALSE))
options(difs.cluster_seed = 777L)
s1 <- seudo_clustering(seu, top, cluster_count = 9, clustering.method = "Refined Louvain")
cmd <- s1@commands
pca <- cmd[[grep("^RunPCA", names(cmd))[1]]]; fc <- cmd[[grep("^FindClusters", names(cmd))[1]]]
check("RunPCA seed.use = option", identical(as.integer(pca$seed.use), 777L), paste("got", pca$seed.use))
check("FindClusters random.seed = option", identical(as.integer(fc$random.seed), 777L), paste("got", fc$random.seed))
difs_rev2_options_reset()
s0 <- seudo_clustering(seu, top, cluster_count = 9, clustering.method = "Refined Louvain")
cmd0 <- s0@commands
check("defaults unchanged (42 / 0)",
      identical(as.integer(cmd0[[grep("^RunPCA", names(cmd0))[1]]]$seed.use), 42L) &&
      identical(as.integer(cmd0[[grep("^FindClusters", names(cmd0))[1]]]$random.seed), 0L))
l101 <- get_run("Seurat", 100, cs = 101); l102 <- get_run("Seurat", 100, cs = 102)
cat(sprintf("     (information) Seurat/Louvain ARI at seed 101: %.4f   102: %.4f   stored: %.4f\n",
            l101$ARI, l102$ARI, stored("Seurat", 100)$ARI))

cat("\n3  HVG stage I\n")
h <- get_run("HVG_stage1only", 100)
gate <- difs_gate_genes(as.matrix(GetAssayData(seu, slot = "data")), log(5))
hv <- HVFInfo(FindVariableFeatures(seu, selection.method = "vst", nfeatures = nrow(seu), verbose = FALSE))
v <- setNames(hv[[grep("variance\\.standardized$", colnames(hv), value = TRUE)[1]]], rownames(hv))
check("100 features returned", h$n_features == 100, paste("got", h$n_features))
check("stage II off", identical(h$difs_use_stage2, FALSE))
cat(sprintf("     gate %d genes; ARI %.4f\n", length(gate), h$ARI))

cat("\n4  stage II only\n")
s2 <- get_run("DIFS_stage2only", 100)
check("stage II found genes", isTRUE(s2$difs_stage2_n > 0), paste("stage2_n", s2$difs_stage2_n))
check("did not fall back to stage I", identical(s2$stage2only_fell_back, FALSE))
check("at most 100 features", s2$n_features <= 100, paste("got", s2$n_features))

cat("\n5  MSE off / ratio only\n")
f1 <- get_run("DIFS_fixed", 100); f2 <- get_run("DIFS_fixratio", 100); full <- get_run("DIFS", 100)
check("DIFS_fixed: stage I seed set = n (no size sweep)", f1$difs_stage1_n == 100, paste("stage1_n", f1$difs_stage1_n))
check("DIFS_fixed: 1:1 mix", identical(f1$difs_ratio, "1:1"), paste("ratio", f1$difs_ratio))
check("DIFS_fixratio: stage I seed set = n, ratio searched", f2$difs_stage1_n == 100,
      paste("stage1_n", f2$difs_stage1_n, "ratio", f2$difs_ratio))
cat(sprintf("     full DIFS: stage1_n %d, ratio %s\n", full$difs_stage1_n, full$difs_ratio))

cat("\n6  k override\n")
ko <- get_run("Seurat", 300, koff = 1)
check("k used = k_true + 1", ko$k_used_clustering == ko$k_true_labels + 1,
      sprintf("k_true %d, used %d", ko$k_true_labels, ko$k_used_clustering))

cat("\n7  stratified 90% subsample\n")
sb <- get_run("Seurat", 100, sub = 1L)
tab <- table(addNA(factor(seu$trueclass))); tab <- tab[tab > 0]
check("cell count = sum of per-class rounded 90%", sb$n_cells_used == sum(pmax(1, round(0.9 * tab))),
      sprintf("used %d of %d", sb$n_cells_used, ncol(seu)))

cat("\n8  SC3 with a seed\n")
sc <- get_run("Seurat", 100, method = "SC3", cs = 101)
check("SC3 labels every cell", sc$n_cells_dropped == 0, paste("dropped", sc$n_cells_dropped))

cat("\n12 HVG full pipeline\n")
hf <- get_run("HVG_full", 100)
check("runs and returns 100", hf$n_features == 100, sprintf("ARI %.4f, stage2_n %s", hf$ARI, hf$difs_stage2_n))

cat("\n", if (ok_all) "ALL CHECKS PASSED" else "*** SOME CHECKS FAILED -- do not submit the arrays", "\n", sep = "")
quit(save = "no", status = if (ok_all) 0 else 1)
