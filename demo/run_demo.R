################################################################################
## run_demo.R -- run DIFS end to end on a tiny dataset, in about ten seconds
##
##   Rscript demo/run_demo.R                    (from the repository root)
##
## TIER 1 always runs.  It needs R >= 4.1 and the 'diptest' package and nothing
## else -- no Bioconductor, no network.  It executes the DIFS stage I criterion
## (expression gate + dip ranking against a precomputed normal null) and checks
## that the genes it puts first are the ones that are informative by
## construction.
##
## TIER 2 runs only if Seurat is installed.  It clusters on the selected
## features and reports the adjusted Rand index against the known labels, which
## is the quantity the manuscript reports.  If Seurat is absent the script says
## so and exits 0 -- a missing optional dependency is not a failure.
##
## What this demo is NOT: it is not one of the paper's datasets and it is not
## evidence for any claim in the paper.  It exists so that a reader can confirm
## the code runs before committing to the full benchmark, which needs a cluster.
################################################################################

ok <- function(msg) cat("  [ok]   ", msg, "\n", sep = "")
bad <- function(msg) { cat("  [FAIL] ", msg, "\n", sep = ""); quit(status = 1) }
hr <- function(s) cat("\n", s, "\n", strrep("-", nchar(s)), "\n", sep = "")

if (!file.exists("R/difs_core.R"))
  stop("run this from the repository root: Rscript demo/run_demo.R")
source("R/difs_core.R")
d <- readRDS("demo/demo_data.rds")

cat("DIFS demo\n")
cat(sprintf("R %s | diptest %s\n", getRversion(),
            as.character(utils::packageVersion("diptest"))))
cat(sprintf("data: %d genes x %d cells, %d cell types (simulated, see demo/make_demo_data.R)\n",
            nrow(d$counts), ncol(d$counts), d$k))

## ------------------------------------------------------------------ TIER 1 --
hr("Tier 1  stage I feature selection  (diptest only)")

t0 <- proc.time()[["elapsed"]]
gate <- difs_gate_genes(d$logcounts, min.expression = log(5),
                        min.expression.proportion = 0.05,
                        gate_rule = "combined", gate_k = d$k)
thr <- difs_gate_threshold(ncol(d$logcounts), k = d$k, rule = "combined")
cat(sprintf("expression gate: a gene must exceed log(5) in >= %.0f of %d cells\n",
            thr, ncol(d$logcounts)))
cat(sprintf("                 %d of %d genes pass\n", length(gate), nrow(d$logcounts)))
if (!length(gate)) bad("the gate passed no gene")

ranked <- difs_stage1_ranking(d$logcounts, min.expression = log(5),
                              gate_rule = "combined", gate_k = d$k,
                              key = "normal_lookup")
el <- proc.time()[["elapsed"]] - t0
cat(sprintf("stage I ranking: %d genes ordered in %.2f s\n", length(ranked), el))

## The demo dataset contains three kinds of gene, and what stage I does with
## each is the point of the demonstration.
rank_of <- function(genes) {
  inr <- intersect(genes, ranked)
  list(n_gated = length(inr), in_top30 = sum(head(ranked, 30) %in% inr),
       median_rank = if (length(inr)) stats::median(match(inr, ranked)) else NA_real_)
}
b <- rank_of(d$bimodal); o <- rank_of(d$onoff); g <- rank_of(d$background)

cat("\n")
cat(sprintf("%-28s %8s %10s %12s\n", "gene class", "gated", "in top 30", "median rank"))
cat(sprintf("%-28s %4d/%-3d %10d %12s\n", "BIMODAL (stage I target)",
            b$n_gated, length(d$bimodal), b$in_top30, format(b$median_rank)))
cat(sprintf("%-28s %4d/%-3d %10d %12s\n", "ONOFF (stage II target)",
            o$n_gated, length(d$onoff), o$in_top30, format(o$median_rank)))
cat(sprintf("%-28s %4d/%-3d %10d %12s\n", "background",
            g$n_gated, length(d$background), g$in_top30, format(g$median_rank)))
cat("\nfirst 10 by stage I: ", paste(head(ranked, 10), collapse = " "), "\n", sep = "")

## Checks.  These are properties of the method, not tuned thresholds.
if (length(ranked) != length(gate)) bad("ranking and gate disagree on gene count")
if (anyDuplicated(ranked))          bad("the ranking contains duplicates")

if (b$in_top30 < 24)
  bad(sprintf("stage I put only %d of the 30 bimodal genes in the top 30", b$in_top30))
ok(sprintf("stage I recovers %d/30 bimodal genes in the top 30 (median rank %s)",
           b$in_top30, format(b$median_rank)))

## This is the interesting one.  An on/off marker is expressed in one cell type
## and absent elsewhere, so the sample stage I actually tests -- the cells above
## the expression gate -- is UNIMODAL.  Stage I therefore ranks these genes no
## better than background, and that is not a defect: it is the reason DIFS has a
## second stage, which finds cell-type-specific genes by a Fisher test against a
## preliminary clustering rather than by multimodality.  If this check ever
## fails, the demo data no longer separates the two regimes.
if (o$in_top30 > 3)
  bad(sprintf("%d on/off genes reached the top 30; the two regimes are not separated",
              o$in_top30))
ok(sprintf("stage I does NOT rank on/off markers (median rank %s vs %s for background)",
           format(o$median_rank), format(g$median_rank)))
cat("        ^ this is by design and is why DIFS has a second stage;\n")
cat("          see the manuscript, Results, 'Overview'.\n")

r2 <- difs_stage1_ranking(d$logcounts, min.expression = log(5),
                          gate_rule = "combined", gate_k = d$k,
                          key = "normal_lookup")
if (!identical(ranked, r2)) bad("the ranking is not deterministic")
ok("the ranking is deterministic across runs")

## The lookup exists to remove the ties the per-gene Monte Carlo produced.
vals <- vapply(gate, function(g) {
  x <- d$logcounts[g, ]; x <- x[x > log(5)]
  diptest::dip(sort(x))
}, numeric(1))
cat(sprintf("\ndip statistic: %d distinct values among %d gated genes\n",
            length(unique(vals)), length(vals)))
ok("tier 1 complete")

## ------------------------------------------------------------------ TIER 2 --
hr("Tier 2  clustering and ARI  (needs Seurat)")

if (!requireNamespace("Seurat", quietly = TRUE)) {
  cat("Seurat is not installed, so the clustering step is skipped.\n")
  cat("Tier 1 above is the method itself and it ran.  To run tier 2:\n")
  cat("  install.packages(\"Seurat\")\n")
  cat("\nDEMO PASSED (tier 1 only)\n")
  quit(status = 0)
}

suppressPackageStartupMessages(library(Seurat))
ari <- function(a, b) {                      # adjusted Rand index, no mclust
  t <- table(a, b); n <- sum(t)
  s <- sum(choose(t, 2)); a_ <- sum(choose(rowSums(t), 2))
  b_ <- sum(choose(colSums(t), 2)); e <- a_ * b_ / choose(n, 2)
  (s - e) / ((a_ + b_) / 2 - e)
}
cluster_on <- function(features) {
  so <- CreateSeuratObject(counts = d$counts)
  so <- NormalizeData(so, verbose = FALSE)
  VariableFeatures(so) <- features
  so <- ScaleData(so, features = features, verbose = FALSE)
  so <- RunPCA(so, features = features, npcs = min(20, length(features) - 1),
               verbose = FALSE)
  so <- FindNeighbors(so, dims = seq_len(min(10, ncol(so@reductions$pca))),
                      verbose = FALSE)
  so <- FindClusters(so, resolution = 0.5, verbose = FALSE)
  as.character(so$seurat_clusters)
}
for (n_feat in c(30, 100)) {
  cl <- cluster_on(head(ranked, n_feat))
  cat(sprintf("DIFS top %-3d features -> %d clusters, ARI %.3f\n",
              n_feat, length(unique(cl)), ari(cl, d$truth[names(d$truth)])))
}
cl_all <- cluster_on(rownames(d$logcounts))
cat(sprintf("all %d genes          -> %d clusters, ARI %.3f\n",
            nrow(d$logcounts), length(unique(cl_all)),
            ari(cl_all, d$truth[names(d$truth)])))
ok("tier 2 complete")
cat("\nDEMO PASSED\n")
