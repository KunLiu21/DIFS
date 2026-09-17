################################################################################
## make_demo_data.R -- build the tiny dataset that demo/run_demo.R uses
##
##   Rscript demo/make_demo_data.R      (from the repository root)
##
## The demo must run on a laptop with no network access and no Bioconductor, so
## it cannot download Kumar or Baron.  This script simulates a small panel with
## the structure DIFS is designed for and commits the result, so the demo is
## reproducible and offline.  It is NOT one of the paper's datasets and no
## result in the paper comes from it; it exists so that a reader can execute the
## method and see it work before deciding whether to run the real benchmark.
##
## Structure: 3 cell types x 60 cells, 600 genes, of which
##   * 30 BIMODAL genes -- expressed in every type but at two clearly separated
##     levels.  Among the cells that express them the distribution has two
##     modes, so these are what STAGE I is built to find.
##   * 30 ONOFF genes -- expressed in one type and essentially absent elsewhere.
##     Among the cells that express them the distribution is UNIMODAL, so
##     stage I cannot rank them highly and is not supposed to.  These are the
##     genes STAGE II exists to recover; the manuscript makes exactly this
##     distinction ("cell type-specific genes also have a unimodal
##     distribution"), and demo/run_demo.R checks that the split comes out this
##     way, which is the clearest possible demonstration of why DIFS has two
##     stages rather than one.
##   * 540 background genes with no type structure.
## Library sizes vary 3-fold so that normalisation is not a no-op.
################################################################################
set.seed(20260917)

K <- 3; per_type <- 60; n_cells <- K * per_type
n_bg <- 540; n_onoff <- 30; n_bimodal <- 30
type <- factor(rep(paste0("type", seq_len(K)), each = per_type))

rate <- matrix(stats::rgamma(n_bg * 1, shape = 0.5, rate = 12), n_bg, 1)
bg <- matrix(rep(rate, n_cells), n_bg, n_cells)

## BIMODAL: on everywhere, but ~6x higher in one type.  Both levels sit well
## above the gate, so the dip test sees a two-mode sample.
bi <- matrix(0.05, n_bimodal, n_cells)
for (i in seq_len(n_bimodal))
  bi[i, type == levels(type)[1 + (i - 1) %% K]] <- 0.30

## ONOFF: present in one type, ~0 elsewhere.  Above the gate only in that type,
## so the sample stage I actually tests is unimodal.
oo <- matrix(0.0005, n_onoff, n_cells)
for (i in seq_len(n_onoff))
  oo[i, type == levels(type)[1 + (i - 1) %% K]] <- 0.45

lambda <- rbind(bi, oo, bg)
rownames(lambda) <- c(sprintf("BIMODAL%02d", seq_len(n_bimodal)),
                      sprintf("ONOFF%02d",   seq_len(n_onoff)),
                      sprintf("BG%03d",      seq_len(n_bg)))
colnames(lambda) <- sprintf("cell%03d", seq_len(n_cells))

depth <- round(stats::runif(n_cells, 3000, 9000))
p <- sweep(lambda, 2, colSums(lambda), "/")
counts <- matrix(stats::rpois(length(p), sweep(p, 2, depth, "*")),
                 nrow(p), ncol(p), dimnames = dimnames(p))

## log-normalise exactly as the pipeline does: log(1 + 1e4 * c / colSum)
logmat <- log1p(sweep(counts, 2, pmax(colSums(counts), 1), "/") * 1e4)

truth <- setNames(as.character(type), colnames(counts))
demo <- list(counts = counts, logcounts = logmat, truth = truth,
             k = K,
             bimodal = grep("^BIMODAL", rownames(lambda), value = TRUE),
             onoff   = grep("^ONOFF",   rownames(lambda), value = TRUE),
             background = grep("^BG",   rownames(lambda), value = TRUE),
             note = paste("Simulated for demonstration only; not a dataset from",
                          "the manuscript. Built by demo/make_demo_data.R."))
saveRDS(demo, "demo/demo_data.rds")
cat(sprintf("wrote demo/demo_data.rds: %d genes x %d cells, %d types\n",
            nrow(counts), ncol(counts), K))
cat(sprintf("library size %d-%d, median genes detected per cell %d\n",
            min(colSums(counts)), max(colSums(counts)),
            as.integer(stats::median(colSums(counts > 0)))))
