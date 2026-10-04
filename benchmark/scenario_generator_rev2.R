# scenario_generator_rev2.R: revision-2 analysis source.
# Source: code/as_run/06_revision2_reps_krange_abl2_runtime/scenario_generator_rev2.R
# Usage, scope and limitations: benchmark/REVISION2.md.
source("Functions_controlled.R")

DATASETS <- c("Kumar", "Trapnell", "SimKumar4easy", "SimKumar4hard", "Zhengmix4eq",
              "Romanov", "Zhengmix8eq", "Lawlor", "Koh", "Darmanis", "Muraro",
              "Fletcher", "Baron")                    # rough order of cost, Baron last
paths <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
nm    <- gsub(".*/([^.]*).*", "\\1", paths)
miss  <- setdiff(DATASETS, nm)
if (length(miss)) stop("not in ../source: ", paste(miss, collapse = ", "))
path_of <- setNames(paths[match(DATASETS, nm)], DATASETS)

base_row <- function(...) data.frame(..., stringsAsFactors = FALSE)
order_rows <- function(d) {
  d <- d[order(match(gsub(".*/([^.]*).*", "\\1", d$data.path), DATASETS),
               d$method, d$n_features), ]
  rownames(d) <- NULL; d
}
report <- function(d, set) {
  save_path <- sprintf("../scenarios/scenarios_%s.Rdata", set)
  scenarios <- d; save(scenarios, file = save_path)
  is_baron <- grepl("/Baron\\.", d$data.path) | grepl("/Baron$", sub("\\.rds$", "", d$data.path))
  b <- which(is_baron); o <- which(!is_baron)
  cat(sprintf("\n=== %s: %d scenarios -> %s\n", set, nrow(d), save_path))
  cat(sprintf("    non-Baron rows %d-%d (%d)   Baron rows %s (%d)\n",
              min(o), max(o), length(o),
              if (length(b)) sprintf("%d-%d", min(b), max(b)) else "none", length(b)))
  if (length(b) && any(diff(b) != 1)) stop("Baron rows are not contiguous")
  invisible(d)
}
prepDir("../scenarios")

ARMS <- rbind(
  base_row(arm = "HVG_stage1only",  feature.selection.method = "DIFS_stage1only",
           n_policy = "original",     stage1_source = "hvg", final_set = "mixed"),
  base_row(arm = "DIFS_fixed",      feature.selection.method = "DIFS",
           n_policy = "fixed_all",    stage1_source = "dip", final_set = "mixed"),
  base_row(arm = "DIFS_fixratio",   feature.selection.method = "DIFS",
           n_policy = "fixed_stage1", stage1_source = "dip", final_set = "mixed"),
  base_row(arm = "HVG_fixed",       feature.selection.method = "DIFS",
           n_policy = "fixed_all",    stage1_source = "hvg", final_set = "mixed"),
  base_row(arm = "HVG_full",        feature.selection.method = "DIFS",
           n_policy = "original",     stage1_source = "hvg", final_set = "mixed"),
  base_row(arm = "DIFS_stage2only", feature.selection.method = "DIFS",
           n_policy = "original",     stage1_source = "dip", final_set = "stage2only"))
g <- expand.grid(a = seq_len(nrow(ARMS)), n_features = c(100, 300, 1000),
                 d = DATASETS, stringsAsFactors = FALSE)
abl2 <- cbind(ARMS[g$a, ], base_row(
  n_features = g$n_features, data.path = path_of[g$d], method = "Refined Louvain",
  k_policy = "true", seed = 1, min.expression = log(5), gate_rule = "submitted",
  cluster_seed = NA_integer_, subsample_rep = 0L, subsample_frac = 1, k_offset = NA_integer_))
report(order_rows(abl2), "abl2")

fam <- rbind(
  base_row(cluster_seed = 101:105, subsample_rep = 0L,  subsample_frac = 1,   seed = 101:105),
  base_row(cluster_seed = NA_integer_, subsample_rep = 1:5, subsample_frac = 0.9, seed = 1))
g <- expand.grid(f = seq_len(nrow(fam)), fm = c("DIFS", "Seurat", "FEAST"),
                 n_features = c(100, 300), method = c("Refined Louvain", "SC3"),
                 d = DATASETS, stringsAsFactors = FALSE)
reps <- cbind(base_row(arm = g$fm, feature.selection.method = g$fm, n_policy = "original",
                       stage1_source = "dip", final_set = "mixed",
                       n_features = g$n_features, data.path = path_of[g$d], method = g$method,
                       k_policy = "true", min.expression = log(5), gate_rule = "submitted",
                       k_offset = NA_integer_), fam[g$f, ])
report(order_rows(reps), "reps")

g <- expand.grid(off = c(-2L, -1L, 1L, 2L), fm = c("DIFS", "Seurat", "FEAST"),
                 method = c("Refined Louvain", "SC3"), d = DATASETS, stringsAsFactors = FALSE)
krange <- base_row(arm = g$fm, feature.selection.method = g$fm, n_policy = "original",
                   stage1_source = "dip", final_set = "mixed", n_features = 300,
                   data.path = path_of[g$d], method = g$method, k_policy = "true", seed = 1,
                   min.expression = log(5), gate_rule = "submitted",
                   cluster_seed = NA_integer_, subsample_rep = 0L, subsample_frac = 1,
                   k_offset = g$off)
report(order_rows(krange), "krange")

cat("\nThe stored runs supply the rest of each comparison: GateOnly, DIFS_stage1only and\n",
    "DIFS for abl2 (results/abl, Baron repaired); seed 1 at the historical seeds for\n",
    "reps and the k_offset = 0 point for krange (results/ncurve, Baron repaired).\n", sep = "")
