
source("Functions_controlled.R")

DATASETS <- c("Kumar", "Trapnell", "SimKumar4easy", "SimKumar4hard",
              "Zhengmix4eq", "Romanov", "Zhengmix8eq", "Lawlor",
              "Koh", "Darmanis", "Muraro", "Fletcher", "Baron")

METHODS <- c("Refined Louvain", "SC3")

FEATURE_METHODS <- c("DIFS",       # the method
                     "GateOnly",   # gate, no ranking -> isolates the ranking
                     "Variance",   # kept as the internal control it has always been
                     "Seurat",
                     "FEAST")

N_FEATURES <- c(100, 200, 300, 500, 750, 1000, 1500, 2500, 0)

K_POLICIES <- c("estimated", "true")
N_POLICIES <- "original"    # n comes from N_FEATURES, not from the policy
SEEDS      <- 1             # deterministic feature selection; ARI_sd was ~0
MIN_EXPRESSION <- log(5)

data.paths <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
nm   <- gsub(".*/([^.]*).*", "\\1", data.paths)
keep <- nm %in% DATASETS
miss <- setdiff(DATASETS, nm)
if (length(miss)) warning("not in ../source: ", paste(miss, collapse = ", "),
                          call. = FALSE)
if (!any(keep)) stop("no requested dataset found in ../source")
data.paths <- data.paths[keep]

scenarios <- expand.grid(
  method                   = METHODS,
  feature.selection.method = FEATURE_METHODS,
  n_features               = N_FEATURES,
  k_policy                 = K_POLICIES,
  n_policy                 = N_POLICIES,
  seed                     = SEEDS,
  data.path                = data.paths,
  stringsAsFactors         = FALSE)
scenarios$min.expression <- MIN_EXPRESSION

prepDir("../scenarios")
save(scenarios, file = "../scenarios/scenarios_ncurve.Rdata")

cat("\nk policy  : ", paste(K_POLICIES, collapse = ", "),
    "   (results are written with k-<policy> in the file name, so this run\n",
    "             does NOT overwrite the k = estimated results already on disk)\n",
    sep = "")
cat("datasets  : ", paste(gsub(".*/([^.]*).*", "\\1", data.paths), collapse = ", "),
    "\n", sep = "")
cat("n grid    : ", paste(ifelse(N_FEATURES == 0, "gate_n", N_FEATURES),
                          collapse = ", "), "\n", sep = "")
cat("scenarios : ", nrow(scenarios), "\n", sep = "")
cat("\n>>> do NOT submit 1-", nrow(scenarios), ": most of it is already on disk.\n",
    ">>> run this instead, and submit the array list it prints:\n",
    ">>>   singularity exec $DIFS_SIF Rscript difs_missing_tasks.R set=ncurve\n",
    sep = "")
cat(">>> walltime must end before 2026-09-18 08:00\n\n")
print(table(scenarios$feature.selection.method, scenarios$n_features))
