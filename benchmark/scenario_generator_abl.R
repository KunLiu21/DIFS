
source("Functions_controlled.R")

DATASETS <- c("Kumar", "Trapnell", "SimKumar4easy", "SimKumar4hard",
              "Zhengmix4eq", "Romanov", "Zhengmix8eq", "Lawlor",
              "Koh", "Darmanis", "Muraro", "Fletcher", "Baron")

ARMS       <- c("GateOnly", "DIFS_stage1only", "DIFS")
N_FEATURES <- c(100, 300, 1000)
METHODS    <- "Refined Louvain"
K_POLICIES <- "true"
N_POLICIES <- "original"
SEEDS      <- 1
MIN_EXPRESSION <- log(5)

data.paths <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
nm   <- gsub(".*/([^.]*).*", "\\1", data.paths)
miss <- setdiff(DATASETS, nm)
if (length(miss)) warning("not in ../source: ", paste(miss, collapse = ", "),
                          call. = FALSE)
data.paths <- data.paths[nm %in% DATASETS]
if (!length(data.paths)) stop("no requested dataset found in ../source")

scenarios <- expand.grid(
  method                   = METHODS,
  feature.selection.method = ARMS,
  n_features               = N_FEATURES,
  k_policy                 = K_POLICIES,
  n_policy                 = N_POLICIES,
  seed                     = SEEDS,
  data.path                = data.paths,
  stringsAsFactors         = FALSE)
scenarios$min.expression <- MIN_EXPRESSION

prepDir("../scenarios")
save(scenarios, file = "../scenarios/scenarios_abl.Rdata")

cat("\ndatasets  : ", paste(sort(gsub(".*/([^.]*).*", "\\1", data.paths)),
                            collapse = ", "), "\n", sep = "")
cat("arms      : ", paste(ARMS, collapse = ", "), "\n", sep = "")
cat("n grid    : ", paste(N_FEATURES, collapse = ", "), "\n", sep = "")
cat("scenarios : ", nrow(scenarios), "   -> ../scenarios/scenarios_abl.Rdata\n",
    sep = "")
cat("\n>>> mkdir -p $DIFS_ROOT/Rout/abl $DIFS_ROOT/results/abl\n")
cat(">>> sbatch --array=1-", nrow(scenarios), "%40 run_abl.sh\n", sep = "")
cat(">>> walltime must end before 2026-09-18 08:00\n\n")
print(table(scenarios$feature.selection.method, scenarios$n_features))
