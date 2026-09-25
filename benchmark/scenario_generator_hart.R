
source("Functions_controlled.R")

DATASETS <- c("Kumar", "Trapnell", "SimKumar4easy", "SimKumar4hard",
              "Zhengmix4eq", "Romanov", "Zhengmix8eq", "Lawlor",
              "Koh", "Darmanis", "Muraro", "Fletcher", "Baron")

METHODS         <- c("Refined Louvain", "SC3")
FEATURE_METHODS <- "DIFS_hartigan"    # DIFS comes from the stored ncurve runs
N_FEATURES      <- c(100, 300)
K_POLICIES      <- "true"
N_POLICIES      <- "original"
SEEDS           <- 1
GATE_RULE       <- "submitted"
MIN_EXPRESSION  <- log(5)

data.paths <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
nm <- gsub(".*/([^.]*).*", "\\1", data.paths)
missing <- setdiff(DATASETS, nm)
if (length(missing))
  stop("missing from ../source: ", paste(missing, collapse = ", "))
data.paths <- data.paths[nm %in% DATASETS]

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
scenarios$gate_rule      <- GATE_RULE

prepDir("../scenarios")
save(scenarios, file = "../scenarios/scenarios_hart.Rdata")

cat("\nReference-distribution check at small budgets\n")
cat("datasets  : ", length(DATASETS), "\n", sep = "")
cat("arm       : ", FEATURE_METHODS, "  (DIFS side comes from ../results/ncurve)\n", sep = "")
cat("budgets   : ", paste(N_FEATURES, collapse = ", "), "\n", sep = "")
cat("k policy  : ", K_POLICIES, "   gate rule: ", GATE_RULE, "\n", sep = "")
cat("scenarios : ", nrow(scenarios),
    "   -> ../scenarios/scenarios_hart.Rdata\n", sep = "")
print(table(scenarios$n_features, scenarios$method))
cat("\n>>> mkdir -p $DIFS_ROOT/Rout/hart $DIFS_ROOT/results/hart\n")
cat(">>> sbatch --array=1-", nrow(scenarios), "%26 run_hart.sh\n", sep = "")
cat("\nAfter it finishes:\n")
cat(">>> singularity exec $DIFS_SIF Rscript difs_hart_compare.R\n\n")
