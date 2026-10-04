# scenario_generator_sc3fix.R: revision-2 analysis source.
# Source: code/as_run/05_repairs_baron_sc3_hartigan/scenario_generator_sc3fix.R
# Usage, scope and limitations: benchmark/REVISION2.md.
source("Functions_controlled.R")

DATASETS   <- "Baron"
K_POLICIES <- "true"
N_POLICIES <- "original"
SEEDS      <- 1
GATE_RULE  <- "submitted"
MIN_EXPRESSION <- log(5)

N_FEATURES <- c(100, 1000)
ARM_SC3     <- c("DIFS", "GateOnly", "Variance", "Seurat", "FEAST")
ARM_LOUVAIN <- "DIFS"

data.paths <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
nm   <- gsub(".*/([^.]*).*", "\\1", data.paths)
if (!any(nm %in% DATASETS))
  stop("Baron.rds not found in ../source")
data.paths <- data.paths[nm %in% DATASETS]

mk <- function(methods, clusterer) expand.grid(
  method                   = clusterer,
  feature.selection.method = methods,
  n_features               = N_FEATURES,
  k_policy                 = K_POLICIES,
  n_policy                 = N_POLICIES,
  seed                     = SEEDS,
  data.path                = data.paths,
  stringsAsFactors         = FALSE)

scenarios <- rbind(mk(ARM_SC3, "SC3"), mk(ARM_LOUVAIN, "Refined Louvain"))
scenarios$min.expression <- MIN_EXPRESSION
scenarios$gate_rule      <- GATE_RULE

prepDir("../scenarios")
save(scenarios, file = "../scenarios/scenarios_sc3fix.Rdata")

cat("\nSC3 repair sensitivity check\n")
cat("dataset   : ", DATASETS, "\n", sep = "")
cat("gate rule : ", GATE_RULE, "   k policy: ", K_POLICIES, "\n", sep = "")
cat("budgets   : ", paste(N_FEATURES, collapse = ", "), "\n", sep = "")
cat("scenarios : ", nrow(scenarios),
    "   -> ../scenarios/scenarios_sc3fix.Rdata\n", sep = "")
print(table(scenarios$feature.selection.method, scenarios$method))
cat("\n>>> mkdir -p $DIFS_ROOT/Rout/sc3fix $DIFS_ROOT/results/sc3fix\n")
cat(">>> sbatch --array=1-", nrow(scenarios), "%12 run_sc3fix.sh\n", sep = "")
cat("\nAfter it finishes:\n")
cat(">>> singularity exec $DIFS_SIF Rscript difs_sc3fix_compare.R\n\n")
