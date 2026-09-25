
source("Functions_controlled.R")

DATASET    <- "Baron"
BUDGETS    <- c(100, 200, 300, 500, 750, 1000, 1500, 2500, 0)   # 0 = whole gate pool
K_BOTH     <- c("true", "estimated")
GATE_RULE  <- "submitted"        # must match the main grid or nothing pairs
MIN_EXPR   <- log(5)

data.paths <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
nm <- gsub(".*/([^.]*).*", "\\1", data.paths)
if (!DATASET %in% nm) stop("Baron.rds not found in ../source")
BARON <- data.paths[nm == DATASET]

mk <- function(fm, cl, n, kp) expand.grid(
  method                   = cl,
  feature.selection.method = fm,
  n_features               = n,
  k_policy                 = kp,
  n_policy                 = "original",
  seed                     = 1,
  data.path                = BARON,
  stringsAsFactors         = FALSE)

A <- mk(c("DIFS", "GateOnly", "Variance", "Seurat", "FEAST"), "SC3", BUDGETS, K_BOTH)
B <- mk("DIFS", "Refined Louvain", BUDGETS, K_BOTH)
C <- mk("DIFS", "Refined Louvain", c(100, 300, 1000), "true")
D <- mk("DIFS_hartigan", c("Refined Louvain", "SC3"), c(100, 300), "true")

C_in_B <- nrow(merge(C, B))
scenarios <- unique(rbind(A, B, D))
scenarios$min.expression <- MIN_EXPR
scenarios$gate_rule      <- GATE_RULE

prepDir("../scenarios")
save(scenarios, file = "../scenarios/scenarios_baronfix.Rdata")

cat("\nBaron rerun with the SC3 repair\n")
cat("group A  SC3 final, 5 methods x 9 budgets x 2 k    : ", nrow(A), "\n", sep = "")
cat("group B  DIFS x Louvain, 9 budgets x 2 k           : ", nrow(B), "\n", sep = "")
cat("group C  ablation DIFS cells                       : ", nrow(C),
    "  (", C_in_B, " already covered by B)\n", sep = "")
cat("group D  DIFS_hartigan, 2 budgets x 2 clusterers   : ", nrow(D), "\n", sep = "")
cat("--------------------------------------------------\n")
cat("unique scenarios to run                            : ", nrow(scenarios), "\n", sep = "")
cat("gate rule: ", GATE_RULE, "   seed: 1   n_policy: original\n", sep = "")
print(table(scenarios$feature.selection.method, scenarios$method))
cat("\n>>> mkdir -p $DIFS_ROOT/Rout/baronfix $DIFS_ROOT/results/baronfix\n")
cat(">>> sbatch --array=1-", nrow(scenarios), "%28 run_baronfix.sh\n", sep = "")
