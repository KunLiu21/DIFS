
source("Functions_controlled.R")

args <- commandArgs(trailingOnly = TRUE)
for (a in args) {
  kv <- regmatches(a, regexec("^([A-Za-z._]+)=(.*)$", a))[[1]]
  if (length(kv) == 3L) assign(kv[2], kv[3])
}
if (!exists("group", inherits = FALSE)) group <- "A"

DATASETS_A <- c(
  "Kumar",          #  246 cells,  3 classes  (random genes reach ARI 0.99)
  "Trapnell",       #  222 cells,  3 classes  (hard for every method)
  "SimKumar4easy",  #  500 cells,  4 classes
  "SimKumar4hard",  #  499 cells,  4 classes
  "Zhengmix4eq",    # 3994 cells,  4 classes  10x
  "Romanov",        # 2881 cells,  7 classes  FEAST used this
  "Zhengmix8eq",    # 3994 cells,  8 classes  10x
  "Lawlor",         #  638 cells,  8 classes  pancreas, Fluidigm C1
  "Koh",            #  531 cells,  9 classes
  "Darmanis",       #  466 cells,  9 classes
  "Muraro",         # 2640 cells, 11 classes  pancreas, CEL-Seq2
  "Fletcher")       #  616 cells, 13 classes

DATASETS_B <- c(
  "Baron",          # 8569 cells, 14 classes  pancreas, inDrop -- R1 asked for it
  "Zhengmix4uneq")  # 6497 cells,  4 classes  10x, unbalanced

DATASETS <- if (identical(group, "B")) DATASETS_B else DATASETS_A

METHODS <- c("monocle", "TSCAN", "Refined Louvain", "SLM", "SC3")

FEATURE_METHODS <- c(
  "DIFS",            # as published: stage I ranked by the MC normal p-value
  "DIFS_hartigan",   # SAME pipeline, stage I ranked by Hartigan's p-value.
  "GateOnly",        # the gate, then NO ranking.  DIFS - GateOnly is what the
  "Seurat",
  "FEAST",
  "SC3",             # SC3's own filter, whole set (NOT 1000 -- record it)
  "SC3_random1000",  # the submitted SC3 baseline, kept as a random control
  "Variance",        # plain variance: the strongest competitor in the pilot
  "Random")

K_POLICIES <- c("estimated")

N_POLICIES <- c("original")

SEEDS <- c(1, 2, 3)

MIN_EXPRESSION <- log(5)

data.paths <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
nm <- gsub(".*/([^.]*).*", "\\1", data.paths)
keep <- nm %in% DATASETS
missing <- setdiff(DATASETS, nm)
if (length(missing))
  warning("not found in ../source: ", paste(missing, collapse = ", "),
          "\n  present: ", paste(sort(nm), collapse = ", "), call. = FALSE)
if (!any(keep)) stop("none of the requested datasets are in ../source")
data.paths <- data.paths[keep]

scenarios <- expand.grid(
  method                   = METHODS,
  feature.selection.method = FEATURE_METHODS,
  k_policy                 = K_POLICIES,
  n_policy                 = N_POLICIES,
  seed                     = SEEDS,
  data.path                = data.paths,
  stringsAsFactors         = FALSE)
scenarios$min.expression <- MIN_EXPRESSION

prepDir("../scenarios")
f <- paste0("../scenarios/scenarios_controlled_", group, ".Rdata")
save(scenarios, file = f)
save(scenarios, file = "../scenarios/scenarios_controlled.Rdata")

cat("\ngroup     : ", group, "\n", sep = "")
cat("datasets  : ", paste(gsub(".*/([^.]*).*", "\\1", data.paths), collapse = ", "),
    "\n", sep = "")
cat("scenarios : ", nrow(scenarios), "\n", sep = "")
cat("saved to  : ", f, " (and scenarios_controlled.Rdata)\n", sep = "")
cat("\n>>> set the SLURM array to  1-", nrow(scenarios), "%50\n", sep = "")
cat(">>> and check the walltime ends before 2026-09-18 08:00\n\n")
print(table(scenarios$feature.selection.method,
            gsub(".*/([^.]*).*", "\\1", scenarios$data.path)))
