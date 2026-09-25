
difs_args <- function(args, known, env = parent.frame()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    if (!kv[2] %in% known)
      stop("unknown argument '", kv[2], "'. This script accepts: ",
           paste(known, collapse = ", "))
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
  invisible(NULL)
}

args <- commandArgs(TRUE)
set <- "ncurve"                      # which scenario table; see below
sc3_repair <- 1
difs_args(args, c("set", "k", "sc3_repair"))

source("Functions_controlled.R")

if (isTRUE(as.numeric(sc3_repair) == 1)) {
  if (!file.exists("difs_sc3_complete.R"))
    stop("difs_sc3_complete.R must be in the code directory, or pass sc3_repair=0 ",
         "to deliberately reproduce the defective behaviour")
  source("difs_sc3_complete.R")
  cat("SC3 missing-label repair: ON (difs_sc3_complete.R). Results on datasets\n",
      "  above 5000 cells will differ from anything stored before 2026-09-20.\n", sep = "")
} else {
  cat("*** SC3 missing-label repair: OFF (sc3_repair=0). On any dataset above\n",
      "*** 5000 cells this reproduces the DEFECT: metrics computed on 5000 cells\n",
      "*** with the remainder unassigned. Do not use for new results.\n", sep = "")
}

resolution <- set
load(sprintf("../scenarios/scenarios_%s.Rdata", set))

parameters <- scenarios[k, ]
task_id <- k
rm(k)

gate_rule <- if ("gate_rule" %in% names(parameters))
  as.character(parameters$gate_rule) else "submitted"

data.name <- gsub(".*/([^.]*).*", "\\1", as.character(parameters$data.path))

results <- Run.method.comparison.controlled(
  data.path                = as.character(parameters$data.path),
  clustering.method        = as.character(parameters$method),
  feature.selection.method = as.character(parameters$feature.selection.method),
  k_policy                 = as.character(parameters$k_policy),
  n_policy                 = as.character(parameters$n_policy),
  n_fixed                  = parameters$n_features,
  gate_rule                = gate_rule,
  min.expression           = parameters$min.expression,
  seed                     = parameters$seed)

results$n_requested <- parameters$n_features   # 0 = gate_n

out <- prepDir(paste0("../results/", resolution))
saveRDS(results, file.path(out, sprintf(
  "%s_%s_%s_n%04d_k-%s_g-%s_seed-%s.rds",
  data.name, parameters$feature.selection.method, parameters$method,
  parameters$n_features, parameters$k_policy, gate_rule, parameters$seed)))

cat("task", task_id, "done:", data.name,
    as.character(parameters$feature.selection.method),
    as.character(parameters$method),
    "| k", as.character(parameters$k_policy),
    "| gate", gate_rule,
    "| requested", parameters$n_features,
    "| used", results$n_features,
    "| gate_n", results$gate_n,
    "| clusters", results$n_clusters_found, "/", results$k_used_clustering,
    "| dropped", results$n_cells_dropped,
    "| ARI =", round(results$ARI, 4), "\n")
