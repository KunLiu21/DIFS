# Method_comparison_sc3fix.R: revision-2 analysis source.
# Source: code/as_run/05_repairs_baron_sc3_hartigan/Method_comparison_sc3fix.R
# Usage, scope and limitations: benchmark/REVISION2.md.
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
set <- "sc3fix"                      # which scenario table; see below
difs_args(args, c("set", "k"))

source("Functions_controlled.R")

if (!file.exists("difs_sc3_complete.R"))
  stop("difs_sc3_complete.R must be in the code directory")
source("difs_sc3_complete.R")

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
