
args <- commandArgs(TRUE)
for (i in seq_along(args)) eval(parse(text = args[i]))

source("Functions_controlled.R")

resolution <- "controlled_comparison"
load("../scenarios/scenarios_controlled.Rdata")

parameters <- scenarios[k, ]
task_id <- k
rm(k)

data.name <- gsub(".*/([^.]*).*", "\\1", as.character(parameters$data.path))

results <- Run.method.comparison.controlled(
  data.path                = as.character(parameters$data.path),
  clustering.method        = as.character(parameters$method),
  feature.selection.method = as.character(parameters$feature.selection.method),
  k_policy                 = as.character(parameters$k_policy),
  n_policy                 = as.character(parameters$n_policy),
  min.expression           = parameters$min.expression,
  seed                     = parameters$seed)

out <- prepDir(paste0("../results/", resolution))
saveRDS(results, file.path(out, sprintf(
  "%s_%s_%s_k-%s_n-%s_seed-%s.rds",
  data.name, parameters$feature.selection.method, parameters$method,
  parameters$k_policy, parameters$n_policy, parameters$seed)))

cat("task", task_id, "done: ARI =", results$ARI,
    "| features =", results$n_features, "\n")
