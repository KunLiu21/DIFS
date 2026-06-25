
data.paths <- list.files("../source", full.names = T)
source("Functions_mixture_auto.R")
methods = c("monocle", "TSCAN","Refined Louvain","SLM", "SC3")
feature.selection.method = c("monocle", "SC3", "Seurat", "Feast", "SCT", "DIFs")

scenarios = expand.grid(method = methods, feature.selection.method = feature.selection.method, data.path = data.paths)
prepDir("../scenarios")
save(scenarios, file="../scenarios/scenarios_parameters.Rdata")

