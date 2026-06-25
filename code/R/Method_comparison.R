args=(commandArgs(TRUE))
for(i in 1:length(args)){
  eval(parse(text=args[i]))
}

source("Functions_mixture_auto.R")
set.seed(2)
resolution = "clustering_mixtureautofeatures_seeds_2"

load("../scenarios/scenarios_parameters.Rdata")

parameters <- scenarios[k, ]
rm(k)

data.name = gsub(".*/([^.]*).*", "\\1", as.character(parameters$data.path))


results <- Run.method.comparison(data.path = as.character(parameters$data.path), 
 clustering.method = parameters$method, feature.selection.method = parameters$feature.selection.method)

prepDir(paste("../results/", resolution, sep =""))
saveRDS(results, paste("../results/", resolution, "/",data.name,"_",parameters$feature.selection.method, "_clustering.method_",
                parameters$method,".rds", sep = ""))




rm(list = ls())