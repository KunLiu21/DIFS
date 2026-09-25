required <- c("Seurat","SeuratObject","diptest","SC3","FEAST","SingleCellExperiment",
              "SummarizedExperiment","S4Vectors","matrixStats","magrittr","DuoClustering2018")
versions <- vapply(required,function(p) tryCatch(as.character(packageVersion(p)),error=function(e) NA_character_),character(1))
print(data.frame(package=required,version=versions),row.names=FALSE)
print(sessionInfo())
if (anyNA(versions)) quit(status=1)
