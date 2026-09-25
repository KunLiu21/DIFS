source("R/difs.R")
input <- "example/input/Kumar.rds"; out <- "example/output"
for (arg in commandArgs(TRUE)) {
  kv <- strsplit(arg,"=",fixed=TRUE)[[1]]
  if (length(kv)!=2L || !kv[1] %in% c("input","out")) stop("Use input=... and/or out=...")
  assign(kv[1],kv[2])
}
if (!file.exists(input)) stop("Input missing. Run example/prepare_kumar.R or supply input=path/to/Kumar.rds")
if (dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE))) stop("Output is not empty; choose a new out= directory")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
log <- file(file.path(out,"run.log"),open="wt"); sink(log,split=TRUE); sink(log,type="message")
difs_require(c("Seurat","SeuratObject"))
seu <- readRDS(input)
if (!inherits(seu,"Seurat")) stop("Input must be a prepared Kumar Seurat object")
counts <- if (utils::packageVersion("SeuratObject") >= "5.0.0") SeuratObject::GetAssayData(seu,assay="RNA",layer="counts") else SeuratObject::GetAssayData(seu,assay="RNA",slot="counts")
truth <- setNames(as.character(seu$trueclass),colnames(seu))
if (ncol(counts)!=246L || anyNA(truth) || length(unique(truth))!=3L) stop("Expected Kumar: 246 cells and 3 annotated classes")
cat("Complete Kumar DIFS example: 100 requested features, seed=1, true-k=3 (oracle).\n")
fit <- difs_fit(counts,k=3L,n_features=100L,seed=1L)
metrics <- difs_metrics(truth,fit$labels)
if (!isTRUE(fit$stage2_executed) || any(!is.finite(metrics)) || length(fit$preliminary_labels)!=246L || anyNA(fit$preliminary_labels)) stop("Full example validation failed")
saveRDS(fit,file.path(out,"difs_fit.rds"))
write.csv(data.frame(gene=fit$features),file.path(out,"selected_genes.csv"),row.names=FALSE)
write.csv(data.frame(cell=names(fit$labels),truth=truth[names(fit$labels)],
  preliminary=as.character(fit$preliminary_labels[names(fit$labels)]),predicted=fit$labels),file.path(out,"cell_labels.csv"),row.names=FALSE)
summary <- data.frame(dataset="Kumar",k_policy="true",k=3,seed=1,n_requested=100,n_features=fit$n_features_realised,
  n_cells=length(fit$labels),stage1_n=fit$stage1_n,stage2_n=fit$stage2_n,ratio=fit$ratio,
  gate_n=fit$gate_n,stage2_executed=fit$stage2_executed,t(metrics),check.names=FALSE)
write.csv(summary,file.path(out,"summary.csv"),row.names=FALSE)
writeLines(capture.output(fit$session),file.path(out,"sessionInfo.txt"))
write.csv(data.frame(input_md5=unname(tools::md5sum(input)),null_table_md5=unname(tools::md5sum("inst/extdata/dip_null_tables.rds"))),file.path(out,"input_checksums.csv"),row.names=FALSE)
ref <- read.csv("example/reference/benchmark_Kumar_n100_true_seed1.csv")
comparison <- data.frame(metric=names(metrics),observed=unname(metrics),historical=as.numeric(unlist(ref[1,names(metrics)])))
comparison$difference <- comparison$observed-comparison$historical
write.csv(comparison,file.path(out,"historical_comparison.csv"),row.names=FALSE)
print(summary);print(comparison)
cat("FULL PIPELINE COMPLETED. Historical numerical reproduction still requires comparison of data, environment and outputs.\n")
sink(type="message");sink();close(log)
