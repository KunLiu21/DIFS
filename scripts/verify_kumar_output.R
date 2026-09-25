# Verify a completed example without rerunning selection or clustering.
source("R/difs.R")
args <- commandArgs(TRUE)
if (length(args)>1L) stop("Use Rscript scripts/verify_kumar_output.R output-directory")
out <- if (length(args)) args[1] else "example/output"
files <- c("difs_fit.rds","selected_genes.csv","cell_labels.csv","summary.csv",
           "sessionInfo.txt","input_checksums.csv","historical_comparison.csv","run.log")
if (!all(file.exists(file.path(out,files)))) stop("Some required output files are absent")
fit <- readRDS(file.path(out,"difs_fit.rds"))
genes <- read.csv(file.path(out,"selected_genes.csv"),stringsAsFactors=FALSE)
labs <- read.csv(file.path(out,"cell_labels.csv"),stringsAsFactors=FALSE,colClasses="character")
s <- read.csv(file.path(out,"summary.csv"),stringsAsFactors=FALSE)
ref <- read.csv("example/reference/benchmark_Kumar_n100_true_seed1.csv")
stopifnot(nrow(s)==1L,nrow(labs)==246L,!anyNA(labs),!anyDuplicated(labs$cell),
          all(nzchar(labs$predicted)),all(nzchar(labs$preliminary)),
          identical(genes$gene,fit$features),isTRUE(fit$stage2_executed),
          identical(names(fit$labels),labs$cell),identical(as.character(fit$labels),labs$predicted),
          identical(as.character(fit$preliminary_labels[labs$cell]),labs$preliminary),
          fit$sc3_audit$n_missing_after==0L,fit$sc3_audit$n_total==246L,
          fit$gate_rule=="submitted",fit$k_input==3L,fit$seed==1L,
          fit$n_features_requested==100L,nrow(genes)==s$n_features)
m <- difs_metrics(setNames(labs$truth,labs$cell),setNames(labs$predicted,labs$cell))
stopifnot(all(abs(m-as.numeric(unlist(s[1,names(m)])))<1e-10))
delta <- m-as.numeric(unlist(ref[1,names(m)]))
print(data.frame(metric=names(m),observed=unname(m),difference_from_historical=unname(delta)))
if (any(abs(delta)>1e-8) || s$n_features!=ref$n_features_realised || s$gate_n!=ref$gate_n)
  stop("Pipeline outputs are coherent, but historical numerical agreement is not established; inspect data and environment")
cat("PASS: output identities, full stage-II execution, all 246 labels, four metrics and historical counts agree.\n")
cat("This verifies this worked configuration, not an independent rerun of the full benchmark.\n")
