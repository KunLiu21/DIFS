# difs_gene_rank.R: revision-2 analysis source.
# Source: code/as_run/07_feature_sets_gate_ranking/difs_gene_rank.R
# Usage, scope and limitations: benchmark/REVISION2.md.
difs_args <- function(args, env = parent.frame()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
}
args <- commandArgs(TRUE)
gene <- "RBM39"; ensembl <- "ENSG00000131051"; datasets <- ""
difs_args(args)

source("Functions_controlled.R")

find_gene <- function(rn, sym, ens) {
  hit <- which(rn == sym);                       if (length(hit)) return(c(hit[1], "symbol"))
  hit <- which(startsWith(rn, paste0(sym, "--")));if (length(hit)) return(c(hit[1], "symbol--chr"))
  hit <- which(rn == ens);                        if (length(hit)) return(c(hit[1], "ensembl"))
  hit <- which(sub("\\..*$", "", rn) == ens);     if (length(hit)) return(c(hit[1], "ensembl.version"))
  hit <- which(toupper(rn) == toupper(sym));      if (length(hit)) return(c(hit[1], "symbol(case)"))
  NULL
}

fs <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
if (nzchar(datasets)) {
  keep <- gsub(".*/([^.]*).*", "\\1", fs) %in% trimws(strsplit(datasets, ",")[[1]])
  fs <- fs[keep]
}

out <- NULL
for (f in fs) {
  nm <- gsub(".*/([^.]*).*", "\\1", f)
  seu <- tryCatch(readRDS(f), error = function(e) NULL); if (is.null(seu)) next
  lm  <- as.matrix(GetAssayData(seu, slot = "data"))
  rn  <- rownames(lm)
  h   <- find_gene(rn, gene, ensembl)
  if (is.null(h)) { cat(sprintf("%-16s %s not present\n", nm, gene)); next }
  idx <- as.integer(h[1]); how <- h[2]

  gate  <- difs_gate_genes(lm, log(5))
  ingate <- rn[idx] %in% gate
  nexpr <- sum(lm[idx, ] > log(5))

  v      <- apply(lm, 1, stats::var)
  pv_var <- 100 * (rank(-v, ties.method = "min")[idx] - 1) / (length(v) - 1)

  pv_vst <- NA_real_
  s2 <- tryCatch(FindVariableFeatures(seu, selection.method = "vst",
                                      nfeatures = nrow(seu), verbose = FALSE),
                 error = function(e) NULL)
  if (!is.null(s2)) {
    hv <- VariableFeatures(s2)
    p  <- match(rn[idx], hv)
    if (!is.na(p)) pv_vst <- 100 * (p - 1) / (length(hv) - 1)
  }

  pv_difs <- NA_real_
  if (ingate) {
    ranked  <- difs_stage1_ranking(lm, min.expression = log(5),
                                   key = "normal_lookup")
    p       <- match(rn[idx], ranked)
    if (!is.na(p)) pv_difs <- 100 * (p - 1) / (length(ranked) - 1)
  }

  cat(sprintf("%-16s matched by %-15s cells expressing %5d / %5d   in gate: %s\n",
              nm, how, nexpr, ncol(lm), ingate))
  cat(sprintf("%-16s   percentile  DIFS %s   variance %6.2f%%   vst %s\n", "",
              if (is.na(pv_difs)) "   n/a  " else sprintf("%6.2f%%", pv_difs),
              pv_var,
              if (is.na(pv_vst)) "   n/a  " else sprintf("%6.2f%%", pv_vst)))
  out <- rbind(out, data.frame(dataset = nm, matched_by = how,
                               n_expressing = nexpr, n_cells = ncol(lm),
                               in_gate = ingate, gate_n = length(gate),
                               pct_DIFS = pv_difs, pct_variance = pv_var,
                               pct_vst = pv_vst, stringsAsFactors = FALSE))
}

if (is.null(out)) { cat("\n", gene, " found in no dataset.\n", sep = ""); quit(save = "no") }
cat("\n================ ", gene, " ================\n", sep = "")
print(out, row.names = FALSE, digits = 4)
dir.create("../results/gene_rank", recursive = TRUE, showWarnings = FALSE)
write.csv(out, file.path("../results/gene_rank", paste0(gene, ".csv")), row.names = FALSE)

cat("\nHOW TO READ THIS (percentile: 0 = top of the ranking)\n")
cat("  variance LOW and DIFS HIGH  -> Figure 2's claim holds and generalises;\n")
cat("     this is the rebuttal, and it now rests on several datasets, not one.\n")
cat("  variance LOW and DIFS LOW   -> the modified dip test selects RBM39 too.\n")
cat("     Figure 2 is then a counterexample to our own method and must be\n")
cat("     REPLACED.  Better to find this now than to have a reviewer find it.\n")
cat("  in_gate FALSE               -> DIFS never tested the gene.  Report that\n")
cat("     as a gate effect, NOT as the dip test rejecting it.  Those are\n")
cat("     different claims and conflating them is exactly what Reviewer 1\n")
cat("     objected to in Section 3.1.3.\n")
