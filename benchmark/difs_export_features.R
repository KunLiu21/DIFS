# difs_export_features.R: revision-2 analysis source.
# Source: code/as_run/07_feature_sets_gate_ranking/difs_export_features.R
# Usage, scope and limitations: benchmark/REVISION2.md.
source("Functions_controlled.R")

difs_parse_args <- function(args, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) { message("ignoring argument: ", a); next }
    key <- kv[2L]; val <- sub('^"(.*)"$', "\\1", sub("^'(.*)'$", "\\1", kv[3L]))
    v <- if (grepl("^(TRUE|FALSE)$", val)) as.logical(val)
         else if (grepl("^-?[0-9.]+$", val)) as.numeric(val) else val
    assign(key, v, envir = env)
  }
}
difs_parse_args(commandArgs(TRUE))

if (!exists("only",    inherits = FALSE)) only    <- "Baron,Muraro,Lawlor,Koh,KohTCC"
if (!exists("out_dir", inherits = FALSE)) out_dir <- "../results/feature_sets"
if (!exists("n_fixed", inherits = FALSE)) n_fixed <- 1000
if (!exists("seed",    inherits = FALSE)) seed    <- 1
if (!exists("min.expression", inherits = FALSE)) min.expression <- log(5)

EXPORT_KEYS <- list(DIFS          = "normal_lookup",
                    DIFS_mc       = "mc_normal",
                    DIFS_hartigan = "hartigan")
FS <- c("DIFS", "DIFS_mc", "DIFS_hartigan", "GateOnly",
        "Seurat", "FEAST", "Variance", "Random")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
wanted <- trimws(strsplit(only, ",")[[1]])
paths  <- list.files("../source", pattern = "\\.rds$", full.names = TRUE)
nmv    <- gsub(".*/([^.]*).*", "\\1", paths)
paths  <- paths[nmv %in% wanted]
if (!length(paths)) stop("none of ", only, " found in ../source")

rows <- list(); sets <- list(); universe <- list(); idspace <- list()

norm_id <- function(x) {
  x <- sub("--.*$", "", x)          # Muraro's "GENE--chrN"
  x <- sub("\\.[0-9]+$", "", x)       # Ensembl version suffix
  toupper(x)
}
id_space <- function(x) {
  f <- function(p) mean(grepl(p, x))
  if (f("^ENS[A-Z]*G[0-9]{6,}") > 0.5) return("ensembl")
  if (f("^(NM_|NR_|XM_|XR_)[0-9]+")  > 0.5) return("refseq")
  if (f("^[0-9]+$")                  > 0.5) return("entrez")
  if (f("^[A-Za-z][A-Za-z0-9._-]*$") > 0.5) return("symbol")
  "other"
}

for (pth in paths) {
  dname <- gsub(".*/([^.]*).*", "\\1", pth)
  message("=== ", dname, " ===")
  data.seu     <- readRDS(pth)
  data.log.mat <- as.matrix(GetAssayData(data.seu, slot = "data"))
  data.mat     <- as.matrix(GetAssayData(data.seu, slot = "counts"))
  trueclass    <- data.seu$trueclass
  k_true       <- length(unique(trueclass))

  sce <- SingleCellExperiment(assays = list(counts = data.mat,
                                            logcounts = log2(data.mat + 1)))
  sce <- sc3_estimate_k(sce)
  k_est <- metadata(sce)$sc3$k_estimation
  cds <- suppressWarnings(sce_to_monocle(sce))
  cds <- suppressWarnings(estimateSizeFactors(cds))
  cds <- suppressWarnings(estimateDispersions(cds))

  universe[[dname]] <- norm_id(rownames(data.log.mat))
  idspace[[dname]]  <- id_space(rownames(data.log.mat))
  gate <- difs_gate_genes(data.log.mat, min.expression)
  message("  ", ncol(data.log.mat), " cells, ", nrow(data.log.mat), " genes, ",
          k_true, " classes, k_est ", k_est, ", gate passes ", length(gate))

  for (m in FS) {
    g <- tryCatch({
      if (m %in% names(EXPORT_KEYS)) {
        key <- EXPORT_KEYS[[m]]
        difs_features(data.log.mat = data.log.mat, data.seu = data.seu,
                      cluster_count = k_est, binary.bound = min.expression,
                      min.expression = min.expression,
                      tune_stage1 = TRUE, tune_ratio = TRUE,
                      stage1_key = key, seed = seed, n_fixed = n_fixed)$markers
      } else {
        r <- feature_ranking(m, data.seu = data.seu, data.log.mat = data.log.mat,
                             data.mat = data.mat, cds = cds,
                             cluster_count = k_est, seed = seed,
                             min.expression = min.expression)
        r[seq_len(min(n_fixed, length(r)))]
      }
    }, error = function(e) { message("  ", m, " FAILED: ", conditionMessage(e))
                             character(0) })
    sets[[paste(dname, m, sep = "|")]] <- norm_id(g)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = dname, method = m, n_cells = ncol(data.log.mat),
      n_genes_total = nrow(data.log.mat), gate_n = length(gate),
      k_true = k_true, k_estimated = k_est, n_selected = length(g),
      stringsAsFactors = FALSE)
    message("  ", m, ": ", length(g), " genes")
  }
}

saveRDS(sets, file.path(out_dir, "feature_sets.rds"))
info <- do.call(rbind, rows)
utils::write.csv(info, file.path(out_dir, "feature_sets_info.csv"), row.names = FALSE)

jacc <- function(a, b) length(intersect(a, b)) / length(union(a, b))
adj  <- function(a, b, P) {
  if (!length(a) || !length(b)) return(NA_real_)
  N <- mean(c(length(a), length(b)))
  ch <- N / (2 * P - N)
  (jacc(a, b) - ch) / (1 - ch)
}

dsets <- unique(info$dataset)
cat("\n========== identifier space per dataset ==========\n")
for (d in dsets)
  cat(sprintf("  %-14s %-8s %6d genes   e.g. %s\n", d, idspace[[d]],
              length(universe[[d]]),
              paste(utils::head(universe[[d]], 4), collapse = ", ")))
for (i in seq_along(dsets)) for (j in seq_along(dsets)) {
  if (j <= i) next
  if (!identical(idspace[[dsets[i]]], idspace[[dsets[j]]])) next
  P <- length(intersect(universe[[dsets[i]]], universe[[dsets[j]]]))
  if (P < 100)
    cat("  !! ", dsets[i], " and ", dsets[j], " are both '", idspace[[dsets[i]]],
        "' but share only ", P, " gene names.  The identifier detection is\n",
        "     wrong for at least one of them -- do NOT read their agreement.\n",
        sep = "")
}

cmp <- list()
for (m in FS) for (i in seq_along(dsets)) for (j in seq_along(dsets)) {
  if (j <= i) next
  d1 <- dsets[i]; d2 <- dsets[j]
  a <- sets[[paste(d1, m, sep = "|")]]
  b <- sets[[paste(d2, m, sep = "|")]]
  if (is.null(a) || is.null(b)) next
  comparable <- identical(idspace[[d1]], idspace[[d2]])
  P <- length(intersect(universe[[d1]], universe[[d2]]))
  cmp[[length(cmp) + 1L]] <- data.frame(
    method = m, d1 = d1, d2 = d2,
    id1 = idspace[[d1]], id2 = idspace[[d2]],
    shared_pool = P,
    shared = if (comparable) length(intersect(a, b)) else NA_integer_,
    jaccard  = if (comparable) round(jacc(a, b), 4) else NA_real_,
    adjusted = if (comparable && P > 0) round(adj(a, b, P), 4) else NA_real_,
    comparable = comparable, stringsAsFactors = FALSE)
}
ov <- do.call(rbind, cmp)
utils::write.csv(ov, file.path(out_dir, "feature_set_overlap.csv"), row.names = FALSE)

cat("\n========== selected set sizes ==========\n")
print(info, row.names = FALSE)
cat("\n========== cross-dataset agreement (chance-corrected) ==========\n")
print(ov[order(ov$method, -ov$adjusted), ], row.names = FALSE)
cat("\nHOW TO READ THIS\n")
cat("  Koh vs KohTCC is the SAME cells quantified twice: whatever agreement a\n")
cat("  method reaches there is its ceiling under a purely technical change.\n")
cat("  Baron/Muraro/Lawlor are the same cell types on three protocols.  If a\n")
cat("  method's pancreas agreement is close to its Koh-vs-KohTCC agreement,\n")
cat("  its genes track biology; if it collapses, they track protocol.\n")
cat("  Compare DIFS against Variance and Seurat -- the claim is only meaningful\n")
cat("  RELATIVE to what the baselines achieve on the same comparison.\n")
cat("\nwritten to: ", normalizePath(out_dir), "\n")
