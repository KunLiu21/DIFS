# difs_gate_scan.R: revision-2 analysis source.
# Source: code/as_run/07_feature_sets_gate_ranking/difs_gate_scan.R
# Usage, scope and limitations: benchmark/REVISION2.md.
source("Functions_controlled.R")

args <- commandArgs(TRUE)
for (a in args) {
  kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
  if (length(kv) == 3L) assign(kv[2],
    if (grepl("^-?[0-9.]+$", kv[3])) as.numeric(kv[3]) else kv[3], envir = globalenv())
}
if (!exists("only",    inherits = FALSE)) only    <- NULL
if (!exists("out_dir", inherits = FALSE)) out_dir <- "../results/gate_scan"
if (!exists("min_expression", inherits = FALSE)) min_expression <- log(5)
if (!exists("src_dir", inherits = FALSE)) src_dir <- "../source"

RULES <- list(
  "max(0.05n,35)  [submitted]" = function(n, k) max(0.05 * n, 35),
  "min(0.05n,35)"              = function(n, k) min(0.05 * n, 35),
  "0.05n"                      = function(n, k) 0.05 * n,
  "0.02n"                      = function(n, k) 0.02 * n,
  "35 fixed"                   = function(n, k) 35,
  "10 fixed"                   = function(n, k) 10,
  "min(0.05n, n/(2k))"         = function(n, k) min(0.05 * n, n / (2 * k)),
  "n/(4k)"                     = function(n, k) n / (4 * k)
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
paths <- sort(list.files(src_dir, pattern = "\\.rds$", full.names = TRUE))
nm <- gsub(".*/([^.]*).*", "\\1", paths)
if (!is.null(only)) paths <- paths[nm %in% trimws(strsplit(only, ",")[[1]])]
if (!length(paths)) stop("no dataset found in ", src_dir)

rows <- list()
for (p in paths) {
  dn <- gsub(".*/([^.]*).*", "\\1", p)
  seu <- tryCatch(readRDS(p), error = function(e) NULL)
  if (is.null(seu)) { message("skip ", dn); next }
  m  <- SeuratObject::GetAssayData(seu, slot = "data")
  n  <- ncol(m); G <- nrow(m)
  tc <- as.character(seu$trueclass)
  sz <- sort(table(tc))
  k  <- length(sz)
  n_expr <- Matrix::rowSums(m > min_expression)   # ONE pass; every rule cuts this
  message(sprintf("%-16s %6d cells  %6d genes  %2d classes  (smallest class %d)",
                  dn, n, G, k, as.integer(sz[1])))
  for (rn in names(RULES)) {
    thr <- RULES[[rn]](n, k)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = dn, n_cells = n, n_genes = G, k = k,
      min_class = as.integer(sz[1]), median_class = as.integer(stats::median(sz)),
      rule = rn, threshold = round(thr, 1),
      pct_cells = round(100 * thr / n, 2),
      gate_n = sum(n_expr >= thr),
      pct_genes = round(100 * sum(n_expr >= thr) / G, 2),
      classes_excluded = sum(sz < thr),
      stringsAsFactors = FALSE)
  }
  rm(m, n_expr); invisible(gc())
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(out_dir, "gate_scan.csv"), row.names = FALSE)

cat("\n================ gate_n by rule ================\n")
w <- reshape(res[, c("dataset", "rule", "gate_n")], idvar = "dataset",
             timevar = "rule", direction = "wide")
names(w) <- sub("gate_n\\.", "", names(w))
meta <- unique(res[, c("dataset", "n_cells", "k", "min_class")])
w <- merge(meta, w, by = "dataset"); w <- w[order(w$n_cells), ]
print(w, row.names = FALSE)

cat("\n======= classes whose specific markers CANNOT pass the gate =======\n")
cat("  (annotated classes smaller than the required number of expressing cells)\n")
w2 <- reshape(res[, c("dataset", "rule", "classes_excluded")], idvar = "dataset",
              timevar = "rule", direction = "wide")
names(w2) <- sub("classes_excluded\\.", "", names(w2))
w2 <- merge(meta, w2, by = "dataset"); w2 <- w2[order(w2$n_cells), ]
print(w2, row.names = FALSE)

cat("\n================ summary over datasets ================\n")
s <- do.call(rbind, lapply(split(res, res$rule), function(z) data.frame(
  rule = z$rule[1],
  datasets_with_empty_pool = sum(z$gate_n == 0),
  median_gate_n = stats::median(z$gate_n),
  median_pct_genes = round(stats::median(z$pct_genes), 2),
  total_classes_excluded = sum(z$classes_excluded),
  stringsAsFactors = FALSE)))
print(s[order(s$datasets_with_empty_pool, -s$total_classes_excluded * 0 +
              s$total_classes_excluded), ], row.names = FALSE)

cat("\nHOW TO READ THIS\n")
cat("  datasets_with_empty_pool > 0 means the method cannot run at all there.\n")
cat("  total_classes_excluded counts cell types that cannot contribute a marker\n")
cat("  under that rule -- a structural blind spot, independent of expression.\n")
cat("  A rule is only worth adopting if it lowers BOTH without inflating the\n")
cat("  pool so far that the ranking has nothing left to do (see the n/gate_n\n")
cat("  relationship: the dip ranking's gain scales with how much it discards).\n")
cat("\nwritten to: ", normalizePath(out_dir), "\n", sep = "")
