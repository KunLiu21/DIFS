# difs_runtime_memory.R: revision-2 analysis source.
# Source: code/as_run/06_revision2_reps_krange_abl2_runtime/difs_runtime_memory.R
# Usage, scope and limitations: benchmark/REVISION2.md.
difs_args <- function(args, known, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    if (!kv[2] %in% known) stop("unknown argument '", kv[2], "'. Accepts: ", paste(known, collapse = ", "))
    assign(kv[2], kv[3], envir = env)
  }
}
options(width = 200)
diag <- "../results/diag"; sacct <- "1"
difs_args(commandArgs(TRUE), c("diag", "sacct"))
METHODS <- c("DIFS", "Seurat", "FEAST", "Variance", "GateOnly")
GENES <- c(Kumar = 44508, Trapnell = 33128, SimKumar4easy = 42517, SimKumar4hard = 42534,
           Zhengmix4eq = 12885, Romanov = 18553, Zhengmix8eq = 13001, Lawlor = 19927,
           Koh = 41410, Darmanis = 19550, Muraro = 17718, Fletcher = 22730, Baron = 16359)

g <- utils::read.csv(file.path(diag, "grid_corrected.csv"), stringsAsFactors = FALSE)
g <- g[g$k_policy == "true" & g$n > 0 & g$fmethod %in% METHODS, ]
if (!all(unique(g$data.name) %in% names(GENES))) stop("dataset without a gene count: ",
  paste(setdiff(unique(g$data.name), names(GENES)), collapse = ", "))

agg <- function(x) c(median = stats::median(x), min = min(x), max = max(x), runs = length(x))
tt <- do.call(data.frame, aggregate(elapsed ~ data.name + fmethod + clust, data = g, FUN = agg))
names(tt) <- sub("^elapsed\\.", "time_", names(tt))
for (v in c("time_median", "time_min", "time_max")) tt[[v]] <- tt[[v]] / 60   # minutes
cells <- tapply(g$total, g$data.name, function(x) x[1])
tt$cells <- as.integer(cells[tt$data.name]); tt$genes <- GENES[tt$data.name]

tt$maxrss_GB <- NA_real_; tt$mem_tasks <- 0L
if (sacct != "0") {
  fj <- file.path(diag, "sacct_jobs.txt"); fs <- file.path(diag, "sacct_steps.txt")
  if (!file.exists(fj) || !file.exists(fs)) stop("sacct dumps not found in ", diag,
    " -- make them with the two commands in the header, or pass sacct=0")
  J <- utils::read.table(fj, sep = "|", header = TRUE, stringsAsFactors = FALSE, quote = "", fill = TRUE)
  S <- utils::read.table(fs, sep = "|", header = TRUE, stringsAsFactors = FALSE, quote = "", fill = TRUE)
  S <- S[grepl("\\.batch$", S$JobID), ]
  S$task <- sub("\\.batch$", "", S$JobID)
  mb <- function(x) { x <- trimws(x); u <- toupper(sub("^[0-9.]+", "", x)); v <- as.numeric(sub("[A-Za-z]+$", "", x))
    v * ifelse(u %in% c("", "M"), 1, ifelse(u == "G", 1024, ifelse(u == "K", 1 / 1024, ifelse(u == "T", 1024^2, NA)))) }
  J <- J[grepl("^[0-9]+_[0-9]+$", J$JobID) & J$State == "COMPLETED", ]
  J$array <- as.integer(sub("_.*", "", J$JobID)); J$index <- as.integer(sub(".*_", "", J$JobID))
  J$maxrss_MB <- mb(S$MaxRSS[match(J$JobID, S$task)])
  J$set <- sub("^difs_", "", J$JobName)
  J <- J[order(J$array, decreasing = TRUE), ]; J <- J[!duplicated(J[, c("set", "index")]), ]
  map <- do.call(rbind, lapply(c("ncurve", "baronfix"), function(s) {
    e <- new.env(); load(sprintf("../scenarios/scenarios_%s.Rdata", s), envir = e); sc <- e$scenarios
    data.frame(set = s, index = seq_len(nrow(sc)), data.name = gsub(".*/([^.]*).*", "\\1", sc$data.path),
               fmethod = sc$feature.selection.method, clust = sc$method, n = sc$n_features,
               k_policy = sc$k_policy, stringsAsFactors = FALSE)
  }))
  M <- merge(J[, c("set", "index", "maxrss_MB", "array")], map, by = c("set", "index"))
  M <- M[order(M$set != "baronfix"), ]
  M <- M[!duplicated(paste(M$data.name, M$fmethod, M$clust, M$n, M$k_policy)), ]
  M <- M[M$k_policy == "true" & M$n > 0 & M$fmethod %in% METHODS & is.finite(M$maxrss_MB), ]
  cat(sprintf("sacct: %d completed tasks with MaxRSS matched to analysis cells\n", nrow(M)))
  mm <- aggregate(maxrss_MB ~ data.name + fmethod + clust, data = M,
                  FUN = function(x) c(med = stats::median(x), n = length(x)))
  mm <- do.call(data.frame, mm)
  i <- match(paste(tt$data.name, tt$fmethod, tt$clust), paste(mm$data.name, mm$fmethod, mm$clust))
  tt$maxrss_GB <- mm$maxrss_MB.med[i] / 1024; tt$mem_tasks <- as.integer(mm$maxrss_MB.n[i])
  tt$mem_tasks[is.na(tt$mem_tasks)] <- 0L
}
tt <- tt[order(tt$clust, tt$cells, match(tt$fmethod, METHODS)), ]
utils::write.csv(tt, file.path(diag, "runtime_scaling.csv"), row.names = FALSE)

for (cl in unique(tt$clust)) {
  s <- tt[tt$clust == cl, ]
  w <- stats::reshape(s[, c("data.name", "cells", "genes", "fmethod", "time_median")],
                      idvar = c("data.name", "cells", "genes"), timevar = "fmethod", direction = "wide")
  names(w) <- sub("^time_median\\.", "", names(w)); w <- w[order(w$cells), ]
  w$DIFS_over_Seurat <- w$DIFS / w$Seurat
  cat(sprintf("\n== %s: median wall time per run, minutes (true k, 8 budgets) ==\n", cl))
  print(w, row.names = FALSE, digits = 3)
  cat("  log-log slope of time on cells (13 datasets; 1 = linear):\n")
  for (m in METHODS) { d <- s[s$fmethod == m, ]
    b <- stats::coef(stats::lm(log(time_median) ~ log(cells), data = d))[2]
    b2 <- stats::coef(stats::lm(log(time_median) ~ log(cells) + log(genes), data = d))
    cat(sprintf("    %-9s %.2f   (with genes in the model: cells %.2f, genes %.2f)\n", m, b, b2[2], b2[3])) }
  if (any(is.finite(s$maxrss_GB))) {
    w2 <- stats::reshape(s[, c("data.name", "cells", "fmethod", "maxrss_GB")],
                         idvar = c("data.name", "cells"), timevar = "fmethod", direction = "wide")
    names(w2) <- sub("^maxrss_GB\\.", "", names(w2)); w2 <- w2[order(w2$cells), ]
    cat(sprintf("\n== %s: median peak memory (MaxRSS) per run, GB ==\n", cl)); print(w2, row.names = FALSE, digits = 3)
  }
}
cat("\nNote: wall times are end to end (normalisation, feature selection, final clustering),\n",
    "on shared nodes with 4 cores; they are comparable across methods within a dataset,\n",
    "not to the second. DIFS's time does not depend on the budget, because the MSE size\n",
    "sweep and the stage II preliminary clustering run the same way at every budget\n",
    "(on Baron under Louvain, 55-77 min across the eight budgets).\n", sep = "")
