
difs_args <- function(args, known, env = parent.frame()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    if (!kv[2] %in% known)
      stop("unknown argument '", kv[2], "'. This script accepts: ",
           paste(known, collapse = ", "))
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
  invisible(NULL)
}
args <- commandArgs(TRUE); out <- "../results/diag"; difs_args(args, c("out"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source("difs_validate.R")
source("difs_signrank.R")   ## exact paired test; see that file for why not wilcox.test

FIXED    <- c(100, 200, 300, 500, 750, 1000, 1500, 2500)
NOT_IND  <- c("SimKumar4easy", "SimKumar4hard", "Zhengmix4eq", "Zhengmix8eq")
THREE    <- c("DIFS", "FEAST", "Seurat")

read_set <- function(set) {
  d  <- file.path("..", "results", set)
  fs <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  if (!length(fs)) stop("no results under ", d)
  g <- function(x, n, def = NA) { v <- x[[n]]; if (is.null(v) || !length(v)) def else v[1] }
  do.call(rbind, lapply(fs, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x)) { cat("*** unreadable: ", f, "\n", sep = ""); return(NULL) }
    data.frame(
      data.name = as.character(g(x, "data.name", "")),
      fmethod   = as.character(g(x, "feature.selection.method", "")),
      clust     = as.character(g(x, "clustering.method", "")),
      n         = as.integer(g(x, "n_requested", NA)),
      k_policy  = as.character(g(x, "k_policy", "")),
      gate_rule = as.character(g(x, "gate_rule", "submitted")),
      gate_rule_src = if (is.null(x[["gate_rule"]]) || !length(x[["gate_rule"]]))
                        "default" else "object",
      seed      = as.integer(g(x, "seed", NA)),
      ARI       = as.numeric(g(x, "ARI", NA)),
      FM        = as.numeric(g(x, "FM", NA)),
      Jaccard   = as.numeric(g(x, "Jaccard", NA)),
      Purity    = as.numeric(g(x, "Purity", NA)),
      nfeat     = as.integer(g(x, "n_features", NA)),
      kfound    = as.integer(g(x, "n_clusters_found", NA)),
      k_used    = as.integer(g(x, "k_used_clustering", NA)),
      stage1_n  = as.integer(g(x, "difs_stage1_n", NA)),
      stage2_n  = as.integer(g(x, "difs_stage2_n", NA)),
      ratio     = as.character(g(x, "difs_ratio", NA)),
      total     = as.integer(g(x, "n_cells_total", NA)),
      dropped   = as.integer(g(x, "n_cells_dropped", NA)),
      elapsed   = as.numeric(g(x, "elapsed_sec", NA)),
      set       = set, stringsAsFactors = FALSE)
  }))
}

fix  <- read_set("baronfix")
grid <- read_set("ncurve")
abl  <- read_set("abl")
expected <- expected_from_scenarios("baronfix")
cat("baronfix rows: ", nrow(fix), "   expected: ", nrow(expected),
    "   main grid: ", nrow(grid), "   ablation: ", nrow(abl), "\n\n", sep = "")

difs_gate(fix, expected, old = NULL, require_full_labels = TRUE,
          required_new = c("ARI", "FM", "Jaccard", "Purity", "nfeat", "total"),
          expect_total = 8569L)

KEY  <- c("data.name", "fmethod", "clust", "n", "k_policy", "gate_rule", "seed")
main <- fix[fix$fmethod != "DIFS_hartigan", ]
kmain <- difs_key_string(main); kgrid <- difs_key_string(grid)
orphan <- setdiff(kmain, kgrid)
if (length(orphan)) {
  cat("*** ", length(orphan), " repaired run(s) have no counterpart in the main grid.\n",
      "    A substitution that adds rows is not a substitution. Stopping.\n", sep = "")
  print(utils::head(orphan, 10)); quit(save = "no", status = 2)
}
before <- grid
after  <- grid
idx <- match(kmain, kgrid)
if (anyNA(idx)) stop("internal: orphan check passed but match() returned NA")
if (anyDuplicated(idx)) stop("a repaired run matches more than one grid row")

CARRY <- c("ARI", "FM", "Jaccard", "Purity", "nfeat", "kfound", "k_used",
           "stage1_n", "stage2_n", "ratio", "total", "dropped", "elapsed")
stopifnot(all(CARRY %in% names(after)), all(CARRY %in% names(main)))
for (cc in CARRY) after[[cc]][idx] <- main[[cc]]
after$set[idx] <- "baronfix"

chk <- after[idx, ]
if (any(!is.finite(chk$ARI) | !is.finite(chk$FM) |
        !is.finite(chk$Jaccard) | !is.finite(chk$Purity)))
  stop("a substituted row has a non-finite metric")
if (any(chk$dropped != 0L) || any(chk$total != 8569L))
  stop("a substituted row does not cover all 8569 cells after substitution -- ",
       "the carried fields and the gated fields disagree")
moved <- abs(after$ARI[idx] - before$ARI[idx]) > 1e-12
cat("substituted ", length(idx), " of ", nrow(grid), " main-grid rows (all Baron); ",
    sum(moved), " changed ARI, range ",
    sprintf("%+.4f .. %+.4f", min(after$ARI[idx] - before$ARI[idx]),
                              max(after$ARI[idx] - before$ARI[idx])), "\n",
    "every carried field replaced: ", paste(CARRY, collapse = ", "), "\n\n", sep = "")

cell_of <- function(d, methods = THREE) {
  s <- d[d$n %in% FIXED & d$fmethod %in% methods, ]
  aggregate(ARI ~ k_policy + data.name + n + fmethod, data = s, FUN = mean)
}
mean_rank <- function(d, kp, drop_baron = FALSE, methods = THREE) {
  s <- cell_of(d, methods); s <- s[s$k_policy == kp, ]
  if (drop_baron) s <- s[s$data.name != "Baron", ]
  w <- reshape(s[, c("data.name","n","fmethod","ARI")], idvar = c("data.name","n"),
               timevar = "fmethod", direction = "wide")
  w <- w[stats::complete.cases(w), ]
  R <- t(apply(-as.matrix(w[, -(1:2)]), 1, rank, ties.method = "average"))
  colnames(R) <- sub("^ARI\\.", "", colnames(w)[-(1:2)])
  list(n = nrow(w), mean = colMeans(R), first = colSums(R == 1))
}
show_rank <- function(lab, a, b) {
  cat(sprintf("  %-34s %s\n", lab,
      paste(sprintf("%s %.3f->%.3f", names(a$mean), a$mean, b$mean[names(a$mean)]),
            collapse = "   ")))
  cat(sprintf("  %-34s first: %s   (of %d)\n", "",
      paste(sprintf("%s %d->%d", names(a$first), a$first, b$first[names(a$first)]),
            collapse = "  "), a$n))
}
cat(strrep("=", 76), "\nMEAN RANK  (before -> after)\n", strrep("=", 76), "\n", sep = "")
for (kp in c("true", "estimated"))
  show_rank(paste0("k = ", kp), mean_rank(before, kp), mean_rank(after, kp))
show_rank("k = true, Baron excluded",
          mean_rank(before, "true", TRUE), mean_rank(after, "true", TRUE))

paired <- function(d, kp, bud, ds = NULL) {
  s <- cell_of(d); s <- s[s$k_policy == kp & s$n == bud, ]
  if (!is.null(ds)) s <- s[s$data.name %in% ds, ]
  w <- reshape(s[, c("data.name","fmethod","ARI")], idvar = "data.name",
               timevar = "fmethod", direction = "wide")
  w <- w[stats::complete.cases(w), ]
  sapply(c("Seurat","FEAST"), function(m) {
    x <- w[["ARI.DIFS"]]; y <- w[[paste0("ARI.", m)]]
    p <- difs_p(x, y)
    c(mean = mean(x - y), better = sum(x > y), n = length(x), p = p)
  })
}
cat("\n", strrep("=", 76), "\nPAIRED AT n = 100, k = true, 13 datasets (before -> after)\n",
    strrep("=", 76), "\n", sep = "")
a <- paired(before, "true", 100); b <- paired(after, "true", 100)
for (m in colnames(a))
  cat(sprintf("  DIFS - %-7s  mean %+.4f -> %+.4f   better %d/%d -> %d/%d   p %.4f -> %.4f\n",
      m, a["mean",m], b["mean",m], a["better",m], a["n",m], b["better",m], b["n",m],
      a["p",m], b["p",m]))

IND <- setdiff(unique(grid$data.name), NOT_IND)
cat("\n", strrep("=", 76), "\nTHE NINE INDEPENDENT REAL DATASETS, k = true (before -> after)\n",
    strrep("=", 76), "\n", sep = "")
for (bud in c(100, 500, 2500)) {
  a <- paired(before, "true", bud, IND); b <- paired(after, "true", bud, IND)
  cat(sprintf("  n=%-5d vs Seurat  %+.4f -> %+.4f (p %.3f -> %.3f)   vs FEAST  %+.4f -> %+.4f (p %.3f -> %.3f)\n",
      bud, a["mean","Seurat"], b["mean","Seurat"], a["p","Seurat"], b["p","Seurat"],
      a["mean","FEAST"],  b["mean","FEAST"],  a["p","FEAST"],  b["p","FEAST"]))
}
grid_mean <- function(d, m, bud, ds = NULL) {
  s <- cell_of(d); s <- s[s$k_policy == "true" & s$n == bud & s$fmethod == m, ]
  if (!is.null(ds)) s <- s[s$data.name %in% ds, ]; mean(s$ARI)
}
cat(sprintf("\n  DIFS@300 vs Seurat@2500, nine independent: %.4f -> %.4f  vs  %.4f -> %.4f\n",
    grid_mean(before,"DIFS",300,IND), grid_mean(after,"DIFS",300,IND),
    grid_mean(before,"Seurat",2500,IND), grid_mean(after,"Seurat",2500,IND)))

cat("\n", strrep("=", 76), "\nCOMPONENT DECOMPOSITION at n = 100 (before -> after)\n",
    strrep("=", 76), "\n", sep = "")
abl_fix <- abl
rep_rows <- fix[fix$data.name == "Baron" & fix$fmethod == "DIFS" &
                fix$clust == "Refined Louvain" & fix$k_policy == "true" &
                fix$n %in% c(100, 300, 1000), ]
if (nrow(rep_rows) != 3L)
  stop("expected 3 repaired Baron DIFS/Louvain/true cells for the ablation, found ",
       nrow(rep_rows))
for (i in seq_len(nrow(rep_rows))) {
  j <- which(abl_fix$data.name == "Baron" & abl_fix$fmethod == "DIFS" &
             abl_fix$n == rep_rows$n[i])
  if (length(j) != 1L)
    stop("ablation Baron DIFS at n = ", rep_rows$n[i], " matched ", length(j), " rows")
  for (cc in CARRY) abl_fix[[cc]][j] <- rep_rows[[cc]][i]
  abl_fix$set[j] <- "baronfix"
}
cat("ablation: 3 Baron DIFS cells refreshed from the repaired runs (whole record)\n")
decomp <- function(d, bud, ds = NULL) {
  s <- d[d$n == bud, ]; if (!is.null(ds)) s <- s[s$data.name %in% ds, ]
  w <- reshape(s[, c("data.name","fmethod","ARI")], idvar = "data.name",
               timevar = "fmethod", direction = "wide")
  w <- w[stats::complete.cases(w), ]
  s1 <- grep("stage1only", colnames(w), value = TRUE)[1]
  r <- w[[s1]] - w[["ARI.GateOnly"]]; t2 <- w[["ARI.DIFS"]] - w[[s1]]
  pp <- function(x) difs_p(x)
  c(rank = mean(r), rank_p = pp(r), st2 = mean(t2), st2_p = pp(t2),
    total = mean(r + t2), total_p = pp(r + t2), n = nrow(w))
}
ABIND <- setdiff(unique(abl$data.name), NOT_IND)
for (lab in c("all 13", "9 independent")) {
  ds <- if (lab == "all 13") NULL else ABIND
  for (bud in c(100, 300, 1000)) {
    A <- decomp(abl, bud, ds); B <- decomp(abl_fix, bud, ds)
    cat(sprintf("  %-14s n=%-5d rank %+.4f->%+.4f  stageII %+.4f->%+.4f  total %+.4f->%+.4f  (%d ds)\n",
        lab, bud, A["rank"], B["rank"], A["st2"], B["st2"], A["total"], B["total"], A["n"]))
  }
}

cat("\n", strrep("=", 76), "\nBARON, BOTH ARMS REPAIRED (for the hart comparison)\n",
    strrep("=", 76), "\n", sep = "")
h  <- fix[fix$fmethod == "DIFS_hartigan", ]
dd <- fix[fix$fmethod == "DIFS" & fix$k_policy == "true" & fix$n %in% c(100, 300), ]
if (nrow(h) != 4L || nrow(dd) != 4L)
  stop("the reference-distribution comparison needs 4 rows per arm; found ",
       nrow(h), " and ", nrow(dd))
mm <- merge(h, dd, by = c("data.name","clust","n","k_policy","gate_rule","seed"),
            suffixes = c("_h","_d"))
mm$dARI <- mm$ARI_d - mm$ARI_h
if (nrow(mm) != 4L) stop("Baron hart pairing produced ", nrow(mm), " rows, expected 4")
if (any(mm$dropped_d != 0L) || any(mm$dropped_h != 0L))
  stop("a Baron hart arm still records dropped cells")
print(mm[order(mm$n, mm$clust), c("clust","n","ARI_d","ARI_h","dARI")],
      row.names = FALSE, digits = 4)
cat("\nFeed these Baron values back into difs_hart_compare.R before restating the\n",
    "reference-distribution result: its significance at 100 features depended on\n",
    "Baron, and until now both of its Baron arms carried the defect.\n", sep = "")

utils::write.csv(after,  file.path(out, "grid_corrected.csv"),  row.names = FALSE)
utils::write.csv(before, file.path(out, "grid_asrun.csv"),      row.names = FALSE)
utils::write.csv(abl_fix, file.path(out, "abl_corrected.csv"),  row.names = FALSE)
utils::write.csv(mm, file.path(out, "baron_hart_repaired.csv"), row.names = FALSE)
cat("\nwritten to ", normalizePath(out), ":\n  grid_corrected.csv  grid_asrun.csv",
    "  abl_corrected.csv  baron_hart_repaired.csv\n", sep = "")
