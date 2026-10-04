# difs_ranking_check.R: revision-2 analysis source.
# Source: code/as_run/07_feature_sets_gate_ranking/difs_ranking_check.R
# Usage, scope and limitations: benchmark/REVISION2.md.
difs_parse_args <- function(args, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) { message("ignoring argument: ", a); next }
    key <- kv[2L]; val <- sub('^"(.*)"$', "\\1", sub("^'(.*)'$", "\\1", kv[3L]))
    v <- if (grepl("^(TRUE|FALSE|T|F)$", val)) as.logical(val)
         else if (grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", val)) as.numeric(val)
         else val
    assign(key, v, envir = env)
  }
}
difs_parse_args(commandArgs(TRUE))

if (!exists("data",       inherits = FALSE)) data       <- "../source/Koh.rds"
if (!exists("out_dir",    inherits = FALSE)) out_dir    <- "../results/ranking_check"
if (!exists("B",          inherits = FALSE)) B          <- 2000
if (!exists("mc_reps",    inherits = FALSE)) mc_reps    <- 2
if (!exists("max_genes",  inherits = FALSE)) max_genes  <- 0      # 0 = 全部基因
if (!exists("min.expression",            inherits = FALSE)) min.expression <- log(5)
if (!exists("min.expression.proportion", inherits = FALSE)) min.expression.proportion <- 0.05
if (!exists("min.cells",  inherits = FALSE)) min.cells  <- 35
if (!exists("seed",       inherits = FALSE)) seed       <- 1

difs_fix_libpaths <- function(lib     = Sys.getenv("DIFS_LIB", ""),
                              damaged = Sys.getenv("DIFS_EXCLUDE_LIB", "")) {
  lp <- .libPaths()
  if (nzchar(damaged)) {
    d <- normalizePath(damaged, mustWork = FALSE)
    lp <- lp[normalizePath(lp, mustWork = FALSE) != d]
  }
  if (nzchar(lib) && dir.exists(lib)) lp <- unique(c(lib, lp))
  .libPaths(lp); Sys.setenv(R_LIBS = paste(.libPaths(), collapse = ":"))
}
difs_fix_libpaths()

suppressPackageStartupMessages({ library(Seurat); library(diptest) })
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("dataset        :", data, "\n")
cat("expression gate:", min.expression,
    if (abs(min.expression - log(5)) < 1e-8) " (log(5) -- as the code actually runs)"
    else if (abs(min.expression - 0.5) < 1e-8) " (0.5 -- as documented)" else "", "\n")
cat("B              :", B, "   Monte-Carlo repeats:", mc_reps, "\n\n")

seu <- readRDS(data)
mat <- as.matrix(GetAssayData(seu, slot = "data"))
n_cells <- ncol(mat)
cat("cells:", n_cells, "  genes:", nrow(mat), "\n")

floor_cells <- max(min.expression.proportion * n_cells, min.cells)
expr_list <- lapply(seq_len(nrow(mat)), function(i) {
  x <- mat[i, ]; x[x > min.expression]
})
n_expr <- vapply(expr_list, length, integer(1))
keep   <- which(n_expr >= floor_cells)
cat("genes passing the gate (>=", floor_cells, "expressing cells):", length(keep), "\n")

set.seed(seed)
if (max_genes > 0 && length(keep) > max_genes) {
  keep <- sort(sample(keep, max_genes))
  cat("NOTE: subsampled to", max_genes,
      "genes for speed -- top-N overlap is computed within this subsample\n")
}
g  <- rownames(mat)[keep]
xs <- expr_list[keep]
ns <- n_expr[keep]
cat("\n")

dipstat <- function(x) diptest::dip(sort(x))

cat("1/3  dip statistic and Hartigan p-value ...\n")
D_raw <- vapply(xs, dipstat, numeric(1))
p_hart <- vapply(xs, function(x) unname(diptest::dip.test(x)$p.value), numeric(1))

cat("2/3  Monte-Carlo normal calibration, B =", B, "x", mc_reps, "runs\n")
cat("     (this is the slow part; roughly", round(length(g) * B * mc_reps / 1.2e6),
    "minutes)\n")
mc_run <- function(run_seed) {
  set.seed(run_seed)
  out <- numeric(length(g)); t0 <- proc.time()[["elapsed"]]
  for (i in seq_along(g)) {
    null <- vapply(seq_len(B), function(b) dipstat(stats::rnorm(ns[i])), numeric(1))
    out[i] <- mean(D_raw[i] <= null)          # exactly as dip.test.v1 computes it
    if (i %% 1000 == 0)
      cat("       ", i, "/", length(g), "  (",
          round(proc.time()[["elapsed"]] - t0), "s )\n", sep = "")
  }
  out
}
p_mc  <- mc_run(seed * 100 + 1)
p_mc2 <- if (mc_reps >= 2) mc_run(seed * 100 + 2) else NULL

cat("3/3  comparing the rankings\n\n")

ord <- function(v, decreasing = FALSE) order(v, decreasing = decreasing)
topset <- function(v, N, decreasing = FALSE) g[ord(v, decreasing)][seq_len(min(N, length(g)))]
jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))

keys <- list(mc_normal = p_mc, hartigan = p_hart, D_raw = -D_raw)  # 全部改成"小的在前"
if (!is.null(p_mc2)) keys$mc_normal_run2 <- p_mc2

cat("========== ties ==========\n")
tie <- data.frame(
  key = c("p_mc_normal", "p_hartigan", "D_raw"),
  distinct_values = c(length(unique(p_mc)), length(unique(p_hart)), length(unique(D_raw))),
  at_boundary = c(sum(p_mc == 0), sum(p_hart >= 1), 0L),
  n_genes = length(g))
tie$distinct_frac <- round(tie$distinct_values / tie$n_genes, 4)
print(tie, row.names = FALSE)
cat("\n  p_mc_normal 在 0 处的并列数 = 无法区分的'最显著'基因数;\n")
cat("  sort() 按行序打破这些并列, 所以 top-N 的这部分由矩阵行序决定。\n\n")

cat("========== rank correlation (Spearman) ==========\n")
nm <- names(keys)
cm <- outer(seq_along(keys), seq_along(keys),
            Vectorize(function(i, j) round(stats::cor(keys[[i]], keys[[j]],
                                                      method = "spearman"), 4)))
dimnames(cm) <- list(nm, nm); print(cm)

cat("\n========== top-N overlap ==========\n")
P <- length(g)
Ns <- c(100, 200, 500, 1000, 2000)
Ns <- Ns[Ns < P]
chance <- function(N) N / (2 * P - N)
adj <- function(obs, N) { c <- chance(N); (obs - c) / (1 - c) }

pairs <- list(c("mc_normal","hartigan"), c("mc_normal","D_raw"), c("hartigan","D_raw"))
if (!is.null(p_mc2)) pairs <- c(pairs, list(c("mc_normal","mc_normal_run2")))

cat("pool of gate-passing genes P =", P, "\n\n")
ov <- do.call(rbind, lapply(Ns, function(N) {
  r <- data.frame(topN = N, chance = round(chance(N), 3))
  for (p in pairs) {
    o <- jac(topset(keys[[p[1]]], N), topset(keys[[p[2]]], N))
    r[[paste0(paste(p, collapse = "~"), " raw")]] <- round(o, 3)
    r[[paste0(paste(p, collapse = "~"), " adj")]] <- round(adj(o, N), 3)
  }
  r
}))
print(ov, row.names = FALSE)
cat("\n  raw = 原始 Jaccard (会被小的基因池抬高, 不要单独引用)\n")
cat("  adj = 扣除随机基线后的一致性, 这才是可比较的数字\n")
cat("        adj ~ 1  两个排序基本相同\n")
cat("        adj ~ 0  两个排序与随机无异\n")

cat("\n========== 表达细胞数 n 在 top-1000 中的分布 ==========\n")
cat("(参照分布只能通过 n 影响排序, 所以这里看它是否偏向 n 小的基因)\n")
N <- min(1000, length(g))
nn <- data.frame(
  key = c("all tested", "p_mc_normal", "p_hartigan", "D_raw"),
  median_n = c(stats::median(ns),
               stats::median(ns[ord(p_mc)][seq_len(N)]),
               stats::median(ns[ord(p_hart)][seq_len(N)]),
               stats::median(ns[ord(-D_raw)][seq_len(N)])))
print(nn, row.names = FALSE)

res <- data.frame(gene = g, n_expressing = ns, D = D_raw,
                  p_mc_normal = p_mc, p_hartigan = p_hart,
                  stringsAsFactors = FALSE)
if (!is.null(p_mc2)) res$p_mc_normal_run2 <- p_mc2
f <- file.path(out_dir, paste0("ranking_check_",
     sub("\\.rds$", "", basename(data)), ".csv"))
utils::write.csv(res, f, row.names = FALSE)
utils::write.csv(ov, file.path(out_dir, "overlap_summary.csv"), row.names = FALSE)

cat("\n每个基因的三个排序键已写入:\n  ", normalizePath(f), "\n")
cat("\n怎么读这个结果 (一律看 adj 列, 不要看 raw)\n")
cat("  mc_normal~hartigan adj\n")
cat("      > 0.9  参照分布几乎不改变排序 -> 在论文中降级为实现细节\n")
cat("      0.5-0.9 排序有实质差异, 但要靠 E4 才能判断哪个更好\n")
cat("      < 0.3  两种参照给出几乎不同的排序, 这个选择必须论证\n")
cat("  mc_normal~D_raw adj 低 -> 按 p 排序 (校正表达细胞数) 与按原始统计量\n")
cat("      排序是两回事; 结合上面 n 的中位数说明为什么前者更合理\n")
cat("  mc_normal~mc_normal_run2 adj 接近 1 -> 现有实现自身可复现\n")
cat("\n注意: 通过门限的基因池只有", P, "个, 而 num_features_decider 的扫描\n")
cat("      范围是 500-2500. 当它选到 >=", P, "时, 排序对最终特征集不再起作用,\n")
cat("      超出的部分是未经检验的基因按矩阵行序填充的.\n")
