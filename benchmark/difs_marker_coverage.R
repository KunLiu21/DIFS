# Marker coverage; numerical ties use 12 decimal places. See benchmark/REVISION2.md.
difs_args <- function(args, known, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    if (!kv[2] %in% known) stop("unknown argument '", kv[2], "'. Accepts: ",
                                paste(known, collapse = ", "))
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
}
task <- "run"; k <- NA; out_dir <- "../results/markers"; seed <- 1
difs_args(commandArgs(TRUE), c("task", "k", "out_dir", "seed"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BUDGETS <- c(100, 200, 300, 500, 750, 1000, 1500, 2500)
DATASETS <- c("Baron", "Darmanis", "Fletcher", "Koh", "Kumar", "Lawlor", "Muraro", "Romanov",
              "SimKumar4easy", "SimKumar4hard", "Trapnell", "Zhengmix4eq", "Zhengmix8eq")
INDEPENDENT <- c("Kumar", "Trapnell", "Romanov", "Lawlor", "Koh", "Darmanis",
                 "Muraro", "Fletcher", "Baron")
MIN_TYPE_CELLS <- 10; Q_MAX <- 0.05; LFC_MIN <- 1; PCT_MIN <- 0.25

source_files <- function() {
  fs <- list.files(file.path("..", "source"), pattern = "\\.rds$", full.names = TRUE)
  setNames(fs, gsub(".*/([^.]*).*", "\\1", fs))
}
task_table <- function() {
  fs <- source_files()
  if (anyDuplicated(names(fs))) stop("two source files share a dataset name")
  if (length(setdiff(DATASETS, names(fs))))
    stop("not in ../source: ", paste(setdiff(DATASETS, names(fs)), collapse = ", "))
  ds <- DATASETS
  rbind(data.frame(task = "prep", dataset = ds, budget = NA_integer_),
        data.frame(task = "difs", dataset = rep(ds, each = length(BUDGETS)),
                   budget = rep(BUDGETS, length(ds))))
}

if (task == "list") {
  tt <- task_table(); tt$index <- seq_len(nrow(tt))
  print(tt, row.names = FALSE)
  b <- tt$index[tt$dataset == "Baron"]
  cat("\nBaron tasks (64G, 12 h): ", paste(b, collapse = ","),
      "\nall others (16G, 4 h): 2-13,", paste(range(setdiff(tt$index[tt$task == "difs"], b)),
      collapse = "-"), "\n", sep = "")
  quit(save = "no")
}

de_markers <- function(X, lab) {
  ok <- !is.na(lab); X <- X[, ok, drop = FALSE]; lab <- as.character(lab[ok])
  N <- ncol(X); tab <- table(lab)
  types <- names(tab)[tab >= MIN_TYPE_CELLS]
  if (length(types) < 2) stop("fewer than two labelled types with >= ", MIN_TYPE_CELLS, " cells")
  R <- if (requireNamespace("matrixStats", quietly = TRUE))
         matrixStats::rowRanks(X, ties.method = "average") else t(apply(X, 1, rank))
  tie <- apply(X, 1, function(x) { t <- tabulate(match(x, unique(x))); sum(t^3 - t) })
  E <- expm1(X); P <- X > 0
  do.call(rbind, lapply(types, function(t) {
    i <- lab == t; n1 <- sum(i); n0 <- N - n1
    U <- rowSums(R[, i, drop = FALSE]) - n1 * (n1 + 1) / 2
    s2 <- n1 * n0 / 12 * ((N + 1) - tie / (N * (N - 1)))
    z <- (U - n1 * n0 / 2 - 0.5) / sqrt(s2)              # continuity-corrected, one-sided
    p <- stats::pnorm(z, lower.tail = FALSE); p[!is.finite(p)] <- 1
    data.frame(gene = rownames(X), type = t, n_type = n1, auc = U / (n1 * n0),
               p = p, q = stats::p.adjust(p, "BH"),
               log2fc = log2(rowMeans(E[, i, drop = FALSE]) + 1) -
                        log2(rowMeans(E[, !i, drop = FALSE]) + 1),
               pct_in = rowMeans(P[, i, drop = FALSE]), pct_out = rowMeans(P[, !i, drop = FALSE]),
               stringsAsFactors = FALSE)
  }))
}
reference_lists <- function(de) {
  de <- de[de$q < Q_MAX & de$log2fc >= LFC_MIN & de$pct_in >= PCT_MIN, ]
  de <- de[order(de$type, -de$auc, -de$log2fc), ]
  de$rank_in_type <- stats::ave(seq_len(nrow(de)), de$type, FUN = seq_along)
  de
}

hurdle_markers <- function(X, lab) {
  ok <- !is.na(lab); X <- X[, ok, drop = FALSE]; lab <- as.character(lab[ok]); N <- ncol(X)
  tab <- table(lab); types <- names(tab)[tab >= MIN_TYPE_CELLS]
  if (length(types) < 2) stop("fewer than two labelled types with >= ", MIN_TYPE_CELLS, " cells")
  E <- X > 0; X2 <- X^2
  kT <- rowSums(E); sT <- rowSums(X); ssT <- rowSums(X2)
  clip <- function(p) pmin(pmax(p, 1e-5), 1 - 1e-5)
  llD <- function(k, n) { a <- clip(k / n); (n - k) * log(1 - a) + k * log(a) }
  llC <- function(k, s, ss) {
    m <- ifelse(k > 0, s / pmax(k, 1), 0); rss <- pmax(ss - k * m^2, 0)
    v <- ifelse(k >= 2, rss / pmax(k - 1, 1), 1); v <- pmax(v, 1e-12)
    ifelse(k > 0, -k / 2 * log(2 * pi * v) - rss / (2 * v), 0)
  }
  totE <- rowSums(expm1(X))
  out <- lapply(types, function(t) {
    i <- lab == t; n1 <- sum(i); n0 <- N - n1
    k1 <- rowSums(E[, i, drop = FALSE]); s1 <- rowSums(X[, i, drop = FALSE]); ss1 <- rowSums(X2[, i, drop = FALSE])
    k0 <- kT - k1; s0 <- sT - s1; ss0 <- ssT - ss1
    LD <- pmax(2 * (llD(k1, n1) + llD(k0, n0) - llD(kT, N)), 0)
    LC <- pmax(2 * (llC(k1, s1, ss1) + llC(k0, s0, ss0) - llC(kT, sT, ssT)), 0)
    pD <- stats::pchisq(LD, 1, lower.tail = FALSE); pC <- stats::pchisq(LC, 2, lower.tail = FALSE)
    pH <- stats::pchisq(LD + LC, 3, lower.tail = FALSE)
    e1 <- rowSums(expm1(X[, i, drop = FALSE]))
    data.frame(gene = rownames(X), type = t, n_type = n1,
               LR_detect = LD, LR_level = LC, LR_hurdle = LD + LC,
               q_detect = stats::p.adjust(pD, "BH"), q_level = stats::p.adjust(pC, "BH"),
               q_hurdle = stats::p.adjust(pH, "BH"),
               pct_in = k1 / n1, pct_out = k0 / n0, n_expr_in = k1,
               mean_expr_in = ifelse(k1 > 0, s1 / pmax(k1, 1), NA_real_),
               mean_expr_out = ifelse(k0 > 0, s0 / pmax(k0, 1), NA_real_),
               log2fc = log2(e1 / n1 + 1) - log2((totE - e1) / n0 + 1), stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}
rank_within <- function(d, key) {
  d <- d[order(d$type, -d[[key]]), , drop = FALSE]
  d$rank_in_type <- stats::ave(seq_len(nrow(d)), d$type, FUN = seq_along); d
}
hurdle_lists <- function(h) list(
  hurdle = rank_within(h[which(h$q_hurdle < Q_MAX & h$log2fc >= LFC_MIN & h$pct_in >= PCT_MIN), ], "LR_hurdle"),
  level = rank_within(h[which(h$q_level < Q_MAX & h$n_expr_in >= MIN_TYPE_CELLS & h$pct_in >= PCT_MIN &
                              (h$mean_expr_in - h$mean_expr_out) >= log(2)), ], "LR_level"),
  detection = rank_within(h[which(h$q_detect < Q_MAX & (h$pct_in - h$pct_out) >= 0.25), ], "LR_detect"))

load_dataset <- function(d) {
  suppressPackageStartupMessages(library(Seurat))
  seu <- readRDS(source_files()[[d]])
  list(seu = seu, lm = as.matrix(GetAssayData(seu, slot = "data")),
       cnt = as.matrix(GetAssayData(seu, slot = "counts")), lab = seu$trueclass)
}
k_true_of <- function(lab) length(unique(lab[!is.na(lab)]))

run_prep <- function(d) {
  t0 <- proc.time()[["elapsed"]]
  D <- load_dataset(d); kt <- k_true_of(D$lab)
  set.seed(seed)
  hu <- hurdle_markers(D$lm, D$lab); refs <- hurdle_lists(hu)
  hu <- hu[which(pmin(hu$q_detect, hu$q_level, hu$q_hurdle) < Q_MAX), ]; gc()
  de <- de_markers(D$lm, D$lab); refs$wilcox <- reference_lists(de); gc()
  ref <- refs$hurdle
  gate <- difs_gate_genes(D$lm, log(5), gate_rule = "submitted", gate_k = kt)
  rk <- list(
    DIFS_stage1only = difs_stage1_ranking(D$lm, min.expression = log(5),
                                          gate_rule = "submitted", gate_k = kt,
                                          key = "normal_lookup", seed = seed),
    Seurat = VariableFeatures(FindVariableFeatures(D$seu, selection.method = "vst",
                                                   nfeatures = nrow(D$seu), verbose = FALSE)),
    Seurat_gated = feature_ranking("Seurat_gated", data.seu = D$seu, data.log.mat = D$lm,
                                   data.mat = D$cnt, cds = NULL, cluster_count = kt,
                                   max_n = nrow(D$lm), seed = seed, min.expression = log(5),
                                   gate_k = kt),
    FEAST = feature_ranking("FEAST", data.seu = D$seu, data.log.mat = D$lm, data.mat = D$cnt,
                            cds = NULL, cluster_count = kt, max_n = nrow(D$lm), seed = seed))
  out <- list(dataset = d, k_true = kt, n_cells = ncol(D$lm), genes = rownames(D$lm),
              gate = gate, rankings = rk, de_all = de[de$q < Q_MAX, ], hurdle_all = hu,
              reference = ref, references = refs,
              elapsed_sec = proc.time()[["elapsed"]] - t0)
  saveRDS(out, file.path(out_dir, sprintf("prep__%s.rds", d)))
  cat(sprintf("%s: k = %d, qualifying rows %s, %d types, gate %d, rankings %s (%.0f s)\n",
              d, kt, paste(sprintf("%s=%d", names(refs), vapply(refs, nrow, 0L)), collapse = " "),
              length(unique(ref$type)), length(gate),
              paste(sprintf("%s=%d", names(rk), lengths(rk)), collapse = " "), out$elapsed_sec))
}

run_difs <- function(d, b) {
  t0 <- proc.time()[["elapsed"]]
  D <- load_dataset(d); kt <- k_true_of(D$lab)
  set.seed(seed)
  r <- difs_features(data.log.mat = D$lm, data.seu = D$seu, cluster_count = kt,
                     binary.bound = binary.bound.decider(D$lm), min.expression = log(5),
                     gate_rule = "submitted", gate_k = kt,
                     tune_stage1 = TRUE, tune_ratio = TRUE, stage1_key = "normal_lookup",
                     use_stage2 = TRUE, seed = seed, tops = DIFS_TOPS_GRID, n_fixed = b)
  out <- list(dataset = d, budget = b, k_true = kt,
              genes = intersect(r$markers, rownames(D$seu)),
              stage1_n = r$stage1_n, stage2_n = r$stage2_n, ratio = r$ratio,
              elapsed_sec = proc.time()[["elapsed"]] - t0)
  saveRDS(out, file.path(out_dir, sprintf("difs__%s__%d.rds", d, b)))
  cat(sprintf("%s b=%d: %d genes, stage I %s, stage II %s, ratio %s (%.0f s)\n", d, b,
              length(out$genes), r$stage1_n, r$stage2_n, r$ratio, out$elapsed_sec))
}

if (task %in% c("run", "smoke")) {
  source("Functions_controlled.R")
  if (!file.exists("difs_sc3_complete.R")) stop("difs_sc3_complete.R is required (SC3 repair)")
  source("difs_sc3_complete.R")
  if (task == "smoke") {
    out_dir <- file.path(out_dir, "smoke"); dir.create(out_dir, showWarnings = FALSE)
    run_prep("Koh"); for (b in c(100, 300)) run_difs("Koh", b)
    task <- "summary"
  } else {
    tt <- task_table()
    if (is.na(k) || k < 1 || k > nrow(tt)) stop("k= must be 1..", nrow(tt), " (task=list)")
    row <- tt[k, ]; cat(sprintf("task %d: %s %s %s\n", k, row$task, row$dataset,
                                ifelse(is.na(row$budget), "", row$budget)))
    if (row$task == "prep") run_prep(row$dataset) else run_difs(row$dataset, row$budget)
    quit(save = "no", status = 0)
  }
}

if (task != "summary") stop("task= must be list, run (with k=), smoke or summary")

hyp_mean <- function(m, P, b) if (P <= 0) 0 else m * min(b, P) / P
hyp_atleast <- function(j, m, P, b) {
  b <- min(b, P); if (m < j) return(0)
  stats::phyper(j - 1, m, P - m, b, lower.tail = FALSE)
}
signrank_exact <- function(d) {        # sign-permutation, zeros dropped, midranks
  d <- round(d[is.finite(d)], 12)        # 2026-10-04: equal differences must tie
  d <- d[d != 0]; m <- length(d)
  if (m < 2) return(NA_real_)
  r <- rank(abs(d)); obs <- sum(r[d > 0])
  S <- as.matrix(expand.grid(rep(list(c(0, 1)), m)))
  null <- S %*% r
  mean(abs(null - sum(r) / 2) >= abs(obs - sum(r) / 2) - 1e-9)
}

LISTS <- c("hurdle", "level", "detection", "wilcox")
pf <- list.files(out_dir, pattern = "^prep__.*\\.rds$", full.names = TRUE)
if (!length(pf)) stop("no prep__*.rds in ", out_dir)
rows <- list(); refs_out <- list()
for (f in pf) {
  P <- readRDS(f); d <- P$dataset
  RL <- if (!is.null(P$references)) P$references else list(wilcox = P$reference)
  G <- length(P$genes); gate <- P$gate
  sets_at <- function(b) {
    s <- lapply(P$rankings, function(r) utils::head(r, b))
    df <- file.path(out_dir, sprintf("difs__%s__%d.rds", d, b))
    if (file.exists(df)) s$DIFS <- readRDS(df)$genes
    s
  }
  for (L in intersect(LISTS, names(RL))) {
    ref <- RL[[L]]; if (!nrow(ref)) next
    top10 <- ref[ref$rank_in_type <= 10, ]; top25 <- ref[ref$rank_in_type <= 25, ]
    M10 <- unique(top10$gene); M25 <- unique(top25$gene); Mall <- unique(ref$gene)
    types <- unique(top10$type)
    refs_out[[paste(d, L)]] <- cbind(dataset = d, list = L, top25[, c("gene", "type", "rank_in_type")])
    for (b in BUDGETS) {
      S <- sets_at(b)
      for (m in names(S)) {
        g <- S[[m]]
        per_type <- vapply(types, function(t) sum(top10$gene[top10$type == t] %in% g), numeric(1))
        rows[[length(rows) + 1]] <- data.frame(dataset = d, list = L, method = m, budget = b,
          n_selected = length(g), cov10 = mean(M10 %in% g), cov25 = mean(M25 %in% g),
          covall = mean(Mall %in% g), types1 = mean(per_type >= 1), types3 = mean(per_type >= 3),
          precision = if (length(g)) mean(g %in% Mall) else NA_real_)
      }
      for (pool in c("GateOnly", "Random")) {
        U <- if (pool == "GateOnly") gate else P$genes; n <- length(U)
        e <- function(M) hyp_mean(sum(M %in% U), n, b) / length(M)
        pt <- vapply(types, function(t) {
          mt <- sum(top10$gene[top10$type == t] %in% U)
          c(hyp_atleast(1, mt, n, b), hyp_atleast(3, mt, n, b)) }, numeric(2))
        rows[[length(rows) + 1]] <- data.frame(dataset = d, list = L, method = pool, budget = b,
          n_selected = min(b, n), cov10 = e(M10), cov25 = e(M25), covall = e(Mall),
          types1 = mean(pt[1, ]), types3 = mean(pt[2, ]), precision = sum(Mall %in% U) / n)
      }
    }
    need <- function(r, share) { h <- cumsum(r %in% M10); w <- which(h >= ceiling(share * length(M10)))
                                 if (length(w)) w[1] else NA_integer_ }
    for (m in names(P$rankings)) rows[[length(rows) + 1]] <- data.frame(dataset = d, list = L, method = m,
      budget = NA, n_selected = length(P$rankings[[m]]), cov10 = NA, cov25 = NA, covall = NA,
      types1 = NA, types3 = NA, precision = NA, need50 = need(P$rankings[[m]], .5), need80 = need(P$rankings[[m]], .8))
    cat(sprintf("%-14s %-9s k=%2d  types %2d  top-10 list %4d genes (in gate %3.0f%%)  all %5d\n",
                d, L, P$k_true, length(types), length(M10), 100 * mean(M10 %in% gate), length(Mall)))
  }
}
fill <- function(x) { for (v in c("need50", "need80")) if (!v %in% names(x)) x[[v]] <- NA; x }
R <- do.call(rbind, lapply(rows, fill))
utils::write.csv(R, file.path(out_dir, "marker_coverage_per_dataset.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, refs_out), file.path(out_dir, "reference_markers.csv"), row.names = FALSE)

C <- R[!is.na(R$budget), ]
missing_difs <- setdiff(paste(rep(unique(C$dataset), each = length(BUDGETS)), BUDGETS),
                        paste(C$dataset, C$budget)[C$method == "DIFS"])
if (length(missing_difs)) cat("\nDIFS result missing for:", paste(missing_difs, collapse = "; "), "\n")
summ <- function(sub, label) {
  a <- stats::aggregate(cbind(cov10, cov25, covall, types1, types3, precision) ~ list + method + budget, sub, mean)
  a$datasets <- label; a
}
ind <- C[C$dataset %in% INDEPENDENT, ]
S <- rbind(summ(ind, "independent9"), summ(C, "all"))
utils::write.csv(S, file.path(out_dir, "marker_coverage_summary.csv"), row.names = FALSE)

tests <- list()
for (L in unique(C$list)) for (lab in c("independent9", "all"))
  for (cmp in c("Seurat", "FEAST", "Seurat_gated", "DIFS_stage1only"))
  for (b in BUDGETS) for (v in c("cov10", "types1", "precision")) {
    sub <- if (lab == "all") C else ind; sub <- sub[sub$list == L, ]
    a <- sub[sub$method == "DIFS" & sub$budget == b, c("dataset", v)]
    z <- sub[sub$method == cmp & sub$budget == b, c("dataset", v)]
    m <- merge(a, z, by = "dataset"); if (nrow(m) < 2) next
    dd <- round(m[[2]] - m[[3]], 12)     # same rounding as signrank_exact
    tests[[length(tests) + 1]] <- data.frame(list = L, datasets = lab, vs = cmp, budget = b, metric = v,
      n = nrow(m), mean_DIFS = mean(m[[2]]), mean_other = mean(m[[3]]), mean_diff = mean(dd),
      wins = sum(dd > 0), losses = sum(dd < 0), p_exact = signrank_exact(dd),
      primary = L == "hurdle" && lab == "independent9" && cmp %in% c("Seurat", "FEAST") && b <= 300)
  }
Tt <- do.call(rbind, tests)
utils::write.csv(Tt, file.path(out_dir, "marker_coverage_tests.csv"), row.names = FALSE)

options(width = 200)
for (L in intersect(LISTS, unique(S$list))) for (v in c("cov10", "types1", "precision")) {
  cat(sprintf("\n== %s list: %s, mean over the nine independent datasets ==\n", L, v))
  w <- stats::reshape(S[S$datasets == "independent9" & S$list == L, c("method", "budget", v)],
                      idvar = "method", timevar = "budget", direction = "wide")
  names(w) <- sub(paste0("^", v, "\\."), "", names(w)); print(w, row.names = FALSE, digits = 3)
}
cat("\n== PRIMARY tests: hurdle list, DIFS minus comparator, nine independent datasets ==\n")
print(Tt[Tt$primary, c("vs", "budget", "metric", "n", "mean_DIFS", "mean_other", "mean_diff",
                       "wins", "losses", "p_exact")], row.names = FALSE, digits = 3)
cat("\n== secondary: level / detection / wilcox lists, cov10 vs Seurat and FEAST, b <= 300 ==\n")
sec <- Tt[Tt$list != "hurdle" & Tt$datasets == "independent9" & Tt$vs %in% c("Seurat", "FEAST") &
          Tt$budget <= 300 & Tt$metric == "cov10", ]
print(sec[, c("list", "vs", "budget", "n", "mean_DIFS", "mean_other", "wins", "losses", "p_exact")],
      row.names = FALSE, digits = 3)
nd <- R[is.na(R$budget) & R$list == "hurdle", c("dataset", "method", "need50")]
cat("\n== genes needed to cover 50% of the hurdle top-10 list (full rankings) ==\n")
print(stats::reshape(nd, idvar = "dataset", timevar = "method", direction = "wide"), row.names = FALSE)
