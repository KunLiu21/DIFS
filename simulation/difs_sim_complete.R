# difs_sim_complete.R: revision-2 analysis source.
# Source: code/as_run/08_simulation/difs_sim_complete.R
# Usage, scope and limitations: benchmark/REVISION2.md.
difs_args <- function(args, known, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    if (!kv[2] %in% known) stop("unknown argument '", kv[2], "'. Accepts: ", paste(known, collapse = ", "))
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
}
task <- "run"; k <- NA; out_dir <- "../results/sim_complete"; seed <- 1
difs_args(commandArgs(TRUE), c("task", "k", "out_dir", "seed"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
BUDGETS <- c(100, 200, 500); K <- 5L; NC <- 600L
PROPS <- c(.45, .25, .15, .10, .05)

TT <- rbind(
  data.frame(regime = "E8",   level = rep(c(2, 4, 6, 8), each = 5), r = rep(1:5, 4)),
  data.frame(regime = "E4nd", level = rep(1:5, each = 5),          r = rep(1:5, 5)),
  data.frame(regime = "E4",   level = rep(c(2, 4, 6, 8), each = 5), r = rep(1:5, 4)),
  data.frame(regime = "E8md", level = rep(c(2, 4, 6, 8), each = 5), r = rep(1:5, 4)))
TT$index <- seq_len(nrow(TT))
if (task == "list") { print(TT, row.names = FALSE); quit(save = "no") }

gen_E8 <- function(effect, seed, md = FALSE) {           # difs_sim_E8_assumption.R
  set.seed(seed); N <- 3000
  type <- sample(seq_along(PROPS), NC, replace = TRUE, prob = PROPS)
  pis <- tabulate(type, length(PROPS)) / NC
  cls <- c(rep("marker", 240), rep("hv_null", 240), rep("null", N - 480))
  sg <- stats::runif(N, .4, .6); mu <- stats::runif(N, 2.5, 4); det <- stats::runif(N, .3, .9)
  hi <- sample(seq_along(PROPS), N, replace = TRUE); delta <- effect * sg
  sd_h <- sqrt(sg^2 + pis[hi] * (1 - pis[hi]) * delta^2)
  X <- matrix(0, N, NC)
  for (g in seq_len(N)) {
    x <- switch(cls[g],
      marker  = mu[g] + delta[g] * (type == hi[g]) + stats::rnorm(NC, 0, sg[g]),
      hv_null = mu[g] + delta[g] * pis[hi[g]] + stats::rnorm(NC, 0, sd_h[g]),
      null    = mu[g] + delta[g] * pis[hi[g]] + stats::rnorm(NC, 0, sg[g]))
    m <- if (cls[g] == "marker") mu[g] + delta[g] * (type == hi[g]) else rep(mu[g] + delta[g] * pis[hi[g]], NC)
    pdet <- if (md) stats::plogis(2.03 * (m - 2.917)) else det[g]
    X[g, ] <- ifelse(stats::runif(NC) < pdet, x, -Inf)
  }
  sf <- stats::rlnorm(NC, 0, .3)
  lam <- sweep(expm1(pmax(X, 0)), 2, sf, "*"); lam[!is.finite(X)] <- 0
  cnt <- matrix(stats::rpois(length(lam), lam), N, NC)
  list(counts = cnt, marker = cls == "marker", hv = cls == "hv_null", scale = stats::median(colSums(cnt)))
}
gen_E4nd <- function(lf, seed, disp = 0.2, depth = 20000, nf = 1200, nb = 10000, L = 6) {
  set.seed(seed)                                         # difs_sim_E4nd_nondropout.R
  lab <- sample(1:5, NC, TRUE, PROPS); libs <- depth * stats::rlnorm(NC, 0, .35)
  cp <- stats::rlnorm(nf, log(L), 0.4); rest <- max(1e4 - sum(cp), 1e3)
  cb <- stats::rlnorm(nb, 0, 1.1); cb <- cb / sum(cb) * rest
  share <- c(cp, cb) / 1e4; ninf <- round(.08 * nf)
  inf <- c(sample(c(rep(TRUE, ninf), rep(FALSE, nf - ninf))), rep(FALSE, nb))
  up <- sample(1:5, nf + nb, TRUE)
  S <- matrix(share, nf + nb, NC); for (g in which(inf)) S[g, lab == up[g]] <- share[g] * 2^lf
  S <- sweep(S, 2, colSums(S), "/"); MU <- sweep(S, 2, libs, "*")
  cnt <- matrix(stats::rnbinom(length(MU), mu = MU, size = 1 / disp), nrow(MU))
  cnt <- cnt * matrix(stats::rbinom(length(MU), 1, 1 - 1 / (1 + exp(log(MU) - 1))), nrow(MU))
  list(counts = cnt, marker = inf, hv = rep(FALSE, length(inf)), scale = 1e4)
}
gen_E4 <- function(lf, seed) {                           # E4 / E4b generator, dispersion 0.2
  s <- sim_panel(n_cells = NC, n_genes = 3000, type_props = PROPS, prop_informative = 0.08,
                 log2fc = lf, disp = 0.2, mu_meanlog = log(10), n_bg_genes = NULL, seed = seed)
  list(counts = s$counts, marker = s$informative, hv = rep(FALSE, nrow(s$counts)), scale = 1e4)
}
make_panel <- function(i) {
  p <- TT[i, ]; sd0 <- 9100000 + i * 101 + p$r
  P <- switch(p$regime, E8 = gen_E8(p$level, sd0), E4nd = gen_E4nd(p$level, sd0), E4 = gen_E4(p$level, sd0),
              E8md = gen_E8(p$level, sd0, md = TRUE))
  rownames(P$counts) <- sprintf("g%05d", seq_len(nrow(P$counts)))
  colnames(P$counts) <- sprintf("c%04d", seq_len(ncol(P$counts)))
  P
}

run_panel <- function(i, budgets = BUDGETS) {
  t0 <- proc.time()[["elapsed"]]
  P <- make_panel(i); g <- rownames(P$counts)
  seu <- CreateSeuratObject(counts = P$counts)
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = P$scale, verbose = FALSE)
  lm_ <- as.matrix(GetAssayData(seu, slot = "data")); cnt <- P$counts
  gate <- difs_gate_genes(lm_, log(5), gate_rule = "submitted", gate_k = K)
  safe <- function(expr, what) tryCatch(expr, error = function(e) {
    message("  ", what, " failed: ", conditionMessage(e)); NULL })
  set.seed(seed)
  rk <- list(
    DIFS_stage1only = safe(difs_stage1_ranking(lm_, min.expression = log(5), gate_rule = "submitted",
                                               gate_k = K, key = "normal_lookup", seed = seed), "stage I"),
    Seurat = safe(VariableFeatures(FindVariableFeatures(seu, selection.method = "vst",
                                                        nfeatures = nrow(seu), verbose = FALSE)), "Seurat"),
    FEAST = NULL)
  feast_seed <- NA_integer_
  for (s in seed + 0:2) {
    set.seed(s)
    r <- safe(feature_ranking("FEAST", data.seu = seu, data.log.mat = lm_, data.mat = cnt, cds = NULL,
                              cluster_count = K, max_n = nrow(lm_), seed = s), sprintf("FEAST seed %d", s))
    if (!is.null(r)) { rk$FEAST <- r; feast_seed <- s; break }
  }
  difs <- lapply(budgets, function(b) {
    set.seed(seed); WHICHPC_FALLBACK <<- 0L
    r <- safe(difs_features(data.log.mat = lm_, data.seu = seu, cluster_count = K,
                            binary.bound = binary.bound.decider(lm_), min.expression = log(5),
                            gate_rule = "submitted", gate_k = K, tune_stage1 = TRUE, tune_ratio = TRUE,
                            stage1_key = "normal_lookup", use_stage2 = TRUE, seed = seed,
                            tops = DIFS_TOPS_GRID, n_fixed = b), sprintf("DIFS b=%d", b))
    if (is.null(r)) return(NULL)
    list(budget = b, genes = intersect(r$markers, g), stage1_n = r$stage1_n,
         stage2_n = r$stage2_n, ratio = r$ratio, whichpc_fallbacks = WHICHPC_FALLBACK)
  })
  out <- list(index = i, panel = TT[i, ], genes = g, marker = g[P$marker], hv = g[P$hv],
              gate = gate, rankings = rk, difs = difs, budgets = budgets, feast_seed = feast_seed,
              elapsed_sec = proc.time()[["elapsed"]] - t0)
  saveRDS(out, file.path(out_dir, sprintf("panel__%03d.rds", i)))
  cat(sprintf("panel %d (%s, level %g, rep %d): gate %d, markers %d (in gate %d), %.0f s\n", i,
              TT$regime[i], TT$level[i], TT$r[i], length(gate), sum(P$marker),
              sum(g[P$marker] %in% gate), out$elapsed_sec))
}

summarise <- function(dir) {
  fs <- list.files(dir, pattern = "^panel__.*\\.rds$", full.names = TRUE)
  if (!length(fs)) stop("no panel__*.rds in ", dir)
  rows <- list()
  for (f in fs) {
    X <- readRDS(f); M <- X$marker; H <- X$hv; gate <- X$gate
    add <- function(m, b, set, extra = list()) rows[[length(rows) + 1]] <<- data.frame(
      X$panel, method = m, budget = b, n_selected = length(set),
      recall = mean(M %in% set), precision = if (length(set)) mean(set %in% M) else NA,
      hv_selected = sum(set %in% H), from_stage2_only = if (is.null(extra$s2)) NA else extra$s2,
      markers_in_gate = mean(M %in% gate), row.names = NULL)
    for (b in X$budgets) {
      for (m in names(X$rankings)) if (!is.null(X$rankings[[m]])) add(m, b, utils::head(X$rankings[[m]], b))
      d <- Filter(function(z) !is.null(z) && z$budget == b, X$difs)
      if (length(d)) {
        d <- d[[1]]; s1 <- utils::head(X$rankings$DIFS_stage1only, d$stage1_n)
        add("DIFS", b, d$genes, list(s2 = sum(d$genes %in% M & !(d$genes %in% s1))))
      }
      n <- length(gate); mg <- sum(M %in% gate)
      rows[[length(rows) + 1]] <- data.frame(X$panel, method = "GateOnly", budget = b,
        n_selected = min(b, n), recall = mg * min(b, n) / n / length(M), precision = mg / n,
        hv_selected = sum(H %in% gate) * min(b, n) / n, from_stage2_only = NA,
        markers_in_gate = mg / length(M), row.names = NULL)
    }
  }
  R <- do.call(rbind, rows)
  utils::write.csv(R, file.path(dir, "sim_complete_per_panel.csv"), row.names = FALSE)
  A <- stats::aggregate(cbind(recall, precision, hv_selected, from_stage2_only, markers_in_gate) ~
                          regime + level + method + budget, R, mean, na.action = stats::na.pass)
  utils::write.csv(A, file.path(dir, "sim_complete_summary.csv"), row.names = FALSE)
  options(width = 200)
  for (rg in unique(R$regime)) for (b in unique(R$budget)) for (v in c("recall", "precision")) {
    s <- A[A$regime == rg & A$budget == b, ]
    w <- stats::reshape(s[, c("level", "method", v)], idvar = "method", timevar = "level", direction = "wide")
    names(w) <- sub(paste0("^", v, "\\."), "", names(w))
    cat(sprintf("\n== %s, %s at %d features, by level ==\n", rg, v, b)); print(w, row.names = FALSE, digits = 3)
  }
  fs_ok <- vapply(fs, function(f) { x <- readRDS(f); if (is.null(x$feast_seed)) NA_real_ else x$feast_seed }, numeric(1))
  cat("\nFEAST: succeeded on", sum(!is.na(fs_ok)), "of", length(fs_ok), "panels;",
      sum(fs_ok > seed, na.rm = TRUE), "needed a second or third seed\n")
  cat("\n== DIFS vs competitors, panel wins on recall (budget 100 / 200) ==\n")
  for (cmp in c("Seurat", "FEAST", "DIFS_stage1only")) for (b in c(100, 200)) {
    a <- R[R$method == "DIFS" & R$budget == b, ]; z <- R[R$method == cmp & R$budget == b, ]
    m <- merge(a, z, by = c("regime", "level", "r"))
    if (!nrow(m)) next
    tab <- stats::aggregate(cbind(win = recall.x > recall.y, loss = recall.x < recall.y) ~ regime, m, sum)
    cat(sprintf("vs %-16s b=%d: %s\n", cmp, b,
                paste(sprintf("%s %d-%d", tab$regime, tab$win, tab$loss), collapse = "; ")))
  }
  dd <- A[A$method == "DIFS", c("regime", "level", "budget", "from_stage2_only")]
  cat("\n== DIFS: markers selected that came from stage II only (mean per panel) ==\n")
  print(stats::reshape(dd, idvar = c("regime", "level"), timevar = "budget", direction = "wide"), row.names = FALSE)
}

if (task == "summary") { summarise(out_dir); quit(save = "no") }

here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) "sim")
source("Functions_controlled.R")
if (!requireNamespace("diptest", quietly = TRUE))
  stop("diptest not found on .libPaths(): ", paste(.libPaths(), collapse = " : "),
       "\n  install diptest in your R library or set DIFS_LIB to its library directory")
core <- new.env(); sys.source(file.path(here, "difs_sim_core.R"), envir = core)
sim_panel <- core$sim_panel                             # E4 generator
if (!file.exists("difs_sc3_complete.R")) stop("difs_sc3_complete.R is required (SC3 repair)")
source("difs_sc3_complete.R")
suppressPackageStartupMessages(library(Seurat))
WHICHPC_FALLBACK <- 0L
whichPC <- function(seurat) {
  pct <- seurat[["pca"]]@stdev / sum(seurat[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  if (is.na(co2)) WHICHPC_FALLBACK <<- WHICHPC_FALLBACK + 1L
  pcs <- suppressWarnings(min(co1, co2, na.rm = TRUE))
  if (!is.finite(pcs)) pcs <- length(pct)
  pcs
}
if (task == "smoke") {
  out_dir <- file.path(out_dir, "smoke"); dir.create(out_dir, showWarnings = FALSE)
  for (i in c(match("E8", TT$regime) + 15, match("E4nd", TT$regime) + 15, match("E4", TT$regime) + 5))
    run_panel(i, budgets = 100)
  summarise(out_dir)
} else {
  if (is.na(k) || k < 1 || k > nrow(TT)) stop("k= must be 1..", nrow(TT), " (task=list)")
  run_panel(as.integer(k))
}
