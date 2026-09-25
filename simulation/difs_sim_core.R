
suppressPackageStartupMessages(library(diptest))


dip_stat <- function(x) diptest::dip(sort(x[is.finite(x)]))

dip_null_sample <- function(n, B = 2000, ref = c("uniform", "normal", "tnormal"),
                            trunc_q = 0.5) {
  ref <- match.arg(ref)
  gen <- switch(ref,
    uniform = function() stats::runif(n),
    normal  = function() stats::rnorm(n),
    tnormal = function() {
      lo <- stats::qnorm(trunc_q)
      stats::qnorm(stats::runif(n, trunc_q, 1))
    })
  vapply(seq_len(B), function(i) dip_stat(gen()), numeric(1))
}

dip_pvalue_mc <- function(D, null, add_one = TRUE) {
  if (add_one) (1 + sum(null >= D)) / (length(null) + 1) else mean(null >= D)
}


DEFAULT_N_GRID <- c(20, 25, 30, 35, 40, 50, 60, 75, 90, 110, 130, 160, 200,
                    250, 300, 400, 500, 650, 800, 1000, 1300, 1600, 2000,
                    2500, 3000, 4000, 5000, 7000, 10000)

build_dip_null_table <- function(n_grid = DEFAULT_N_GRID, B = 5000,
                                 ref = "normal", trunc_q = 0.5,
                                 probs = c(seq(0.0005, 0.9990, length.out = 700),
                                           0.9995, 0.9999),
                                 verbose = TRUE) {
  Q <- matrix(NA_real_, nrow = length(n_grid), ncol = length(probs),
              dimnames = list(as.character(n_grid), NULL))
  for (i in seq_along(n_grid)) {
    if (verbose) message(sprintf("  null table [%s] n = %d (%d/%d)",
                                 ref, n_grid[i], i, length(n_grid)))
    d <- dip_null_sample(n_grid[i], B = B, ref = ref, trunc_q = trunc_q)
    Q[i, ] <- stats::quantile(d, probs = probs, names = FALSE, type = 7)
  }
  structure(list(n_grid = n_grid, probs = probs, Q = Q, ref = ref, B = B),
            class = "dip_null_table")
}

dip_pvalue_lookup <- function(D, n, tab) {
  dip_pvalue_lookup_vec(D, n, tab)
}

dip_pvalue_lookup_vec <- function(D, n, tab) {
  stopifnot(inherits(tab, "dip_null_table"))
  lg    <- log(tab$n_grid)
  upper <- 1 - tab$probs
  logQ  <- log(tab$Q)
  n_cl  <- pmin(pmax(n, min(tab$n_grid)), max(tab$n_grid))
  out   <- numeric(length(D))
  cache <- new.env(parent = emptyenv())
  qcurve_for <- function(nn) {
    key <- format(nn, digits = 10)
    if (!is.null(cache[[key]])) return(cache[[key]])
    ln <- log(nn)
    i2 <- findInterval(ln, lg, rightmost.closed = TRUE)
    i2 <- min(max(i2, 1L), length(lg) - 1L); i1 <- i2; i2 <- i2 + 1L
    w  <- if (lg[i2] > lg[i1]) (ln - lg[i1]) / (lg[i2] - lg[i1]) else 0
    q  <- cummax(exp((1 - w) * logQ[i1, ] + w * logQ[i2, ]))
    k   <- length(q); idx <- (k - 39):k
    cf  <- stats::coef(stats::lm(log(upper[idx]) ~ q[idx]))
    val <- list(q = q, cf = cf)
    assign(key, val, envir = cache)
    val
  }
  for (i in seq_along(D)) {
    if (!is.finite(D[i])) { out[i] <- NA_real_; next }
    v <- qcurve_for(n_cl[i]); q <- v$q
    if (D[i] <= q[1]) { out[i] <- 1; next }
    if (D[i] >= q[length(q)]) {
      out[i] <- max(.Machine$double.xmin,
                    exp(v$cf[1] + v$cf[2] * D[i])); next
    }
    out[i] <- stats::approx(q, upper, xout = D[i], rule = 2)$y
  }
  out
}


DIFS_GATE_AS_IMPLEMENTED <- log(5)
DIFS_GATE_AS_DOCUMENTED  <- 0.5

difs_stage1_gene <- function(x, n_cells, method = c("lookup", "mc", "hartigan"),
                             ref = "normal", tab = NULL, B = 2000,
                             min.expression = 0.5,
                             min.expression.proportion = 0.05,
                             min.cells = 35, add_one = TRUE,
                             jitter_sd = 0, min.distinct.frac = 0) {
  method <- match.arg(method)
  if (!length(x) || max(x, na.rm = TRUE) <= min.expression)
    return(c(p = 1.1, D = NA_real_, n = 0))
  xe <- x[x > min.expression]
  if (length(xe) < max(min.expression.proportion * n_cells, min.cells))
    return(c(p = 1.1, D = NA_real_, n = length(xe)))
  if (min.distinct.frac > 0 &&
      length(unique(xe)) < min.distinct.frac * length(xe))
    return(c(p = 1.1, D = NA_real_, n = length(xe)))
  if (jitter_sd > 0) xe <- xe + stats::rnorm(length(xe), 0, jitter_sd * stats::sd(xe))
  D <- dip_stat(xe)
  p <- switch(method,
    hartigan = unname(diptest::dip.test(xe)$p.value),
    mc       = dip_pvalue_mc(D, dip_null_sample(length(xe), B, ref), add_one),
    lookup   = dip_pvalue_lookup(D, length(xe), tab))
  c(p = p, D = D, n = length(xe))
}

difs_stage1_matrix <- function(mat, ...) {
  n_cells <- ncol(mat)
  res <- t(vapply(seq_len(nrow(mat)),
                  function(i) difs_stage1_gene(mat[i, ], n_cells, ...),
                  numeric(3)))
  rownames(res) <- rownames(mat)
  as.data.frame(res)
}


difs_stage2 <- function(mat, clusters, binary.bound = log(5),
                        min.mean.expr = 0.1) {
  bin  <- (mat > binary.bound) * 1L
  cl   <- as.factor(clusters)
  rate <- vapply(levels(cl), function(k) rowMeans(bin[, cl == k, drop = FALSE]),
                 numeric(nrow(bin)))
  if (is.null(dim(rate))) rate <- matrix(rate, nrow = nrow(bin))
  best <- max.col(rate, ties.method = "first")
  p <- vapply(seq_len(nrow(bin)), function(i) {
    if (max(rate[i, ]) < min.mean.expr) return(1.1)
    inn <- bin[i, cl == levels(cl)[best[i]]]
    out <- bin[i, cl != levels(cl)[best[i]]]
    tabm <- matrix(c(sum(inn == 0), sum(out == 0),
                     sum(inn == 1), sum(out == 1)), nrow = 2)
    stats::fisher.test(tabm)$p.value
  }, numeric(1))
  data.frame(gene = rownames(mat), p = p, cluster = best,
             stringsAsFactors = FALSE)
}

difs_interleave <- function(vec1, vec2, ratio) {
  loops <- sum(ratio); n <- max(length(vec1), length(vec2))
  v1 <- rep(vec1, length.out = n * ratio[1])
  v2 <- rep(vec2, length.out = n * ratio[2])
  out <- character(0)
  for (b in seq_len(n)) {
    out <- c(out,
             v1[((b - 1) * ratio[1] + 1):(b * ratio[1])],
             v2[((b - 1) * ratio[2] + 1):(b * ratio[2])])
  }
  unique(out[!is.na(out)])
}


lognormalize <- function(counts, scale.factor = 1e4) {
  ls <- pmax(colSums(counts), 1)
  log1p(sweep(counts, 2, ls, "/") * scale.factor)
}

normalize_counts <- function(counts, method = c("lognorm", "sctransform"),
                             scale.factor = 1e4) {
  method <- match.arg(method)
  if (method == "lognorm") return(lognormalize(counts, scale.factor))
  if (!requireNamespace("Seurat", quietly = TRUE))
    stop("method = 'sctransform' requires the Seurat package")
  seu <- Seurat::CreateSeuratObject(counts)
  seu <- Seurat::SCTransform(seu, verbose = FALSE)
  as.matrix(Seurat::GetAssayData(seu, assay = "SCT", slot = "data"))
}

CONTINUOUS_NULL_TYPES <- c("normal", "lognormal", "gamma")
COUNT_NULL_TYPES      <- c("nb", "zinb", "zinb_high")

sim_null_gene <- function(type, n, mu = 3, disp = 0.5, size_factor = NULL,
                          dropout_mid = 1.0, dropout_shape = -1.0) {
  if (is.null(size_factor)) size_factor <- rep(1, n)
  switch(type,
    normal    = pmax(stats::rnorm(n, mean = log1p(mu), sd = 0.6), 0),
    lognormal = stats::rlnorm(n, meanlog = log(log1p(mu)), sdlog = 0.5),
    gamma     = stats::rgamma(n, shape = 2, scale = log1p(mu) / 2),
    nb        = stats::rnbinom(n, mu = mu * size_factor, size = 1 / disp),
    zinb      = {
      y <- stats::rnbinom(n, mu = mu * size_factor, size = 1 / disp)
      pd <- 1 / (1 + exp(-dropout_shape * (log(mu * size_factor) - dropout_mid)))
      y * stats::rbinom(n, 1, 1 - pd)
    },
    zinb_high = {
      y <- stats::rnbinom(n, mu = mu * size_factor, size = 1 / disp)
      pd <- 1 / (1 + exp(-dropout_shape * (log(mu * size_factor) - (dropout_mid + 1.5))))
      y * stats::rbinom(n, 1, 1 - pd)
    },
    stop("unknown null type: ", type))
}

sim_alt_gene <- function(type, n, pi, log2fc, mu = 3, disp = 0.5,
                         size_factor = NULL, dropout_mid = 1.0,
                         dropout_shape = -1.0) {
  if (is.null(size_factor)) size_factor <- rep(1, n)
  n2 <- max(1L, round(n * pi))
  grp <- rep(0L, n); grp[sample.int(n, n2)] <- 1L
  mu_c <- mu * (2 ^ (log2fc * grp))
  switch(type,
    normal = pmax(stats::rnorm(n, mean = log1p(mu_c), sd = 0.6), 0),
    zinb   = {
      y <- stats::rnbinom(n, mu = mu_c * size_factor, size = 1 / disp)
      pd <- 1 / (1 + exp(-dropout_shape * (log(mu_c * size_factor) - dropout_mid)))
      y * stats::rbinom(n, 1, 1 - pd)
    },
    nb     = stats::rnbinom(n, mu = mu_c * size_factor, size = 1 / disp),
    stop("unknown alt type: ", type))
}

sim_panel <- function(n_cells = 600, n_genes = 4000, type_props = c(.45,.25,.15,.10,.05),
                      prop_informative = 0.08, log2fc = 3, disp = 0.5,
                      mu_meanlog = log(2), mu_sdlog = 1.1,
                      n_bg_genes = NULL,
                      dropout_mid = 1.0, dropout_shape = -1.0,
                      libsize_sd = 0.35, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  type_props <- type_props / sum(type_props)
  labels <- factor(sample(seq_along(type_props), n_cells, TRUE, type_props))
  sf     <- stats::rlnorm(n_cells, 0, libsize_sd)
  mu_g   <- stats::rlnorm(n_genes, meanlog = mu_meanlog, sdlog = mu_sdlog)
  n_inf  <- round(n_genes * prop_informative)
  inf_id <- sort(sample(n_genes, n_inf))
  up_in  <- sample(seq_along(type_props), n_inf, TRUE)
  MU <- matrix(mu_g, n_genes, n_cells)
  for (j in seq_len(n_inf)) {
    g <- inf_id[j]
    MU[g, labels == up_in[j]] <- mu_g[g] * 2 ^ log2fc
  }
  MU <- sweep(MU, 2, sf, "*")
  counts <- matrix(stats::rnbinom(length(MU), mu = MU, size = 1 / disp),
                   nrow = n_genes)
  pd <- 1 / (1 + exp(-dropout_shape * (log(MU) - dropout_mid)))
  counts <- counts * matrix(stats::rbinom(length(MU), 1, 1 - pd), nrow = n_genes)
  if (is.null(n_bg_genes)) n_bg_genes <- max(0L, 6000L - n_genes)
  if (n_bg_genes > 0) {
    mu_b <- stats::rlnorm(n_bg_genes, meanlog = mu_meanlog, sdlog = mu_sdlog)
    MUb  <- sweep(matrix(mu_b, n_bg_genes, n_cells), 2, sf, "*")
    cb   <- matrix(stats::rnbinom(length(MUb), mu = MUb, size = 1 / disp),
                   nrow = n_bg_genes)
    pdb  <- 1 / (1 + exp(-dropout_shape * (log(MUb) - dropout_mid)))
    cb   <- cb * matrix(stats::rbinom(length(MUb), 1, 1 - pdb), nrow = n_bg_genes)
    counts <- rbind(counts, cb)
    inf_flag <- c(seq_len(n_genes) %in% inf_id, rep(FALSE, n_bg_genes))
  } else {
    inf_flag <- seq_len(n_genes) %in% inf_id
  }
  rownames(counts) <- sprintf("gene%05d", seq_len(nrow(counts)))
  colnames(counts) <- sprintf("cell%05d", seq_len(n_cells))
  up_vec <- rep(NA_integer_, nrow(counts)); up_vec[inf_id] <- up_in
  names(up_vec) <- rownames(counts)
  list(counts = counts, labels = labels, informative = inf_flag,
       up_in = up_vec, type_props = type_props,
       n_focal = n_genes, n_background = nrow(counts) - n_genes)
}


rank_seurat_vst <- function(counts, clip = NULL) {
  m <- rowMeans(counts); v <- apply(counts, 1, stats::var)
  keep <- m > 0 & v > 0
  fit  <- stats::loess(log10(v[keep]) ~ log10(m[keep]), span = 0.3)
  exp_sd <- rep(NA_real_, length(m))
  exp_sd[keep] <- sqrt(10 ^ stats::predict(fit))
  if (is.null(clip)) clip <- sqrt(ncol(counts))
  z <- sweep(counts, 1, m, "-"); z <- sweep(z, 1, exp_sd, "/")
  z[!is.finite(z)] <- 0; z[z > clip] <- clip
  sv <- apply(z, 1, stats::var); sv[!keep] <- NA_real_
  sv
}

rank_variance    <- function(logmat) apply(logmat, 1, stats::var)
rank_dropoutrate <- function(counts) rowMeans(counts == 0)
rank_oracle_F <- function(logmat, labels) {
  labels <- as.factor(labels); K <- nlevels(labels); N <- ncol(logmat)
  grand <- rowMeans(logmat)
  ssb <- matrix(0, nrow(logmat), 1); ssw <- matrix(0, nrow(logmat), 1)
  for (k in levels(labels)) {
    sub <- logmat[, labels == k, drop = FALSE]; nk <- ncol(sub)
    mk  <- rowMeans(sub)
    ssb <- ssb + nk * (mk - grand) ^ 2
    ssw <- ssw + rowSums((sub - mk) ^ 2)
  }
  f <- (ssb / (K - 1)) / (ssw / (N - K))
  as.vector(ifelse(is.finite(f), f, 0))
}


auc_score <- function(score, truth) {
  ok <- is.finite(score); score <- score[ok]; truth <- truth[ok]
  if (!any(truth) || all(truth)) return(NA_real_)
  r <- rank(score); n1 <- sum(truth); n0 <- sum(!truth)
  (sum(r[truth]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

precision_at_n <- function(order_idx, truth, N) mean(truth[utils::head(order_idx, N)])

genes_needed_for_recall <- function(order_idx, truth, recall = 0.5) {
  hits <- cumsum(truth[order_idx]); target <- ceiling(recall * sum(truth))
  w <- which(hits >= target)
  if (!length(w)) NA_integer_ else w[1]
}

enrichment_at_n <- function(order_idx, truth, N)
  precision_at_n(order_idx, truth, N) / mean(truth)

order_by <- function(score, decreasing) {
  score[!is.finite(score)] <- if (decreasing) -Inf else Inf
  order(score, decreasing = decreasing)
}

difs_discreteness_diagnostic <- function(logmat, min.expression = 0.5,
                                         min.expression.proportion = 0.05,
                                         min.cells = 35,
                                         risky_distinct_frac = 0.9) {
  n_cells <- ncol(logmat)
  res <- t(vapply(seq_len(nrow(logmat)), function(i) {
    x  <- logmat[i, ]
    xe <- x[x > min.expression]
    c(n_expressing = length(xe),
      n_distinct   = length(unique(xe)),
      mean_expr    = if (length(xe)) mean(xe) else 0)
  }, numeric(3)))
  res <- as.data.frame(res); rownames(res) <- rownames(logmat)
  res$tested <- res$n_expressing >= pmax(min.expression.proportion * n_cells,
                                         min.cells)
  res$distinct_frac <- res$n_distinct / pmax(res$n_expressing, 1)
  res$risky <- res$tested & res$distinct_frac < risky_distinct_frac
  res
}
