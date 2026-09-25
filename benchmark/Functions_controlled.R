
DIFS_OPTIONAL_PKGS <- c("scran", "glmGamPoi")

difs_fix_libpaths <- function(lib     = Sys.getenv("DIFS_LIB", ""),
                              damaged = Sys.getenv("DIFS_EXCLUDE_LIB", "")) {
  lp <- .libPaths()
  if (nzchar(damaged)) {
    d  <- normalizePath(damaged, mustWork = FALSE)
    lp <- lp[normalizePath(lp, mustWork = FALSE) != d]
  }
  if (nzchar(lib) && dir.exists(lib)) lp <- unique(c(lib, lp))
  .libPaths(lp)
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = ":"))
  invisible(.libPaths())
}
difs_fix_libpaths()

difs_source_base <- function(path = "Functions_mixture_auto.R",
                             optional = DIFS_OPTIONAL_PKGS) {
  txt <- readLines(path, warn = FALSE)
  for (p in optional) {
    pat <- sprintf("^(\\s*)library\\(%s\\)\\s*$", p)
    rep <- sprintf("\\1if (!inherits(try(suppressPackageStartupMessages(library(%s)), silent = TRUE), \"try-error\")) invisible(TRUE) else message(\"note: %s unavailable, continuing without it\")", p, p)
    txt <- sub(pat, rep, txt)
  }
  eval(parse(text = txt), envir = globalenv())
  invisible(TRUE)
}
difs_source_base()

difs_nb_family <- function() {
  for (ns in c("VGAM", "monocle")) {
    f <- tryCatch(getExportedValue(ns, "negbinomial.size"), error = function(e) NULL)
    if (is.function(f)) return(f())
  }
  if (exists("negbinomial.size")) return(get("negbinomial.size")())
  stop("cannot find negbinomial.size(); VGAM does not seem to be available")
}

difs_adf <- function(df) {
  if (requireNamespace("Biobase", quietly = TRUE))
    return(Biobase::AnnotatedDataFrame(data = df))
  methods::new("AnnotatedDataFrame", data = df)
}

num_features_decider <- function(data.seu, cluster_count,
                                 tops = seq(500, 2500, 100), top_genes) {
  top_genes <- top_genes[!is.na(top_genes)]
  n_avail <- length(top_genes)
  if (n_avail < 1L) stop("num_features_decider: no usable genes in top_genes")
  tops <- sort(unique(pmin(as.integer(tops), n_avail)))
  tops <- tops[tops >= 1L]
  data.log.mat <- as.matrix(GetAssayData(data.seu, slot = "data"))
  mse <- numeric(length(tops)); sets <- vector("list", length(tops))
  for (i in seq_along(tops)) {
    markers <- top_genes[seq_len(tops[i])]
    data.seu <- seudo_clustering(data.seu, markers = markers,
                                 cluster_count = cluster_count,
                                 clustering.method = "Refined Louvain")
    sets[[i]] <- markers
    mse[i] <- cal_MSE(data.log.mat, data.seu$seurat_clusters)
  }
  names(mse) <- tops
  best <- which.min(mse)
  list(markers = sets[[best]], mse = mse[best],
       n_available = n_avail, clamped = any(as.integer(tops) == n_avail))
}

sce_to_monocle <- function(sce) {
  if (requireNamespace("scran", quietly = TRUE) &&
      "convertTo" %in% getNamespaceExports("scran")) {
    ok <- try(suppressWarnings(scran::convertTo(sce, type = "monocle")),
              silent = TRUE)
    if (!inherits(ok, "try-error")) return(ok)
  }
  cnt <- as.matrix(SummarizedExperiment::assay(sce, "counts"))
  pd  <- as.data.frame(SummarizedExperiment::colData(sce))
  if (!ncol(pd)) pd <- data.frame(row.names = colnames(cnt))
  rownames(pd) <- colnames(cnt)
  fd  <- data.frame(gene_short_name = rownames(cnt), row.names = rownames(cnt),
                    stringsAsFactors = FALSE)
  monocle::newCellDataSet(
    cnt,
    phenoData        = difs_adf(pd),
    featureData      = difs_adf(fd),
    expressionFamily = difs_nb_family())
}

DIFS_RES_GRID  <- c(seq(0.01, 0.1, 0.01), seq(0.1, 1.5, 0.05), seq(1.5, 3, 0.1))
DIFS_TOPS_GRID <- seq(500, 2500, 100)
DIFS_MIX_RATIO <- rbind(c(1, 1), c(2, 3), c(1, 2), c(1, 3), c(3, 2), c(2, 1), c(3, 1))

DIFS_GATE_RULES <- c("submitted", "combined", "prop_only", "n_over_4k",
                     "fixed35", "fixed10", "min_prop_n2k")

difs_gate_threshold <- function(n, k = NULL, prop = 0.05, rule = "submitted") {
  rule <- match.arg(rule, DIFS_GATE_RULES)
  need_k <- function() {
    if (is.null(k) || is.na(k) || k < 1)
      stop("gate_rule '", rule, "' needs k, and k is ", format(k))
    as.numeric(k)
  }
  switch(rule,
    submitted    = max(prop * n, 35),
    combined     = max(min(prop * n, n / (4 * need_k())), 30),
    prop_only    = prop * n,
    fixed35      = 35,
    fixed10      = 10,
    n_over_4k    = n / (4 * need_k()),
    min_prop_n2k = min(prop * n, n / (2 * need_k())))
}

difs_gate_genes <- function(data.log.mat, min.expression = log(5),
                            min.expression.proportion = 0.05,
                            gate_rule = "submitted", gate_k = NULL) {
  thresh <- difs_gate_threshold(ncol(data.log.mat), k = gate_k,
                                prop = min.expression.proportion,
                                rule = gate_rule)
  n_expr <- rowSums(data.log.mat > min.expression)
  rownames(data.log.mat)[n_expr >= thresh]
}

difs_null_path <- function() {
  p <- Sys.getenv("DIFS_NULL_TABLES", "")
  if (nzchar(p) && file.exists(p)) return(p)
  for (cand in c("../results/dip_null_tables.rds", "dip_null_tables.rds",
                 file.path(Sys.getenv("DIFS_ROOT", ".."),
                           "results/dip_null_tables.rds")))
    if (file.exists(cand)) return(cand)
  stop("dip null tables not found. Build them once with\n",
       "  singularity exec $DIFS_SIF Rscript difs_null_tables.R\n",
       "or set DIFS_NULL_TABLES to the .rds path.")
}

.difs_null_cache <- new.env(parent = emptyenv())
difs_null_tables <- function(ref) {
  if (is.null(.difs_null_cache$tabs))
    .difs_null_cache$tabs <- readRDS(difs_null_path())
  t <- .difs_null_cache$tabs[[ref]]
  if (is.null(t)) stop("null table has no reference '", ref, "'; it has: ",
                       paste(names(.difs_null_cache$tabs), collapse = ", "))
  t
}

dip_pvalue_lookup_vec <- function(D, n, tab) {
  lg <- log(tab$n_grid); upper <- 1 - tab$probs; logQ <- log(tab$Q)
  n_cl <- pmin(pmax(n, min(tab$n_grid)), max(tab$n_grid))
  out <- numeric(length(D)); cache <- new.env(parent = emptyenv())
  qcurve_for <- function(nn) {
    key <- format(nn, digits = 10)
    if (!is.null(cache[[key]])) return(cache[[key]])
    ln <- log(nn)
    i1 <- findInterval(ln, lg, rightmost.closed = TRUE)
    i1 <- min(max(i1, 1L), length(lg) - 1L); i2 <- i1 + 1L
    w <- if (lg[i2] > lg[i1]) (ln - lg[i1]) / (lg[i2] - lg[i1]) else 0
    q <- cummax(exp((1 - w) * logQ[i1, ] + w * logQ[i2, ]))
    k <- length(q); idx <- (k - 39):k
    cf <- stats::coef(stats::lm(log(upper[idx]) ~ q[idx]))
    val <- list(q = q, cf = cf); assign(key, val, envir = cache); val
  }
  for (i in seq_along(D)) {
    if (!is.finite(D[i])) { out[i] <- NA_real_; next }
    v <- qcurve_for(n_cl[i]); q <- v$q
    if (D[i] <= q[1L]) { out[i] <- 1; next }
    if (D[i] >= q[length(q)]) {
      out[i] <- max(.Machine$double.xmin, exp(v$cf[1] + v$cf[2] * D[i])); next
    }
    out[i] <- stats::approx(q, upper, xout = D[i], rule = 2)$y
  }
  out
}

DIFS_STAGE1_KEYS <- c("normal_lookup", "tnormal_lookup", "uniform_lookup",
                      "hartigan", "mc_normal", "D_raw")

difs_stage1_ranking <- function(data.log.mat, min.expression = log(5),
                                min.expression.proportion = 0.05,
                                gate_rule = "submitted", gate_k = NULL,
                                key = "normal_lookup", B = 2000, seed = 1) {
  key <- match.arg(key, DIFS_STAGE1_KEYS)
  g <- difs_gate_genes(data.log.mat, min.expression, min.expression.proportion,
                       gate_rule = gate_rule, gate_k = gate_k)
  if (!length(g)) stop("no gene passed the stage I expression gate")

  if (key == "mc_normal") {
    set.seed(seed)
    vals <- vapply(g, function(gn) {
      x <- data.log.mat[gn, ]; dip.test.v1(x[x > min.expression], B = B)
    }, numeric(1))
  } else if (key == "hartigan") {
    vals <- vapply(g, function(gn) {
      x <- data.log.mat[gn, ]
      suppressWarnings(diptest::dip.test(x[x > min.expression])$p.value)
    }, numeric(1))
  } else {
    Dn <- vapply(g, function(gn) {
      x <- data.log.mat[gn, ]; x <- x[x > min.expression]
      c(diptest::dip(sort(x)), length(x))
    }, numeric(2))
    D <- Dn[1, ]; nn <- Dn[2, ]
    vals <- if (key == "D_raw") -D else
      dip_pvalue_lookup_vec(D, nn, difs_null_tables(sub("_lookup$", "", key)))
  }
  names(vals) <- g
  names(sort(vals))
}

difs_features <- function(data.log.mat, data.seu, cluster_count,
                          binary.bound = log(5),
                          min.expression = binary.bound,   # as actually implemented
                          min.expression.proportion = 0.05,
                          gate_rule = "submitted",
                          gate_k = cluster_count,
                          tune = TRUE,
                          tune_stage1 = tune,   # the 500-2500 sweep: OPTIONAL
                          tune_ratio  = tune,   # the 7-ratio search at fixed n
                          stage1_key = "normal_lookup",
                          use_stage2 = TRUE,
                          seed = 1,
                          tops = DIFS_TOPS_GRID,
                          n_fixed = 1000,
                          stage2_alpha = 0.05,
                          min.mean.expr = 0.1) {

  ranked <- difs_stage1_ranking(
    data.log.mat              = data.log.mat,
    min.expression            = min.expression,
    min.expression.proportion = min.expression.proportion,
    gate_rule = gate_rule, gate_k = gate_k,
    key = stage1_key, seed = seed)

  stage1 <- if (tune_stage1) {
    num_features_decider(data.seu, cluster_count = cluster_count,
                         tops = tops, top_genes = ranked)$markers
  } else {
    ranked[seq_len(min(n_fixed, length(ranked)))]
  }

  if (!use_stage2) {
    only1 <- ranked[seq_len(min(n_fixed, length(ranked)))]
    return(list(markers = only1, stage1_n = length(only1),
                stage2_n = NA_integer_, ratio = NA_character_))
  }

  data.seu <- seudo_clustering(data.seu, markers = stage1,
                               cluster_count = cluster_count,
                               res = DIFS_RES_GRID, clustering.method = "SC3")
  filtered.seu <- if (any(is.na(data.seu$seurat_clusters))) {
    data.seu[, !is.na(data.seu$seurat_clusters)]
  } else data.seu

  sorted_p <- fisher_test(filtered.seu, "seurat_clusters",
                          binary.bound = binary.bound,
                          min.mean.expr = min.mean.expr)
  stage2 <- names(sorted_p[sorted_p < stage2_alpha])
  if (!length(stage2)) return(list(markers = stage1, stage1_n = length(stage1),
                                   stage2_n = 0, ratio = NA_character_))

  if (tune_ratio) {
    mixture_lists <- apply(DIFS_MIX_RATIO, 1, function(x)
      num_features_decider(filtered.seu, cluster_count = cluster_count,
                           tops = n_fixed,
                           top_genes = interleave(stage2, stage1, ratio = x)))
    best <- which.min(vapply(mixture_lists, function(z) z$mse, numeric(1)))
    return(list(markers = mixture_lists[[best]]$markers,
                stage1_n = length(stage1), stage2_n = length(stage2),
                ratio = paste(DIFS_MIX_RATIO[best, ], collapse = ":")))
  }
  mixed <- interleave(stage2, stage1, ratio = c(1, 1))
  list(markers = mixed[seq_len(min(n_fixed, length(mixed)))],
       stage1_n = length(stage1), stage2_n = length(stage2), ratio = "1:1")
}

difs_binomial_deviance <- function(cnt, chunk = 1000L) {
  nj <- colSums(cnt)
  N  <- sum(nj)
  if (!is.finite(N) || N <= 0)
    stop("binomial deviance: the count matrix sums to ", N,
         " -- this is not a counts slot")
  G   <- nrow(cnt)
  out <- numeric(G)
  for (start in seq(1L, G, by = chunk)) {
    ii <- start:min(start + chunk - 1L, G)
    y  <- cnt[ii, , drop = FALSE]
    pg <- rowSums(y) / N
    mu <- outer(pg, nj)                 # n_j * pi_g
    ny <- sweep(-y,  2L, nj, "+")       # n_j - y_gj
    nm <- sweep(-mu, 2L, nj, "+")       # n_j - mu_gj = n_j (1 - pi_g)
    t1 <- y  * log(y  / mu)
    t2 <- ny * log(ny / nm)
    t1[!is.finite(t1)] <- 0             # 0 log 0 = 0, and all-zero genes
    t2[!is.finite(t2)] <- 0
    out[ii] <- 2 * (rowSums(t1) + rowSums(t2))
  }
  out
}

feature_ranking <- function(method, data.seu, data.log.mat, data.mat, cds,
                            cluster_count, max_n = 2500,
                            monocle_decreasing = TRUE, seed = 1,
                            min.expression = log(5),
                            min.expression.proportion = 0.05,
                            gate_rule = "submitted", gate_k = NULL) {
  set.seed(seed)
  genes <- rownames(data.log.mat)
  switch(method,

    Seurat = {
      data.seu <- FindVariableFeatures(data.seu, nfeatures = max_n, verbose = FALSE)
      VariableFeatures(data.seu)
    },

    FEAST = {
      ixs <- if (ncol(data.log.mat) >= 1000)
        FEAST_fast(data.log.mat, k = cluster_count, batch_size = 1000)
      else FEAST(data.log.mat, k = cluster_count)
      genes[ixs]
    },

    SC3 = SC3_filter_genes(data.log.mat, lower_threshold = 0.06,
                           upper_threshold = 0.94),

    SC3_random1000 = {
      f <- SC3_filter_genes(data.log.mat, lower_threshold = 0.06,
                            upper_threshold = 0.94)
      sample(f, min(1000, length(f)), replace = FALSE)
    },

    Variance = genes[order(apply(data.log.mat, 1, stats::var), decreasing = TRUE)],

    Deviance = {
      if (!identical(rownames(data.mat), genes))
        stop("Deviance: counts and log matrices disagree on gene order (",
             nrow(data.mat), " vs ", length(genes), " genes)")
      genes[order(difs_binomial_deviance(data.mat), decreasing = TRUE)]
    },

    monocle = {
      disp <- dispersionTable(cds)
      disp <- disp[order(disp$dispersion_fit, decreasing = monocle_decreasing), ]
      as.character(head(disp, max_n)$gene_id)
    },

    Seurat_gated = {
      g <- difs_gate_genes(data.log.mat, min.expression, min.expression.proportion,
                           gate_rule = gate_rule, gate_k = gate_k)
      if (length(g) < 10L)
        stop("Seurat_gated: the gate passes only ", length(g),
             " genes; vst cannot be fitted")
      sub <- subset(data.seu, features = g)
      sub <- FindVariableFeatures(sub, selection.method = "vst",
                                  nfeatures = min(max_n, length(g)),
                                  verbose = FALSE)
      VariableFeatures(sub)
    },

    GateOnly = {
      g <- difs_gate_genes(data.log.mat, min.expression, min.expression.proportion,
                           gate_rule = gate_rule, gate_k = gate_k)
      sample(g, length(g), replace = FALSE)
    },

    Random = sample(genes, length(genes), replace = FALSE),

    stop("unknown feature selection method: ", method))
}

Run.method.comparison.controlled <- function(
    data.path, clustering.method, feature.selection.method,
    k_policy = c("original", "estimated", "true"),
    n_policy = c("original", "tuned_all", "fixed_stage1", "fixed_all"),
    n_fixed = 1000, tops = DIFS_TOPS_GRID,
    gate_rule = "submitted",
    min.expression = log(5), seed = 1) {

  k_policy <- match.arg(k_policy); n_policy <- match.arg(n_policy)
  gate_rule <- match.arg(gate_rule, DIFS_GATE_RULES)
  set.seed(seed)
  t0 <- proc.time()[["elapsed"]]

  data.name  <- gsub(".*/([^.]*).*", "\\1", data.path)
  data.seu   <- readRDS(file = data.path)
  trueclass  <- data.seu$trueclass
  data.log.mat <- as.matrix(GetAssayData(data.seu, slot = "data"))
  data.mat   <- as.matrix(GetAssayData(data.seu, slot = "counts"))

  sce <- SingleCellExperiment(assays = list(counts = data.mat,
                                            logcounts = log2(data.mat + 1)))
  n_unlabelled <- sum(is.na(trueclass))
  if (n_unlabelled)
    message("  ", n_unlabelled, " of ", length(trueclass),
            " cells have no label; excluded from k_true and from scoring")
  k_true <- length(unique(trueclass[!is.na(trueclass)]))

  k_estimated <- tryCatch({
    metadata(sc3_estimate_k(sce))$sc3$k_estimation
  }, error = function(e) {
    message("  sc3_estimate_k failed (", conditionMessage(e),
            "); k_estimated = NA")
    NA_integer_
  })
  if (k_policy == "estimated" && is.na(k_estimated))
    stop("k_policy = 'estimated' but sc3_estimate_k() failed on ", data.name)

  needs_cds <- identical(clustering.method, "monocle") ||
               identical(feature.selection.method, "monocle")
  cds <- if (!needs_cds) NULL else tryCatch({
    z <- suppressWarnings(sce_to_monocle(sce))
    z <- suppressWarnings(estimateSizeFactors(z))
    suppressWarnings(estimateDispersions(z))
  }, error = function(e)
    stop("monocle setup failed on ", data.name, ": ", conditionMessage(e)))

  k_fs      <- switch(k_policy, original = NA, estimated = k_estimated, true = k_true)
  k_cluster <- switch(k_policy, original = k_true, estimated = k_estimated, true = k_true)
  k_difs    <- if (k_policy == "original") k_true      else k_fs
  k_feast   <- if (k_policy == "original") k_estimated else k_fs

  difs_tune_stage1 <- switch(n_policy,
    original = TRUE, tuned_all = TRUE, fixed_stage1 = FALSE, fixed_all = FALSE)
  difs_tune_ratio  <- switch(n_policy,
    original = TRUE, tuned_all = TRUE, fixed_stage1 = TRUE,  fixed_all = FALSE)
  tune_this <- identical(n_policy, "tuned_all")

  difs_info <- list(stage1_n = NA_integer_, stage2_n = NA_integer_, ratio = NA_character_)
  k_gate <- if (is.na(k_difs)) k_true else k_difs
  gate_n <- length(difs_gate_genes(data.log.mat, min.expression,
                                   gate_rule = gate_rule, gate_k = k_gate))
  if (!gate_n)
    stop("gate_rule '", gate_rule, "' passes no gene on ", data.name,
         " (", ncol(data.log.mat), " cells, k = ", k_gate, "). ",
         "Nothing downstream is meaningful; exclude this dataset or change ",
         "the rule.")
  if (n_fixed <= 0) n_fixed <- gate_n

  DIFS_VARIANTS <- c(DIFS          = "normal_lookup",   # production
                     DIFS_mc       = "mc_normal",       # the submitted version
                     DIFS_tnorm    = "tnormal_lookup",
                     DIFS_unif     = "uniform_lookup",
                     DIFS_hartigan = "hartigan",
                     DIFS_Draw     = "D_raw",
                     DIFS_stage1only = "normal_lookup")
  use_stage2 <- !identical(feature.selection.method, "DIFS_stage1only")
  if (feature.selection.method %in% names(DIFS_VARIANTS)) {
    markers <- difs_features(
      data.log.mat = data.log.mat, data.seu = data.seu,
      cluster_count = k_difs, binary.bound = binary.bound.decider(data.log.mat),
      min.expression = min.expression,
      gate_rule = gate_rule, gate_k = k_gate,
      tune_stage1 = difs_tune_stage1, tune_ratio = difs_tune_ratio,
      stage1_key = DIFS_VARIANTS[[feature.selection.method]],
      use_stage2 = use_stage2, seed = seed,
      tops = tops, n_fixed = n_fixed)
    difs_info <- markers[c("stage1_n", "stage2_n", "ratio")]
    markers <- markers$markers
  } else {
    ranking <- feature_ranking(
      feature.selection.method, data.seu = data.seu,
      data.log.mat = data.log.mat, data.mat = data.mat, cds = cds,
      cluster_count = if (feature.selection.method == "FEAST") k_feast else k_fs,
      max_n = max(tops), seed = seed,
      min.expression = min.expression,
      gate_rule = gate_rule, gate_k = k_gate)
    markers <- if (feature.selection.method %in% c("SC3", "SC3_random1000")) {
      ranking
    } else if (tune_this) {
      num_features_decider(data.seu, cluster_count = k_cluster,
                           tops = tops, top_genes = ranking)$markers
    } else {
      ranking[seq_len(min(n_fixed, length(ranking)))]
    }
  }
  markers <- intersect(markers, rownames(data.seu))

  data.seu <- seudo_clustering(data.seu, markers, cluster_count = k_cluster,
                               clustering.method = clustering.method,
                               res = DIFS_RES_GRID, cds = cds)

  cl  <- data.seu$seurat_clusters
  ok  <- !is.na(cl) & !is.na(trueclass)
  n_dropped <- sum(!ok)
  if (n_dropped) message("  ", n_dropped, " of ", length(cl),
                         " cells had no cluster assignment; excluded from scoring")
  if (sum(ok) < 2L || length(unique(cl[ok])) < 2L) {
    perf <- rep(NA_real_, 4L)
    n_clusters_found <- length(unique(cl[ok]))
  } else {
    perf <- eval_Cluster(droplevels(factor(cl[ok])), trueclass[ok])
    n_clusters_found <- length(unique(cl[ok]))
  }

  k_hit <- isTRUE(n_clusters_found == k_cluster)

  list(data.name = data.name,
       clustering.method = clustering.method,
       feature.selection.method = feature.selection.method,
       k_policy = k_policy, n_policy = n_policy, seed = seed,
       gate_rule = gate_rule, gate_k = k_gate,
       min.expression = min.expression,
       k_true = k_true, k_estimated = k_estimated,
       k_used_feature_selection = if (feature.selection.method %in% names(DIFS_VARIANTS)) k_difs else
         if (feature.selection.method == "FEAST") k_feast else k_fs,
       k_used_clustering = k_cluster,
       n_features = length(markers),
       gate_n = gate_n,
       n_clusters_found = n_clusters_found,
       k_hit = k_hit,
       n_cells_total = length(cl),
       n_cells_dropped = n_dropped,
       n_cells_unlabelled = n_unlabelled,
       difs_use_stage2 = if (feature.selection.method %in% names(DIFS_VARIANTS))
                           use_stage2 else NA,
       difs_stage1_n = difs_info$stage1_n,
       difs_stage2_n = difs_info$stage2_n,
       difs_ratio = difs_info$ratio,
       ARI = perf[1], Purity = perf[2], Jaccard = perf[3], FM = perf[4],
       elapsed_sec = proc.time()[["elapsed"]] - t0)
}
