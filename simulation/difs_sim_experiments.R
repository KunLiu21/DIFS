

E1_invariance <- function(n = c(50, 200, 1000), reps = 200, seed = 1,
                          tab_norm = NULL,
                          out_dir = "../results/simulation") {
  set.seed(seed)
  res <- do.call(rbind, lapply(n, function(nn) {
    do.call(rbind, lapply(seq_len(reps), function(r) {
      x <- stats::rnorm(nn)
      a <- stats::runif(1, -10, 10); b <- stats::runif(1, 0.1, 10)
      data.frame(n = nn, rep = r,
                 dip_x    = dip_stat(x),
                 dip_affine = dip_stat(a + b * x),
                 dip_neg  = dip_stat(-b * x))
    }))
  }))
  res$max_abs_diff <- pmax(abs(res$dip_x - res$dip_affine),
                           abs(res$dip_x - res$dip_neg))
  save_rds(res, "E1_invariance", out_dir)

  mono <- do.call(rbind, lapply(n, function(nn) {
    x  <- lapply(seq_len(200), function(i)
      if (i %% 2) stats::rnorm(nn) else c(stats::rnorm(nn %/% 2),
                                          stats::rnorm(nn - nn %/% 2, 4))
    )
    D  <- vapply(x, dip_stat, numeric(1))
    pu <- vapply(x, function(v) unname(diptest::dip.test(v)$p.value), numeric(1))
    nu <- dip_null_sample(nn, 2000, "normal")
    pn_mc <- vapply(D, function(d) dip_pvalue_mc(d, nu), numeric(1))
    pn_lk <- if (!is.null(tab_norm)) dip_pvalue_lookup_vec(D, rep(nn, length(D)),
                                                           tab_norm) else NA_real_
    disc <- sum(outer(D, D, ">") & outer(pn_mc, pn_mc, ">"))
    data.frame(n = nn,
               spearman_D_vs_pnormal_mc     = stats::cor(D, pn_mc, method = "spearman"),
               spearman_D_vs_pnormal_lookup = if (length(pn_lk) > 1)
                 stats::cor(D, pn_lk, method = "spearman") else NA_real_,
               spearman_puniform_vs_pnormal = stats::cor(pu, pn_mc, method = "spearman"),
               n_strictly_discordant_pairs  = disc,
               n_tied_pvalues_mc            = sum(duplicated(pn_mc)))
  }))
  save_rds(mono, "E1_monotonicity", out_dir)
  list(invariance = res, monotonicity = mono)
}

E2_null_calibration <- function(n_cells = c(200, 500, 1000),
                                types = c("normal", "lognormal", "gamma",
                                          "nb", "zinb", "zinb_high"),
                                mus = c(1, 3, 10, 30, 100),
                                n_genes = 1000, disp = 0.5,
                                libsize_sd = 0.35,
                                alphas = c(0.01, 0.05, 0.10),
                                min.expression = DIFS_GATE_AS_DOCUMENTED,
                                norm_method = "lognorm",
                                tab_norm, tab_unif, seed = 2,
                                out_dir = "../results/simulation") {
  set.seed(seed)
  out <- list()
  for (nc in n_cells) {
    sf <- stats::rlnorm(nc, 0, libsize_sd)
    for (mu in mus) for (ty in types) {
      cnt <- t(vapply(seq_len(n_genes), function(g)
        sim_null_gene(ty, nc, mu = mu, disp = disp, size_factor = sf), numeric(nc)))
      if (ty %in% CONTINUOUS_NULL_TYPES) {
        lm_ <- cnt
      } else {
        bg  <- matrix(stats::rnbinom(2000 * nc, mu = 5, size = 2), nrow = 2000)
        lm_ <- normalize_counts(rbind(cnt, bg), norm_method)[seq_len(n_genes), ,
                                                             drop = FALSE]
      }
      st  <- difs_stage1_matrix(lm_, method = "lookup", tab = tab_norm,
                                min.expression = min.expression)
      st_u <- difs_stage1_matrix(lm_, method = "hartigan",
                                 min.expression = min.expression)
      tested <- st$p <= 1
      if (sum(tested) < 10) {                      # too few testable genes
        message(sprintf("  E2 n=%4d mu=%5.0f %-10s tested=%4d  (skipped)",
                        nc, mu, ty, sum(tested)))
        next
      }
      for (a in alphas) {
        out[[length(out) + 1]] <- data.frame(
          n_cells = nc, null_type = ty, mu = mu, alpha = a,
          min_expression = min.expression, norm_method = norm_method,
          n_tested = sum(tested),
          frac_tested = mean(tested),
          median_n_expressing = stats::median(st$n[tested]),
          median_distinct_frac = stats::median(
            vapply(which(tested), function(i) {
              xe <- lm_[i, ][lm_[i, ] > min.expression]
              length(unique(xe)) / max(length(xe), 1)
            }, numeric(1))),
          fpr_modified = mean(st$p[tested] < a),
          fpr_original = mean(st_u$p[tested] < a))
      }
      ks_m <- suppressWarnings(stats::ks.test(st$p[tested], "punif")$statistic)
      ks_o <- suppressWarnings(stats::ks.test(st_u$p[tested], "punif")$statistic)
      for (i in seq(length(out) - length(alphas) + 1, length(out))) {
        out[[i]]$ks_modified <- unname(ks_m); out[[i]]$ks_original <- unname(ks_o)
      }
      message(sprintf("  E2 n=%4d mu=%5.0f %-10s tested=%4d  FPR.05 mod=%.3f orig=%.3f",
                      nc, mu, ty, sum(tested), mean(st$p[tested] < .05),
                      mean(st_u$p[tested] < .05)))
    }
  }
  res <- do.call(rbind, out)
  save_rds(res, "E2_null_calibration", out_dir)
  res
}

E3_power <- function(n_cells = c(200, 500, 1000),
                     pis = c(0.05, 0.10, 0.20, 0.35, 0.50),
                     log2fcs = c(1, 2, 3, 4),
                     mus = c(3, 10, 30),
                     gene_type = "zinb", n_genes = 400, disp = 0.5,
                     libsize_sd = 0.35, alpha = 0.05,
                     min.expression = DIFS_GATE_AS_DOCUMENTED,
                     norm_method = "lognorm",
                     tab_norm, seed = 3,
                     out_dir = "../results/simulation") {
  set.seed(seed)
  out <- list()
  for (nc in n_cells) for (mu in mus) {
    sf <- stats::rlnorm(nc, 0, libsize_sd)
    bg <- matrix(stats::rnbinom(2000 * nc, mu = 5, size = 2), nrow = 2000)
    cnt0 <- t(vapply(seq_len(n_genes), function(g)
      sim_null_gene(gene_type, nc, mu = mu, disp = disp, size_factor = sf),
      numeric(nc)))
    lm0 <- normalize_counts(rbind(cnt0, bg), norm_method)[seq_len(n_genes), ,
                                                          drop = FALSE]
    p0m <- difs_stage1_matrix(lm0, method = "lookup", tab = tab_norm,
                              min.expression = min.expression)
    p0u <- difs_stage1_matrix(lm0, method = "hartigan",
                              min.expression = min.expression)
    tm <- p0m$p <= 1; tu <- p0u$p <= 1
    crit_m <- if (sum(tm)) stats::quantile(p0m$p[tm], alpha, names = FALSE) else 0
    crit_u <- if (sum(tu)) stats::quantile(p0u$p[tu], alpha, names = FALSE) else 0
    for (pi in pis) for (fc in log2fcs) {
      cnt <- t(vapply(seq_len(n_genes), function(g)
        sim_alt_gene(gene_type, nc, pi = pi, log2fc = fc, mu = mu,
                     disp = disp, size_factor = sf), numeric(nc)))
      lm_ <- normalize_counts(rbind(cnt, bg), norm_method)[seq_len(n_genes), ,
                                                           drop = FALSE]
      pm <- difs_stage1_matrix(lm_, method = "lookup", tab = tab_norm,
                               min.expression = min.expression)
      pu <- difs_stage1_matrix(lm_, method = "hartigan",
                               min.expression = min.expression)
      km <- pm$p <= 1; ku <- pu$p <= 1
      out[[length(out) + 1]] <- data.frame(
        n_cells = nc, mu = mu, pi = pi, log2fc = fc, alpha = alpha,
        min_expression = min.expression,
        frac_testable = mean(km),
        power_nominal_modified = if (sum(km)) mean(pm$p[km] < alpha) else NA,
        power_nominal_original = if (sum(ku)) mean(pu$p[ku] < alpha) else NA,
        power_matched_modified = if (sum(km)) mean(pm$p[km] < crit_m) else NA,
        power_matched_original = if (sum(ku)) mean(pu$p[ku] < crit_u) else NA,
        crit_modified = crit_m, crit_original = crit_u)
      message(sprintf("  E3 n=%4d mu=%3.0f pi=%.2f fc=%d  matched power mod=%.3f orig=%.3f",
                      nc, mu, pi, fc,
                      out[[length(out)]]$power_matched_modified,
                      out[[length(out)]]$power_matched_original))
    }
  }
  res <- do.call(rbind, out)
  save_rds(res, "E3_power", out_dir)
  res
}

E4_ranking <- function(n_cells = c(300, 600, 1200), n_genes = 3000,
                       prop_informative = 0.08,
                       log2fcs = c(2, 3, 4, 5, 6, 8),
                       disps = c(0.2, 0.5),
                       mu_meanlog = log(10),
                       type_props = c(.45, .25, .15, .10, .05),
                       reps = 5, Ns = c(500, 1000, 2000), n_bg_genes = NULL,
                       min.expression = DIFS_GATE_AS_DOCUMENTED,
                       norm_method = "lognorm",
                       tab_norm, seed = 4,
                       out_dir = "../results/simulation") {
  out <- list()
  for (nc in n_cells) for (fc in log2fcs) for (dp in disps) for (r in seq_len(reps)) {
    sim <- sim_panel(n_cells = nc, n_genes = n_genes, type_props = type_props,
                     prop_informative = prop_informative, log2fc = fc,
                     disp = dp, mu_meanlog = mu_meanlog, n_bg_genes = n_bg_genes,
                     seed = seed * 1000 + nc * 97 + fc * 13 + round(dp * 10) + r)
    lm_ <- normalize_counts(sim$counts, norm_method)
    truth <- sim$informative
    st_m <- difs_stage1_matrix(lm_, method = "lookup", tab = tab_norm,
                               min.expression = min.expression)
    st_o <- difs_stage1_matrix(lm_, method = "hartigan",
                               min.expression = min.expression)
    scores <- list(
      DIFS_stage1_modified = list(s = st_m$p,                  dec = FALSE),
      dip_original         = list(s = st_o$p,                  dec = FALSE),
      dip_statistic_raw    = list(s = st_m$D,                  dec = TRUE),
      n_expressing         = list(s = st_m$n,                  dec = TRUE),
      Seurat_vst           = list(s = rank_seurat_vst(sim$counts), dec = TRUE),
      variance_lognorm     = list(s = rank_variance(lm_),      dec = TRUE),
      dropout_rate         = list(s = rank_dropoutrate(sim$counts), dec = TRUE),
      random               = list(s = stats::runif(nrow(lm_)), dec = TRUE),
      oracle_F             = list(s = rank_oracle_F(lm_, sim$labels), dec = TRUE))
    for (nm in names(scores)) {
      s <- scores[[nm]]$s; dec <- scores[[nm]]$dec
      if (!dec) s[!is.finite(s)] <- Inf else s[!is.finite(s)] <- -Inf
      ord <- order_by(s, dec)
      out[[length(out) + 1]] <- data.frame(
        n_cells = nc, log2fc = fc, disp = dp, rep = r,
        min_expression = min.expression, method = nm,
        auc = auc_score(if (dec) s else -s, truth),
        precision_500  = precision_at_n(ord, truth, min(500, n_genes)),
        precision_1000 = precision_at_n(ord, truth, min(1000, n_genes)),
        precision_2000 = precision_at_n(ord, truth, min(2000, n_genes)),
        enrichment_1000 = enrichment_at_n(ord, truth, min(1000, n_genes)),
        n_for_recall50 = genes_needed_for_recall(ord, truth, 0.5),
        n_for_recall80 = genes_needed_for_recall(ord, truth, 0.8))
    }
    message(sprintf("  E4 n=%4d log2fc=%d disp=%.1f rep=%d done", nc, fc, dp, r))
  }
  res <- do.call(rbind, out)
  save_rds(res, "E4_ranking", out_dir)
  res
}

E5_ties_runtime <- function(n_cells = 600, n_genes = 1500, B = 2000,
                            topN = c(500, 1000), tab_norm, seed = 5,
                            mu_meanlog = log(10), n_bg_genes = 0,
                            min.expression = DIFS_GATE_AS_DOCUMENTED,
                            norm_method = "lognorm",
                            out_dir = "../results/simulation") {
  sim <- sim_panel(n_cells = n_cells, n_genes = n_genes, seed = seed,
                   mu_meanlog = mu_meanlog, n_bg_genes = n_bg_genes)
  lm_ <- normalize_counts(sim$counts, norm_method)
  set.seed(seed + 1)
  t1 <- system.time(mc1 <- difs_stage1_matrix(lm_, method = "mc", ref = "normal",
                                              B = B, add_one = FALSE,
                                              min.expression = min.expression))[["elapsed"]]
  set.seed(seed + 2)
  t2 <- system.time(mc2 <- difs_stage1_matrix(lm_, method = "mc", ref = "normal",
                                              B = B, add_one = FALSE,
                                              min.expression = min.expression))[["elapsed"]]
  t3 <- system.time(lk  <- difs_stage1_matrix(lm_, method = "lookup",
                                              tab = tab_norm,
                                              min.expression = min.expression))[["elapsed"]]
  tested <- mc1$p <= 1
  jac <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  rows <- lapply(topN, function(N) {
    o1 <- rownames(mc1)[order(mc1$p)][seq_len(N)]
    o2 <- rownames(mc2)[order(mc2$p)][seq_len(N)]
    l1 <- rownames(lk)[order(lk$p)][seq_len(N)]
    data.frame(topN = N,
               jaccard_mc_run1_vs_run2 = jac(o1, o2),
               jaccard_mc_vs_lookup    = jac(o1, l1))
  })
  res <- do.call(rbind, rows)
  res$n_genes            <- n_genes
  res$n_cells            <- n_cells
  res$B                  <- B
  res$n_tested           <- sum(tested)
  res$n_tied_at_zero_mc  <- sum(mc1$p[tested] == 0)
  res$frac_tied_at_zero  <- mean(mc1$p[tested] == 0)
  res$n_distinct_p_mc    <- length(unique(mc1$p[tested]))
  res$n_distinct_p_lookup<- length(unique(lk$p[tested]))
  res$secs_mc_run        <- mean(c(t1, t2))
  res$secs_lookup        <- t3
  res$speedup            <- mean(c(t1, t2)) / max(t3, 1e-9)
  save_rds(res, "E5_ties_runtime", out_dir)
  res
}

E6_gate_reference <- function(n_cells = 600, n_genes = 2000,
                              gates = c(0, 0.25, 0.5, 1.0, DIFS_GATE_AS_IMPLEMENTED,
                                        2.0),
                              tab_norm, tab_tnorm, topN = 1000, seed = 6,
                              mu_meanlog = log(10), n_bg_genes = NULL,
                              norm_method = "lognorm",
                              out_dir = "../results/simulation") {
  sim <- sim_panel(n_cells = n_cells, n_genes = n_genes, seed = seed,
                   mu_meanlog = mu_meanlog, n_bg_genes = n_bg_genes)
  lm_ <- normalize_counts(sim$counts, norm_method); truth <- sim$informative
  ref_set <- NULL; out <- list()
  for (g in gates) for (rf in c("normal", "tnormal")) {
    tb <- if (rf == "normal") tab_norm else tab_tnorm
    st <- difs_stage1_matrix(lm_, method = "lookup", tab = tb, min.expression = g)
    s <- st$p; s[!is.finite(s)] <- Inf
    s[s > 1] <- Inf
    ord <- order_by(s, FALSE)
    n_ok <- sum(is.finite(s))
    top <- rownames(st)[ord][seq_len(min(topN, n_ok))]
    if (is.null(ref_set) && isTRUE(all.equal(g, DIFS_GATE_AS_IMPLEMENTED)) &&
        rf == "normal") ref_set <- top
    out[[length(out) + 1]] <- data.frame(
      min_expression = g, reference = rf, norm_method = norm_method,
      n_testable = n_ok,
      auc = auc_score(-s, truth),
      precision_at_topN = precision_at_n(ord, truth, min(topN, n_ok)),
      n_for_recall50 = genes_needed_for_recall(ord, truth, 0.5),
      top_set = I(list(top)))
  }
  res <- do.call(rbind, out)
  res$jaccard_vs_as_implemented <- vapply(res$top_set, function(s)
    length(intersect(s, ref_set)) / length(union(s, ref_set)), numeric(1))
  res$top_set <- NULL
  save_rds(res, "E6_gate_reference", out_dir)
  res
}

E7_stage2_robustness <- function(n_cells = 800, n_genes = 2000,
                                 type_props = c(.40, .25, .17, .12, .06),
                                 perturb = c(0, 0.05, 0.10, 0.20, 0.40),
                                 reps = 3, alpha = 0.05, seed = 7,
                                 mu_meanlog = log(10), n_bg_genes = NULL,
                                 norm_method = "lognorm",
                                 out_dir = "../results/simulation") {
  out <- list()
  for (r in seq_len(reps)) {
    sim <- sim_panel(n_cells = n_cells, n_genes = n_genes,
                     type_props = type_props, seed = seed * 100 + r,
                     mu_meanlog = mu_meanlog, n_bg_genes = n_bg_genes)
    lm_ <- normalize_counts(sim$counts, norm_method); truth <- sim$informative
    truth_lab <- sim$labels
    scen <- list()
    for (p in perturb) {
      lab <- truth_lab
      if (p > 0) {
        idx <- sample(length(lab), round(p * length(lab)))
        lab[idx] <- sample(levels(lab), length(idx), TRUE)
      }
      scen[[paste0("perturb_", p)]] <- lab
    }
    lab <- truth_lab; lab[lab == levels(lab)[2]] <- levels(lab)[1]
    scen[["merge_two_types"]] <- droplevels(lab)
    rare <- levels(truth_lab)[which.min(table(truth_lab))]
    lab <- truth_lab; lab[lab == rare] <- levels(lab)[1]
    scen[["absorb_rarest_type"]] <- droplevels(lab)
    for (nm in names(scen)) {
      s2 <- difs_stage2(lm_, scen[[nm]])
      sel <- s2$p < alpha
      out[[length(out) + 1]] <- data.frame(
        rep = r, scenario = nm,
        n_clusters = nlevels(droplevels(as.factor(scen[[nm]]))),
        ari_vs_truth = adjusted_rand(scen[[nm]], truth_lab),
        n_selected = sum(sel),
        precision = if (sum(sel)) mean(truth[sel]) else NA_real_,
        recall = if (sum(truth)) sum(truth & sel) / sum(truth) else NA_real_,
        enrichment = if (sum(sel)) mean(truth[sel]) / mean(truth) else NA_real_)
    }
    message(sprintf("  E7 rep %d done", r))
  }
  res <- do.call(rbind, out)
  save_rds(res, "E7_stage2_robustness", out_dir)
  res
}

save_rds <- function(x, name, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file.path(out_dir, paste0(name, ".rds")))
  utils::write.csv(x, file.path(out_dir, paste0(name, ".csv")), row.names = FALSE)
  invisible(x)
}

adjusted_rand <- function(a, b) {
  tb <- table(a, b); n <- sum(tb)
  ci <- function(x) sum(choose(x, 2))
  s  <- ci(as.vector(tb)); sa <- ci(rowSums(tb)); sb <- ci(colSums(tb))
  exp_ <- sa * sb / choose(n, 2); mx <- (sa + sb) / 2
  (s - exp_) / (mx - exp_)
}
