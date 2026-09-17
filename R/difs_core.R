################################################################################
## difs_core.R -- the DIFS stage I criterion, with no Bioconductor dependency
##
##   source("R/difs_core.R")            # needs only R >= 4.1 and diptest
##
## WHY THIS FILE EXISTS.  The reviewers asked for code that runs end to end.
## The full benchmark needs Seurat, SC3, FEAST, monocle and TSCAN, a cluster,
## and about a day of compute -- that is not something a reader can check.  The
## METHOD, though, is small and has no such dependencies: an expression gate and
## a dip-statistic ranking read off a precomputed null table.  This file is that
## method and nothing else, so it can be run on a laptop in seconds and read in
## one sitting.
##
## The benchmark harness that produced the paper's tables lives in
## benchmark/ and is documented separately; it is not needed to use DIFS.
##
## Every function here is copied verbatim from the version used for the
## manuscript (benchmark/Functions_controlled.R) so that the two cannot drift.
## If you change one, change both.
################################################################################

`%||%` <- function(a, b) if (is.null(a)) b else a

if (!requireNamespace("diptest", quietly = TRUE))
  stop("DIFS needs the 'diptest' package: install.packages(\"diptest\")")

################################################################################
## DIFS feature selection with every hidden choice made explicit
################################################################################
################################################################################
## Stage I: the gate, and the three candidate ranking keys
################################################################################
## Reproduces dip_test_pip.v1's gate exactly: a CELL counts as expressing when
## its log-normalised value exceeds min.expression (= log(5) in the pipeline),
## and a GENE is testable when at least max(0.05*n, 35) cells do.
##
## THE GATE RULE IS NOW A PARAMETER.  The submitted rule is max(0.05n, 35).
## The 21-dataset gate scan showed it is the WORST of the eight rules tried on
## the one criterion that matters structurally: it excludes 59 annotated cell
## types (across 10 datasets) whose specific markers cannot physically pass the
## threshold, because the class is smaller than the number of expressing cells
## demanded.  n/(4k) excludes 23 at essentially the same pool size (median
## gate_n 1562 vs 1569; median 7.4% vs 8.3% of genes).  Concretely: Baron has a
## 7-cell class against a gate of 428, Muraro a 3-cell class against a gate of
## 132.  That is a blind spot independent of expression level, and it lines up
## with E7's weakest perturbation (absorb_rarest_type, recall 0.74-0.83).
##
## The scan is a STRUCTURAL argument only -- it says which classes CAN be seen,
## not whether clustering improves.  gate_rule exists so the downstream question
## can be answered with the same machinery, on the same data, changing the pool
## and nothing else.  That makes it the one clean manipulation of gate_n we
## have: the budget-ratio relation gain = 0.148 - 0.156 * (n/gate_n) predicts
## that ENLARGING the pool at fixed n RAISES the dip ranking's contribution,
## and the rules below move gate_n over a 2-3x range on the same matrix.
##
## n/(4k) and min(0.05n, n/(2k)) consume k.  That is not new label leakage:
## difs_features already takes cluster_count as a user-supplied argument, and
## FEAST takes k as well.  Whatever k_policy hands the method is what the gate
## sees -- no rule may reach past it to the labels.
##
## "combined" IS THE RULE WE PROPOSE, and it is built from two separate
## arguments rather than from a scan of what performs best -- performance may
## only VETO it (see the pre-registration), never select it.
##
##   gate = max( min(0.05n, n/(4k)), 30 )
##
##   min(0.05n, n/(4k))  -- at the large-n end the submitted 0.05n grows without
##       bound (Baron 428, Shekhar 2246) and that is where the structural damage
##       is: 59 annotated cell types across 21 datasets cannot contribute a
##       marker.  n/(4k) scales with the expected size of ONE cluster instead of
##       with the sample, so it takes over exactly there (Baron 428 -> 153,
##       Shekhar 2246 -> 591) and roughly doubles the candidate pool on most
##       datasets (+33% to +128%).
##   floor 30 -- justified by the POWER analysis, not by observation.  E3: power
##       is ~0.03 when the differing subpopulation has <= 25 cells, and E2 shows
##       the modified test is already miscalibrated at low expression (type-I
##       error 0.09-0.22 under NB / ZINB, from count discreteness).  Below ~30
##       points stage I is both powerless and invalid, so the floor is where the
##       statistic stops working, not a round number.  The submitted 35 is inside
##       this range; the change from 35 to 30 is cosmetic and is NOT the point.
##
## WHAT THIS DOES NOT FIX, and must be stated as such.  Checked per dataset: the
## SMALLEST annotated class still fails the gate on all eight datasets where it
## failed before -- Muraro's 3-cell class against a new threshold of 66, Baron's
## 7-cell class against 153, Koh's 22-cell class against 30.  Any threshold that
## leaves the dip test the ~30 observations it needs is necessarily larger than
## these classes.  The rare-type blind spot is therefore a structural property of
## testing a distribution for bimodality, NOT a badly chosen constant, and no
## amount of threshold tuning removes it.  The rule reduces the number of
## excluded classes; it cannot reach zero.
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

## Three ways to order the gate-passing genes.  The ranking_check on Koh found:
##   mc_normal  874 distinct values out of 1156  -> 282 genes tie, broken by
##              matrix row order; and two runs of it disagree with each other
##              (chance-corrected top-500 agreement 0.97, not 1)
##   hartigan  1156 distinct values out of 1156  -> NO ties, deterministic
##   D_raw     1156 distinct, but not adjusted for the number of expressing
##             cells (median n in the top 1000: 79.5 vs 96.0 for the p-values)
## The manuscript motivates the modified test by saying the original cannot
## rank because most p-values are 1.  Among genes that actually enter the test
## that is false, so this is the ablation that settles which key to keep.
################################################################################
## Precomputed null tables: the same normal calibration, computed accurately
################################################################################
## The normal reference is a deliberate modelling choice -- DIFS assumes
## non-informative genes are unimodal and approximately normal, so it calibrates
## against that rather than against Hartigan's least-favourable uniform.  That
## choice is defensible.  Estimating the reference by B = 2000 draws PER GENE is
## what is not:
##   * p-values land on a 1/B grid and tie.  Koh: 1156 genes, 874 distinct
##     values, and 3 genes at exactly 0 -- at the selection end, broken by row
##     order.  Hartigan's tabulated p gives 1156 distinct values and 0 ties
##     there, so the manuscript's "the original test cannot rank" is backwards.
##   * two runs disagree with each other (chance-corrected top-500 agreement
##     0.97).  A feature selection method should be deterministic.
##   * it is the runtime Reviewer 1 objected to.
## Reading the p-value off a precomputed table of null quantiles removes all
## three and changes nothing about the statistical assumption.
difs_null_path <- function() {
  p <- Sys.getenv("DIFS_NULL_TABLES", "")
  if (nzchar(p) && file.exists(p)) return(p)
  ## In this repository the table ships alongside the code, so a fresh clone
  ## works with no build step.  The other candidates are the layout used on the
  ## cluster, kept so that one copy of this file serves both.
  for (cand in c("inst/extdata/dip_null_tables.rds",
                 "../inst/extdata/dip_null_tables.rds",
                 file.path(dirname(sys.frames()[[1]]$ofile %||% "."),
                           "..", "inst", "extdata", "dip_null_tables.rds"),
                 "../results/dip_null_tables.rds", "dip_null_tables.rds",
                 file.path(Sys.getenv("DIFS_ROOT", ".."),
                           "results/dip_null_tables.rds")))
    if (!is.null(cand) && !is.na(cand) && file.exists(cand)) return(cand)
  stop("dip null tables not found.  In a clone of this repository they are at\n",
       "  inst/extdata/dip_null_tables.rds\n",
       "and are found automatically when R's working directory is the repository\n",
       "root.  Otherwise set DIFS_NULL_TABLES to the file's path, or rebuild the\n",
       "table from scratch with benchmark/difs_null_tables.R.")
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

## Upper-tail p from the table.  Interpolates across n on the log scale (the
## null dip scales roughly as n^-1/2) and fits an exponential to the far upper
## tail, so the most extreme genes keep distinct positive p-values instead of
## all collapsing onto zero -- which is exactly the failure the Monte Carlo
## version has.
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

## Ranking keys for stage I.  All are sorted ASCENDING, so smaller = stronger
## evidence of multimodality, for every key.
##   normal_lookup   the recommended production key: normal reference, read off
##                   the precomputed table.  Deterministic, tie-free, ~2000x
##                   faster than the per-gene Monte Carlo.
##   tnormal_lookup  normal LEFT-TRUNCATED at the gate.  Stage I tests
##                   x[x > log(5)], so the sample it sees is truncated; this is
##                   the reference that actually matches it.
##   uniform_lookup  Hartigan's reference via the same machinery, so the
##                   comparison isolates the REFERENCE and not the estimator.
##   hartigan        diptest::dip.test(), the published implementation.
##   mc_normal       the submitted per-gene Monte Carlo.  Kept only to
##                   reproduce the submitted numbers; do not use for new work.
##   D_raw           the dip statistic with no adjustment for n at all.
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
    ## one pass for the statistic and the sample size, then a vectorised lookup
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

################################################################################
## Baseline feature rankings -- each returns genes ordered best-first
################################################################################
## monocle: NOTE that monocle_gene_filter() in Functions_mixture_auto.R orders
## dispersionTable() by dispersion_fit with decreasing = FALSE, i.e. it keeps the
## genes with the LOWEST dispersion, which is the opposite of what a
## dispersion-based feature selection should do.  The version used here orders
## by decreasing dispersion; set monocle_decreasing = FALSE to reproduce the
## original behaviour.
## BINOMIAL DEVIANCE, Townes et al., Genome Biology 2019 (the criterion behind
## scry::devianceFeatureSelection).  Added 2026-09-15.
##
## WHY IT IS HERE.  Two Genome Biology papers independently recommend deviance
## as the feature-selection criterion of choice -- Townes et al. 2019, and
## pipeComp (Germain et al. 2020), whose finding is
##   "deviance and unstandardized estimates of variance ... provide the best
##    results across datasets and normalization methods"
## and
##   "standardized measures of variability were systematically worse than their
##    non-standardized counterparts".
## Our own grid reproduces the second half of that (Seurat vst loses to a plain
## variance ranking by 0.071 mean ARI, p = 3e-4), which makes the ABSENCE of
## deviance the first thing a reviewer who knows this literature will ask about.
## It is also the honest test of the method: if DIFS loses to deviance, we need
## to know before the claim is written, not after it is published.
##
## For gene g with counts y_gj in cell j and library size n_j:
##   pi_g = sum_j y_gj / sum_j n_j        (the gene's share of all counts)
##   mu_gj = n_j * pi_g                   (expected counts under the null)
##   D_g  = 2 * sum_j [ y_gj log(y_gj / mu_gj)
##                    + (n_j - y_gj) log((n_j - y_gj) / (n_j - mu_gj)) ]
## with the convention 0 log 0 = 0, which is what the !is.finite() sweeps do.
## Genes are ranked by DECREASING deviance: a large deviance means the gene's
## counts are poorly explained by a constant rate across cells.
##
## Chunked over genes on purpose.  The runner hands this a DENSE matrix, and on
## Baron (16359 x 8569) a single un-chunked temporary is 1.1 GB -- five of them
## at once is how a 64G job dies without an error message in the .Rout.
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
