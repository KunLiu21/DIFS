# difs_sim_E8_assumption.R: revision-2 analysis source.
# Source: code/as_run/08_simulation/difs_sim_E8_assumption.R
# Usage, scope and limitations: benchmark/REVISION2.md.
lib <- Sys.getenv("DIFS_LIB", "")
if (nzchar(lib) && dir.exists(lib)) .libPaths(c(lib, .libPaths()))
if (!requireNamespace("diptest", quietly = TRUE)) stop("diptest not found on .libPaths()")
here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), error = function(e) ".")
source(file.path(here, "difs_sim_core.R"))
out_dir <- Sys.getenv("E8_OUT", "../results/simulation")
tf <- c(file.path(out_dir, "dip_null_tables.rds"), "../results/dip_null_tables.rds",
        Sys.getenv("DIFS_TABLES", ""), "../inst/extdata/dip_null_tables.rds")
tf <- tf[nzchar(tf) & file.exists(tf)][1]
if (is.na(tf)) stop("dip_null_tables.rds not found")
TABS <- readRDS(tf)
cores <- as.integer(Sys.getenv("E8_CORES", "2"))

PROPS   <- c(.45, .25, .15, .10, .05)     # as E4
N_GENES <- 3000; P_MARK <- 0.08; P_HV <- 0.08   # 240 markers, 240 high-variance nulls
SIGMA   <- c(0.4, 0.6)                    # per-gene component sd, uniform on this range
MU      <- c(2.5, 4.0)                    # per-gene log level of the expressed cells
DETECT  <- c(0.3, 0.9)                    # per-gene detection rate, same for every class
SF_SD   <- 0.3                            # cell size factors, log-normal
G <- expand.grid(r = 1:5, effect = c(2, 4, 6, 8), nc = c(300, 600, 1200))
gate <- DIFS_GATE_AS_IMPLEMENTED          # ln 5, as the pipeline uses

sim_E8 <- function(nc, effect, seed) {
  set.seed(seed)
  type <- sample(seq_along(PROPS), nc, replace = TRUE, prob = PROPS)
  pis  <- tabulate(type, length(PROPS)) / nc
  n_m <- round(P_MARK * N_GENES); n_h <- round(P_HV * N_GENES)
  cls <- c(rep("marker", n_m), rep("hv_null", n_h), rep("null", N_GENES - n_m - n_h))
  sg  <- stats::runif(N_GENES, SIGMA[1], SIGMA[2])
  mu  <- stats::runif(N_GENES, MU[1], MU[2])
  det <- stats::runif(N_GENES, DETECT[1], DETECT[2])
  hi  <- sample(seq_along(PROPS), N_GENES, replace = TRUE)   # the up type (markers) /
  delta <- effect * sg                                       # the matched type (hv nulls)
  sd_h <- sqrt(sg^2 + pis[hi] * (1 - pis[hi]) * delta^2)
  X <- matrix(0, N_GENES, nc)
  for (g in seq_len(N_GENES)) {
    x <- switch(cls[g],
      marker  = mu[g] + delta[g] * (type == hi[g]) + stats::rnorm(nc, 0, sg[g]),
      hv_null = mu[g] + delta[g] * pis[hi[g]] + stats::rnorm(nc, 0, sd_h[g]),
      null    = mu[g] + delta[g] * pis[hi[g]] + stats::rnorm(nc, 0, sg[g]))
    X[g, ] <- ifelse(stats::runif(nc) < det[g], x, -Inf)     # -Inf = not detected
  }
  sf <- stats::rlnorm(nc, 0, SF_SD)
  lam <- sweep(expm1(pmax(X, 0)), 2, sf, "*"); lam[!is.finite(X)] <- 0
  counts <- matrix(stats::rpois(length(lam), lam), N_GENES, nc)
  rownames(counts) <- paste0("g", seq_len(N_GENES))
  list(counts = counts, cls = cls, type = type)
}

one <- function(i) {
  nc <- G$nc[i]; ef <- G$effect[i]; r <- G$r[i]
  sim <- sim_E8(nc, ef, seed = 8000000 + nc * 97 + ef * 13 + r)
  cnt <- sim$counts
  lm_ <- normalize_counts(cnt, "lognorm", scale.factor = stats::median(colSums(cnt)))
  truth <- sim$cls == "marker"; hv <- sim$cls == "hv_null"
  st <- difs_stage1_matrix(lm_, method = "lookup", tab = TABS$normal, min.expression = gate)
  pass <- is.finite(st$p) & st$p <= 1
  sc <- list(DIFS_stage1 = -st$p, dip_statistic_raw = st$D, n_expressing = st$n,
             Seurat_vst = rank_seurat_vst(cnt), random = stats::runif(nrow(lm_)))
  do.call(rbind, lapply(names(sc), function(m) {
    s <- sc[[m]]
    s_all <- s; if (m == "DIFS_stage1") s_all[!pass] <- -1.1   # gate failures ranked last
    s_all[!is.finite(s_all)] <- -1e300
    o <- order(s_all, decreasing = TRUE); top <- utils::head(o, 500)
    vs_hv <- truth | hv
    data.frame(n_cells = nc, effect = ef, rep = r, method = m,
               auc_all = auc_score(s_all, truth),
               auc_within_gate = auc_score(s[pass], truth[pass]),
               auc_vs_hvnull = auc_score(s_all[vs_hv], truth[vs_hv]),
               genes_for_half = genes_needed_for_recall(o, truth, 0.5),
               top500_sens = sum(truth[top]) / sum(truth),
               top500_spec = 1 - sum(!truth[top]) / sum(!truth),
               top500_hv = sum(hv[top]), top500_markers = sum(truth[top]),
               markers_in_gate = sum(truth & pass), hv_in_gate = sum(hv & pass),
               genes_in_gate = sum(pass), markers = sum(truth))
  }))
}

t0 <- Sys.time()
res <- do.call(rbind, parallel::mclapply(seq_len(nrow(G)), one, mc.cores = cores))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(res, file.path(out_dir, "E8_assumption.csv"), row.names = FALSE)
cat("done in", format(Sys.time() - t0), "\n")
agg <- function(v) round(tapply(res[[v]], list(res$method, res$effect), mean), 3)
for (v in c("auc_all", "auc_within_gate", "auc_vs_hvnull", "genes_for_half",
            "top500_sens", "top500_spec", "top500_hv"))
  { cat("\n", v, " by effect (delta / sigma):\n", sep = ""); print(agg(v)) }
cat("\nmarkers passing the gate (share):\n")
print(round(tapply(res$markers_in_gate / res$markers, res$effect, mean), 3))
