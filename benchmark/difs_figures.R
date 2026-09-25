
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
args <- commandArgs(TRUE); out <- "../results/figures"; n_focus <- 300
grid <- "../results/diag/grid_corrected.csv"
abl  <- "../results/diag/abl_corrected.csv"
allow_defective <- 0
difs_args(args, c("out", "n_focus", "grid", "abl", "allow_defective"))
NOT_INDEPENDENT <- c("SimKumar4easy", "SimKumar4hard", "Zhengmix4eq", "Zhengmix8eq")
stamp <- function(f) cat(sprintf("  source: %s  (%s, %d bytes)\n", f,
  format(file.info(f)$mtime, "%Y-%m-%d %H:%M"), file.info(f)$size))

suppressPackageStartupMessages({ library(ggplot2) })
dir.create(out, recursive = TRUE, showWarnings = FALSE)

PAL <- c(DIFS = "#0072B2", FEAST = "#D55E00", Seurat = "#009E73",
         GateOnly = "#CC79A7", DIFS_stage1only = "#56B4E9",
         Variance = "#E69F00")
LTY <- c(DIFS = "solid", FEAST = "longdash", Seurat = "dotdash",
         GateOnly = "dotted", DIFS_stage1only = "twodash", Variance = "dashed")
MAIN <- c("DIFS", "FEAST", "Seurat")

theme_difs <- function() {
  theme_bw(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
          panel.border = element_rect(colour = "grey60", linewidth = 0.4),
          strip.background = element_rect(fill = "grey96", colour = "grey60"),
          strip.text = element_text(face = "bold", size = 8),
          legend.key.width = unit(1.4, "lines"),
          legend.position = "bottom", legend.title = element_blank())
}

long <- function(csv, pol) {
  d <- utils::read.csv(csv, stringsAsFactors = FALSE)
  m <- intersect(names(PAL), names(d))
  z <- do.call(rbind, lapply(m, function(k) data.frame(
    data.name = d$data.name, n = as.integer(d$n_requested),
    method = k, ARI = suppressWarnings(as.numeric(d[[k]])),
    stringsAsFactors = FALSE)))
  z$k_policy <- pol
  z[!is.na(z$ARI), ]
}

if (nzchar(grid) && file.exists(grid)) {
  cat("Figures A and B read the CORRECTED grid.\n"); stamp(grid)
  g <- utils::read.csv(grid, stringsAsFactors = FALSE)
  g <- g[is.finite(g$ARI), ]
  z <- stats::aggregate(ARI ~ data.name + n + fmethod + k_policy, data = g, FUN = mean)
  D <- data.frame(data.name = z$data.name, n = as.integer(z$n),
                  method = z$fmethod, ARI = z$ARI,
                  k_policy = ifelse(z$k_policy == "true", "k = true", "k = estimated"),
                  stringsAsFactors = FALSE)
} else {
  f_true <- file.path("..", "results", "ncurve", "summary", "ncurve_by_dataset_k-true.csv")
  f_est  <- file.path("..", "results", "ncurve", "summary", "ncurve_by_dataset_k-estimated.csv")
  if (!file.exists(f_true)) stop("run results_summary_ncurve.R first: ", f_true, " not found")
  if (!isTRUE(as.numeric(allow_defective) == 1))
    stop("the stored summaries carry the uncorrected Baron runs. Point grid= at\n",
         "  grid_corrected.csv, or pass allow_defective=1 if you really mean to\n",
         "  reproduce the pre-repair figure.")
  cat("*** Figures A and B read the UNCORRECTED summaries (allow_defective=1).\n")
  D <- rbind(long(f_true, "k = true"), long(f_est, "k = estimated"))
}
D <- D[D$method %in% MAIN, ]
D$method <- factor(D$method, levels = MAIN)
D <- D[D$n > 0, ]

A <- aggregate(ARI ~ method + n + k_policy, D, mean)
lab <- A[A$n == max(A$n) & A$k_policy == "k = true", ]
pA <- ggplot(A, aes(n, ARI, colour = method, linetype = method)) +
  geom_line(linewidth = 0.6) + geom_point(size = 1.3) +
  geom_text(data = lab, aes(label = method), hjust = -0.12, size = 2.6,
            show.legend = FALSE) +
  facet_wrap(~ k_policy) +
  scale_x_log10(breaks = c(100, 300, 1000, 2500),
                minor_breaks = c(200, 500, 750, 1500)) +
  scale_colour_manual(values = PAL) + scale_linetype_manual(values = LTY) +
  coord_cartesian(clip = "off") +
  labs(x = "number of selected features (log scale)",
       y = "mean adjusted Rand index over 13 datasets") +
  theme_difs() + theme(plot.margin = margin(6, 34, 6, 6))
ggsave(file.path(out, "figA_feature_count_curve.pdf"), pA, width = 7.2, height = 3.2)

B <- D[D$n == n_focus & D$k_policy == "k = true", ]
ord <- stats::aggregate(ARI ~ data.name, B, mean)
B$data.name <- factor(B$data.name, levels = ord$data.name[order(ord$ARI)])
pB <- ggplot(B, aes(ARI, data.name, colour = method, shape = method)) +
  geom_line(aes(group = data.name), colour = "grey80", linewidth = 0.4) +
  geom_point(size = 2) +
  scale_colour_manual(values = PAL) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = sprintf("adjusted Rand index at %d features (k = true)", n_focus),
       y = NULL) +
  theme_difs()
ggsave(file.path(out, "figB_per_dataset.pdf"), pB, width = 5.6, height = 4.2)

fa  <- file.path("..", "results", "abl")
rds <- list.files(fa, pattern = "\\.rds$", full.names = TRUE)
A   <- NULL
if (nzchar(abl) && file.exists(abl)) {
  cat("Figure C reads the CORRECTED ablation.\n"); stamp(abl)
  x <- utils::read.csv(abl, stringsAsFactors = FALSE)
  x <- x[is.finite(x$ARI), ]
  A <- data.frame(data.name = x$data.name, arm = x$fmethod,
                  n = as.integer(x$n), ARI = x$ARI, stringsAsFactors = FALSE)
} else if (length(rds)) {
  if (!isTRUE(as.numeric(allow_defective) == 1))
    stop("the raw ablation .rds carry the uncorrected Baron DIFS runs. Point abl=\n",
         "  at abl_corrected.csv, or pass allow_defective=1.")
  cat("*** Figure C reads the UNCORRECTED .rds (allow_defective=1).\n")
  A <- do.call(rbind, lapply(rds, function(f) {
    x <- readRDS(f)
    data.frame(data.name = x$data.name, arm = x$feature.selection.method,
               n = x$n_requested, ARI = x$ARI, stringsAsFactors = FALSE)
  }))
}
if (is.null(A)) {
  message("ablation results not found (abl=", abl, ", ", fa, "); skipping Figure C")
} else {
  W <- stats::reshape(A, idvar = c("data.name", "n"), timevar = "arm",
                      direction = "wide")
  names(W) <- sub("^ARI\\.", "", names(W))
  need <- c("DIFS", "DIFS_stage1only", "GateOnly")
  if (!all(need %in% names(W))) stop("ablation lacks arms: ",
                                     paste(setdiff(need, names(W)), collapse = ", "))
  ok <- stats::complete.cases(W[, need])
  message("Figure C uses ", sum(ok), " of ", nrow(W), " complete triples, ",
          length(unique(W$data.name[ok])), " datasets")
  W <- W[ok, ]
  W$ranking <- W$DIFS_stage1only - W$GateOnly
  W$stage2  <- W$DIFS - W$DIFS_stage1only
  W$indep   <- !W$data.name %in% NOT_INDEPENDENT
  utils::write.csv(W, file.path(out, "figC_decomposition.csv"), row.names = FALSE)

  PARTS <- c(ranking = "dip ranking", stage2 = "stage II")
  PC    <- c("dip ranking" = "#0072B2", "stage II" = "#D55E00")   # validated pair
  PS    <- c("dip ranking" = 16,        "stage II" = 15)

  lng <- function(d) do.call(rbind, lapply(names(PARTS), function(k)
    data.frame(data.name = d$data.name, n = d$n, part = PARTS[[k]],
               value = d[[k]], indep = d$indep, stringsAsFactors = FALSE)))

  a <- lng(W[W$n == min(W$n), ])
  tot <- stats::aggregate(value ~ data.name, a, sum)
  a$data.name <- factor(a$data.name, levels = tot$data.name[order(tot$value)])
  a$part <- factor(a$part, levels = PARTS)
  seg <- stats::reshape(a[, c("data.name", "part", "value")], idvar = "data.name",
                        timevar = "part", direction = "wide")
  names(seg) <- c("data.name", "lo", "hi")
  pA <- ggplot(a, aes(value, data.name)) +
    geom_vline(xintercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_segment(data = seg, aes(x = lo, xend = hi, y = data.name, yend = data.name),
                 colour = "grey85", linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(aes(colour = part, shape = part), size = 2.1) +
    scale_colour_manual(values = PC) + scale_shape_manual(values = PS) +
    labs(x = sprintf("contribution to adjusted Rand index (at %d features)",
                     min(W$n)), y = NULL,
         title = "A   Per dataset") +
    theme_difs() + theme(plot.title = element_text(face = "bold", hjust = 0))

  b1 <- stats::aggregate(value ~ n + part, lng(W), mean); b1$set <- "all 13 datasets"
  b2 <- stats::aggregate(value ~ n + part, lng(W[W$indep, ]), mean)
  b2$set <- "9 independent real datasets"
  b <- rbind(b1, b2)
  b$part <- factor(b$part, levels = PARTS)
  pB <- ggplot(b, aes(n, value, colour = part, linetype = set)) +
    geom_hline(yintercept = 0, colour = "grey55", linewidth = 0.4) +
    geom_line(linewidth = 0.6) +
    geom_point(aes(shape = part), size = 1.4, show.legend = FALSE) +
    scale_x_log10(breaks = sort(unique(W$n))) +
    scale_colour_manual(values = PC) + scale_shape_manual(values = PS) +
    scale_linetype_manual(values = c("all 13 datasets" = "solid",
                                     "9 independent real datasets" = "dashed")) +
    labs(x = "number of selected features (log scale)",
         y = "mean contribution to ARI",
         title = "B   Mean, and its leave-out sensitivity") +
    theme_difs() + theme(plot.title = element_text(face = "bold", hjust = 0))

  ggsave(file.path(out, "figC_ablation_perdataset.pdf"), pA, width = 5.2, height = 3.6)
  ggsave(file.path(out, "figC_ablation_mean.pdf"),       pB, width = 5.0, height = 3.4)
  message("Figure C written as two files; combine them side by side in LaTeX ",
          "(subfigure) rather than with patchwork, which is not installed here.")
}

cat("\nwritten to: ", normalizePath(out), "\n", sep = "")
cat("\nSTILL TO MAKE (needs the simulation output, not the benchmark):\n")
cat("  Fig D  empirical type-I error by null shape x expression level (E2).\n")
cat("         This is the figure Reviewer 2's major comment 1 asks for by name,\n")
cat("         and it is the only one of the set we cannot draw from the grid.\n")
