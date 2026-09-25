
args <- commandArgs(TRUE)
for (a in args) {
  kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
  if (length(kv) == 3L) assign(kv[2],
    if (grepl("^-?[0-9.]+$", kv[3])) as.numeric(kv[3]) else kv[3], envir = globalenv())
}
if (!exists("B",       inherits = FALSE)) B       <- 20000
if (!exists("out",     inherits = FALSE)) out     <- "../results/dip_null_tables.rds"
if (!exists("seed",    inherits = FALSE)) seed    <- 20260910
if (!exists("trunc_q", inherits = FALSE)) trunc_q <- 0.5
if (!exists("refs",    inherits = FALSE)) refs    <- "normal,tnormal,uniform"

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
cat("libPaths:\n"); cat(paste0("  ", .libPaths(), collapse = "\n"), "\n")

if (!requireNamespace("diptest", quietly = TRUE))
  stop("diptest not found on:\n  ", paste(.libPaths(), collapse = "\n  "),
       "\nDIFS_LIB is '", Sys.getenv("DIFS_LIB", "<unset>"), "'. ",
       "Did you `source ../tools/difs_env.sh` first?")
suppressPackageStartupMessages(library(diptest))
cat("diptest ", as.character(utils::packageVersion("diptest")), " from ",
    dirname(system.file(package = "diptest")), "\n\n", sep = "")

DIFS_N_GRID <- c(35, 40, 45, 50, 60, 70, 85, 100, 120, 150, 180, 220, 270,
                 330, 400, 500, 620, 780, 1000, 1300, 1700, 2200, 3000,
                 4000, 5500, 7500, 10000)

dip_stat <- function(x) diptest::dip(sort(x[is.finite(x)]))

dip_null_sample <- function(n, B, ref, trunc_q = 0.5) {
  gen <- switch(ref,
    uniform = function() stats::runif(n),
    normal  = function() stats::rnorm(n),
    tnormal = function() stats::qnorm(stats::runif(n, trunc_q, 1)),
    stop("unknown ref: ", ref))
  vapply(seq_len(B), function(i) dip_stat(gen()), numeric(1))
}

build_dip_null_table <- function(n_grid = DIFS_N_GRID, B = 10000, ref = "normal",
                                 trunc_q = 0.5,
                                 probs = c(seq(0.0005, 0.9990, length.out = 700),
                                           0.9995, 0.9999)) {
  Q <- matrix(NA_real_, nrow = length(n_grid), ncol = length(probs),
              dimnames = list(as.character(n_grid), NULL))
  for (i in seq_along(n_grid)) {
    t0 <- proc.time()[["elapsed"]]
    d <- dip_null_sample(n_grid[i], B = B, ref = ref, trunc_q = trunc_q)
    Q[i, ] <- stats::quantile(d, probs = probs, names = FALSE, type = 7)
    message(sprintf("  [%s] n = %6d  (%2d/%2d)  %.1fs", ref, n_grid[i], i,
                    length(n_grid), proc.time()[["elapsed"]] - t0))
  }
  structure(list(n_grid = n_grid, probs = probs, Q = Q, ref = ref, B = B,
                 trunc_q = if (ref == "tnormal") trunc_q else NA_real_,
                 built = Sys.time()),
            class = "dip_null_table")
}

set.seed(seed)
tabs <- list()
for (r in trimws(strsplit(refs, ",")[[1]])) {
  message("=== ", r, " (B = ", B, ") ===")
  tabs[[r]] <- build_dip_null_table(B = B, ref = r, trunc_q = trunc_q)
}
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
saveRDS(tabs, out)
cat("\nwritten: ", normalizePath(out), "\n", sep = "")
cat("references: ", paste(names(tabs), collapse = ", "), "\n", sep = "")
cat("n grid    : ", paste(DIFS_N_GRID, collapse = " "), "\n", sep = "")
cat("B         : ", B, "\n", sep = "")
cat("\nPoint DIFS_NULL_TABLES at this file (difs_null_path() reads it).\n")
