
args <- commandArgs(TRUE)
for (a in args) {
  kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
  if (length(kv) != 3L) stop("invalid argument: ", a)
  key <- kv[2L]; val <- kv[3L]
  if (!key %in% c("preset","quick","out_dir","seed","exp","gate","norm_method","k")) stop("unknown argument: ", key)
  val <- sub('^"(.*)"$', "\\1", val); val <- sub("^'(.*)'$", "\\1", val)
  v <- if (grepl("^(TRUE|FALSE|T|F)$", val)) as.logical(val)
       else if (grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", val)) as.numeric(val)
       else val
  assign(key, v, envir = globalenv())
}

if (!exists("preset", inherits = FALSE))  preset  <- NULL
if (!exists("quick", inherits = FALSE))   quick   <- FALSE
if (!is.null(preset) && preset == "quick") quick <- TRUE
if (is.null(preset)) preset <- if (quick) "quick" else "full"
if (!exists("out_dir", inherits = FALSE)) out_dir <- "../results/simulation"
if (!exists("seed", inherits = FALSE))    seed    <- 20260905
if (!exists("exp", inherits = FALSE))     exp     <- NULL
if (!exists("gate", inherits = FALSE))    gate    <- "implemented"
if (!exists("norm_method", inherits = FALSE)) norm_method <- "lognorm"
if (exists("k", inherits = FALSE)) exp <- paste0("E", k)

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
       "\nDIFS_LIB is '", Sys.getenv("DIFS_LIB", "<unset>"),
       "'.  Did you `source tools/difs_env.sh` first?")
cat("diptest ", as.character(utils::packageVersion("diptest")), "\n\n", sep = "")

here <- tryCatch(dirname(normalizePath(sub("^--file=", "",
         grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))),
         error = function(e) ".")
if (is.na(here) || !nzchar(here)) here <- "."
source(file.path(here, "difs_sim_core.R"))
source(file.path(here, "difs_sim_experiments.R"))

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GATE <- if (is.numeric(gate)) gate else
  switch(gate, documented = DIFS_GATE_AS_DOCUMENTED,
               implemented = DIFS_GATE_AS_IMPLEMENTED,
               stop("gate must be \"documented\", \"implemented\", or a number"))
message(sprintf("stage I expression gate = %.4f (%s);  normalisation = %s",
                GATE, if (is.numeric(gate)) "explicit" else gate, norm_method))

tab_file <- file.path(out_dir, "dip_null_tables.rds")
B_tab   <- switch(preset, quick = 800, medium = 3000, 5000)
n_grid  <- switch(preset,
  quick  = c(20, 35, 60, 110, 200, 400, 800, 1600, 3000),
  medium = DEFAULT_N_GRID[DEFAULT_N_GRID <= 5000],
  DEFAULT_N_GRID)
if (file.exists(tab_file)) {
  message("loading cached null tables: ", tab_file)
  TABS <- readRDS(tab_file)
} else {
  message("building null tables (B = ", B_tab, ") ...")
  set.seed(seed)
  TABS <- list(
    normal  = build_dip_null_table(n_grid, B_tab, ref = "normal"),
    tnormal = build_dip_null_table(n_grid, B_tab, ref = "tnormal"),
    uniform = build_dip_null_table(n_grid, B_tab, ref = "uniform"))
  saveRDS(TABS, tab_file)
}

cfg <- if (preset == "medium") {
  list(E2 = list(n_cells = c(200, 500, 1000), n_genes = 300,
                 mus = c(1, 3, 10, 30, 100)),
       E3 = list(n_cells = c(200, 500), n_genes = 200,
                 pis = c(0.05, 0.10, 0.20, 0.35, 0.50),
                 log2fcs = c(2, 3, 4, 6), mus = c(3, 10, 30)),
       E4 = list(n_cells = c(300, 600), n_genes = 1000, n_bg_genes = 2000,
                 reps = 3, log2fcs = c(2, 3, 4, 5, 6, 8), disps = c(0.2, 0.5)),
       E5 = list(n_cells = 600, n_genes = 1500, n_bg_genes = 1500, B = 2000),
       E6 = list(n_cells = 600, n_genes = 1500, n_bg_genes = 1500, topN = 1000),
       E7 = list(n_cells = 800, n_genes = 1500, n_bg_genes = 1500, reps = 3))
} else if (quick) {
  list(E2 = list(n_cells = c(200, 500), n_genes = 150, mus = c(1, 10, 100)),
       E3 = list(n_cells = c(200, 500), n_genes = 100,
                 pis = c(0.05, 0.20, 0.50), log2fcs = c(2, 4), mus = c(10)),
       E4 = list(n_cells = c(600), n_genes = 800, reps = 1,
                 log2fcs = c(3, 6), disps = c(0.5), n_bg_genes = 400),
       E5 = list(n_cells = 300, n_genes = 300, B = 500, n_bg_genes = 400),
       E6 = list(n_cells = 300, n_genes = 600, n_bg_genes = 400, topN = 300),
       E7 = list(n_cells = 400, n_genes = 600, reps = 2, n_bg_genes = 400))
} else {
  list(E2 = list(n_cells = c(200, 500, 1000), n_genes = 1000),
       E3 = list(n_cells = c(200, 500, 1000), n_genes = 400,
                 pis = c(0.05, 0.10, 0.20, 0.35, 0.50), log2fcs = c(1, 2, 3, 4)),
       E4 = list(n_cells = c(300, 600, 1200), n_genes = 3000, reps = 5),
       E5 = list(n_cells = 600, n_genes = 1500, B = 2000),
       E6 = list(n_cells = 600, n_genes = 2000, topN = 1000),
       E7 = list(n_cells = 800, n_genes = 2000, reps = 3))
}

for (nm in names(cfg)) {
  if (nm %in% c("E2", "E3", "E4", "E5")) cfg[[nm]]$min.expression <- GATE
  cfg[[nm]]$norm_method <- norm_method
}

run <- function(id) is.null(exp) || !is.character(exp) || id %in% strsplit(paste(exp, collapse = ","), ",")[[1]]

if (run("E1")) { message("== E1 location-scale invariance =="); print(utils::head(
  E1_invariance(out_dir = out_dir, tab_norm = TABS$normal)$monotonicity)) }

if (run("E2")) { message("== E2 null calibration =="); print(
  do.call(E2_null_calibration, c(cfg$E2, list(tab_norm = TABS$normal,
    tab_unif = TABS$uniform, out_dir = out_dir)))) }

if (run("E3")) { message("== E3 power =="); print(
  do.call(E3_power, c(cfg$E3, list(tab_norm = TABS$normal, out_dir = out_dir)))) }

if (run("E4")) { message("== E4 ranking quality =="); print(
  do.call(E4_ranking, c(cfg$E4, list(tab_norm = TABS$normal, out_dir = out_dir)))) }

if (run("E5")) { message("== E5 ties and runtime =="); print(
  do.call(E5_ties_runtime, c(cfg$E5, list(tab_norm = TABS$normal, out_dir = out_dir)))) }

if (run("E6")) { message("== E6 gate and reference =="); print(
  do.call(E6_gate_reference, c(cfg$E6, list(tab_norm = TABS$normal,
    tab_tnorm = TABS$tnormal, out_dir = out_dir)))) }

if (run("E7")) { message("== E7 stage II robustness =="); print(
  do.call(E7_stage2_robustness, c(cfg$E7, list(out_dir = out_dir)))) }

message("done. results in ", normalizePath(out_dir))
