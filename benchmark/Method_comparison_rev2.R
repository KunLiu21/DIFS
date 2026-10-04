# Method_comparison_rev2.R: revision-2 analysis source.
# Source: code/as_run/06_revision2_reps_krange_abl2_runtime/Method_comparison_rev2.R
# Usage, scope and limitations: benchmark/REVISION2.md.
difs_args <- function(args, known, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    if (!kv[2] %in% known)
      stop("unknown argument '", kv[2], "'. This script accepts: ", paste(known, collapse = ", "))
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
  invisible(NULL)
}
set <- ""; k <- NA
difs_args(commandArgs(TRUE), c("set", "k"))
if (!nzchar(set) || is.na(k)) stop("set= and k= are both required")

source("Functions_controlled.R")
for (f in c("difs_sc3_complete.R", "difs_rev2_overlay.R"))
  if (!file.exists(f)) stop(f, " must be in the code directory")
source("difs_sc3_complete.R")
source("difs_rev2_overlay.R")

load(sprintf("../scenarios/scenarios_%s.Rdata", set))
if (k < 1 || k > nrow(scenarios)) stop("k = ", k, " outside 1..", nrow(scenarios))
P <- scenarios[k, ]
task_id <- k; rm(k)
col <- function(nm, def) if (nm %in% names(P) && !is.na(P[[nm]])) P[[nm]] else def
data.path <- as.character(P$data.path)
data.name <- gsub(".*/([^.]*).*", "\\1", data.path)
arm       <- as.character(col("arm", P$feature.selection.method))
stage1    <- as.character(col("stage1_source", "dip"))
final_set <- as.character(col("final_set", "mixed"))
cseed     <- col("cluster_seed", NA)
sub_rep   <- as.integer(col("subsample_rep", 0L))
sub_frac  <- as.numeric(col("subsample_frac", 1))
k_offset  <- col("k_offset", NA)
gate_rule <- as.character(col("gate_rule", "submitted"))

out <- prepDir(paste0("../results/", set))
fname <- sprintf("%s_%s_%s_n%04d_k-%s%s_cs-%s_sub-%d_seed-%s.rds",
                 data.name, arm, gsub(" ", "", P$method), P$n_features, P$k_policy,
                 if (is.na(k_offset)) "" else sprintf("%+d", as.integer(k_offset)),
                 if (is.na(cseed)) "hist" else as.character(cseed), sub_rep, P$seed)

seu <- readRDS(data.path)
lab <- seu$trueclass
if (sub_rep > 0L) {
  if (!(sub_frac > 0 && sub_frac < 1)) stop("subsample_frac must be in (0, 1)")
  set.seed(10000L + sub_rep)
  idx <- unlist(lapply(split(seq_along(lab), addNA(factor(lab)), drop = TRUE), function(ii)
    ii[sample.int(length(ii), max(1L, round(sub_frac * length(ii))))]))
  seu <- seu[, sort(idx)]
  tmp <- file.path(tempdir(), paste0(data.name, ".rds"))   # keeps the dataset name
  saveRDS(seu, tmp); data.path <- tmp
  lab <- seu$trueclass
}
n_cells_used <- ncol(seu)
k_true_labels <- length(unique(lab[!is.na(lab)]))
assign("DIFS_CURRENT_SEU", if (identical(stage1, "hvg")) seu else NULL, envir = globalenv())
rm(seu); invisible(gc())

k_policy_run <- as.character(P$k_policy)
k_target <- NA_integer_
if (!is.na(k_offset)) {
  k_target <- k_true_labels + as.integer(k_offset)
  if (k_target < 2L) {
    saveRDS(list(data.name = data.name, arm = arm, status = "skipped_k_below_2",
                 clustering.method = as.character(P$method), n_requested = P$n_features,
                 cluster_seed = if (is.na(cseed)) NA_integer_ else as.integer(cseed),
                 subsample_rep = sub_rep, seed = P$seed,
                 k_true = k_true_labels, k_offset = as.integer(k_offset)), file.path(out, fname))
    cat("task", task_id, "skipped:", data.name, "k_true", k_true_labels,
        "offset", k_offset, "gives k < 2\n"); quit(save = "no", status = 0)
  }
  options(difs.k_override = k_target); k_policy_run <- "estimated"
}
options(difs.stage1_source = stage1, difs.final_set = final_set,
        difs.cluster_seed = if (is.na(cseed)) NULL else as.integer(cseed))

res <- Run.method.comparison.controlled(
  data.path = data.path, clustering.method = as.character(P$method),
  feature.selection.method = as.character(P$feature.selection.method),
  k_policy = k_policy_run, n_policy = as.character(P$n_policy),
  n_fixed = P$n_features, gate_rule = gate_rule,
  min.expression = P$min.expression, seed = P$seed)
if (!is.na(k_target) && !identical(as.integer(res$k_used_clustering), k_target))
  stop("k override did not reach the clustering: used ", res$k_used_clustering,
       ", intended ", k_target)

res$runner_fmethod   <- res$feature.selection.method
res$feature.selection.method <- arm
res$arm <- arm; res$stage1_source <- stage1; res$final_set <- final_set
res$cluster_seed <- if (is.na(cseed)) NA_integer_ else as.integer(cseed)
res$subsample_rep <- sub_rep; res$subsample_frac <- if (sub_rep > 0L) sub_frac else 1
res$n_cells_used <- n_cells_used; res$k_true_labels <- k_true_labels
res$k_offset <- if (is.na(k_offset)) NA_integer_ else as.integer(k_offset)
if (!is.na(k_offset)) res$k_policy <- "true_offset"
res$n_requested <- P$n_features
res$stage2only_fell_back <- identical(final_set, "stage2only") &&
  isTRUE(res$difs_stage2_n == 0)
res$status <- "ok"
saveRDS(res, file.path(out, fname))

cat("task", task_id, "done:", data.name, arm, as.character(P$method),
    "| n", P$n_features, "used", res$n_features,
    "| k", res$k_used_clustering, "found", res$n_clusters_found,
    "| cseed", if (is.na(cseed)) "hist" else cseed, "| sub", sub_rep,
    "| cells", n_cells_used, "| stage2_n", res$difs_stage2_n,
    if (isTRUE(res$stage2only_fell_back)) "| STAGE II EMPTY -> fell back to stage I" else "",
    "| ARI =", round(res$ARI, 4), "\n")
