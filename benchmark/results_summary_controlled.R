
args <- commandArgs(TRUE); for (a in args) eval(parse(text = a))
if (!exists("resolution", inherits = FALSE)) resolution <- "controlled_comparison"
path <- file.path("../results", resolution)
out  <- file.path(path, "summary")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

files <- list.files(path, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("no result files in ", path)
cat("reading", length(files), "result files\n")

flat <- function(x) {
  keep <- c("data.name", "clustering.method", "feature.selection.method",
            "k_policy", "n_policy", "seed", "min.expression",
            "k_true", "k_estimated", "k_used_feature_selection",
            "k_used_clustering", "n_features", "difs_stage1_n", "difs_stage2_n",
            "difs_ratio", "ARI", "Purity", "Jaccard", "FM", "elapsed_sec")
  as.data.frame(lapply(keep, function(f) {
    v <- x[[f]]; if (is.null(v) || !length(v)) NA else unname(v[1])
  }), col.names = keep, stringsAsFactors = FALSE)
}
d <- do.call(rbind, lapply(files, function(f) tryCatch(flat(readRDS(f)),
       error = function(e) { message("skipped ", basename(f), ": ", e$message); NULL })))
d$cell <- paste(d$k_policy, d$n_policy, sep = " / ")
utils::write.csv(d, file.path(out, "raw_results.csv"), row.names = FALSE)

msg <- function(t) cat("\n========== ", t, " ==========\n")

s1 <- stats::aggregate(cbind(ARI, Purity, Jaccard, FM, n_features) ~
                         feature.selection.method + k_policy + n_policy, d, mean)
s1 <- s1[order(s1$k_policy, s1$n_policy, -s1$ARI), ]
utils::write.csv(s1, file.path(out, "S1_headline.csv"), row.names = FALSE)
msg("S1  mean ARI by feature selection method")
w <- stats::reshape(
  transform(s1, cell = paste(k_policy, n_policy, sep = "/"))[
    , c("feature.selection.method", "cell", "ARI")],
  idvar = "feature.selection.method", timevar = "cell", direction = "wide")
print(w[order(-w[[2]]), ], row.names = FALSE, digits = 3)

controls <- c("Random", "SC3_random1000")
cells <- unique(d$cell)
s2 <- do.call(rbind, lapply(cells, function(cc) {
  sub <- subset(s1, paste(k_policy, n_policy, sep = " / ") == cc)
  difs <- sub$ARI[sub$feature.selection.method == "DIFS"]
  comp <- subset(sub, !feature.selection.method %in% c("DIFS", controls))
  best <- if (nrow(comp)) comp$feature.selection.method[which.max(comp$ARI)] else NA
  data.frame(cell = cc,
             ARI_DIFS = if (length(difs)) difs else NA,
             best_competitor = best,
             ARI_best_competitor = if (nrow(comp)) max(comp$ARI) else NA,
             margin = if (length(difs) && nrow(comp)) difs - max(comp$ARI) else NA,
             DIFS_rank = if (length(difs))
               rank(-sub$ARI)[sub$feature.selection.method == "DIFS"] else NA)
}))
base <- s2$margin[s2$cell == "original / original"]
s2$margin_change_vs_submitted <- if (length(base)) s2$margin - base else NA
utils::write.csv(s2, file.path(out, "S2_attribution.csv"), row.names = FALSE)
msg("S2  DIFS margin over the best competing method, by configuration")
print(s2, row.names = FALSE, digits = 3)
cat("\n'original / original' reproduces the submitted setup.\n",
    "A margin that stays positive under 'estimated / tuned_all' is the result\n",
    "the revision can claim.\n", sep = "")

s3 <- stats::aggregate(ARI ~ data.name + feature.selection.method + k_policy +
                         n_policy + clustering.method, d,
                       function(z) c(mean = mean(z), sd = stats::sd(z), n = length(z)))
s3 <- do.call(data.frame, s3)
names(s3)[(ncol(s3) - 2):ncol(s3)] <- c("ARI_mean", "ARI_sd", "n_seeds")
utils::write.csv(s3, file.path(out, "S3_per_dataset.csv"), row.names = FALSE)
msg("S3  per dataset, mean +/- sd over seeds (first 20 rows)")
print(utils::head(s3[order(s3$data.name, -s3$ARI_mean), ], 20), row.names = FALSE, digits = 3)

s4 <- do.call(rbind, lapply(split(d, list(d$data.name, d$clustering.method, d$cell),
                                  drop = TRUE), function(g) {
  a <- stats::aggregate(ARI ~ feature.selection.method, g, mean)
  a$rank <- rank(-a$ARI, ties.method = "min")
  a$data.name <- g$data.name[1]; a$clustering.method <- g$clustering.method[1]
  a$cell <- g$cell[1]; a
}))
utils::write.csv(s4, file.path(out, "S4_ranking.csv"), row.names = FALSE)
msg("S4  how often each method ranks first, by configuration")
print(with(subset(s4, rank == 1), table(feature.selection.method, cell)))

s5 <- subset(d, feature.selection.method == "DIFS",
             select = c(data.name, clustering.method, k_policy, n_policy, seed,
                        k_used_feature_selection, n_features, difs_stage1_n,
                        difs_stage2_n, difs_ratio, elapsed_sec))
utils::write.csv(s5, file.path(out, "S5_difs_internals.csv"), row.names = FALSE)
msg("S5  DIFS internals -- feature counts and chosen mixing ratio")
print(utils::head(s5[order(s5$data.name, s5$k_policy), ], 20), row.names = FALSE)

msg("feature set sizes actually used (mean)")
print(stats::reshape(
  stats::aggregate(n_features ~ feature.selection.method + n_policy, d, mean),
  idvar = "feature.selection.method", timevar = "n_policy", direction = "wide"),
  row.names = FALSE, digits = 4)

cat("\nsummary written to ", normalizePath(out), "\n")
