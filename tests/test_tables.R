# Validate bundled summaries without claiming recomputation from cell labels.
s <- read.csv("results/supp/S0_per_run.csv",stringsAsFactors=FALSE)
key_cols <- c("dataset","feature_method","clustering_method","k_policy","n_requested")
key <- function(x,cols=key_cols) do.call(paste,c(x[,cols,drop=FALSE],sep="|"))
stopifnot(nrow(s)==2601L,sum(s$is_analysis_row)==2489L,
  sum(s$result_set=="ncurve")==2320L,sum(s$result_set=="abl")==117L,
  sum(s$result_set=="hart")==52L,sum(s$result_set=="baronfix")==112L,
  !anyDuplicated(key(s,c("result_set",key_cols,"gate_rule","seed"))))
repaired <- s[s$baron_repaired,]
stopifnot(nrow(repaired)==115L,all(repaired$n_cells_total==8569L),
  all(repaired$n_cells_dropped==0L),all(repaired$n_cells_unlabelled==0L),
  all(repaired$source_result_set=="baronfix"),
  all(is.finite(as.matrix(repaired[,c("ARI","FM","Jaccard","Purity")]))))
main <- s[s$result_set=="ncurve",]
for (i in seq_along(c("ARI","FM","Jaccard","Purity"))) {
  metric <- c("ARI","FM","Jaccard","Purity")[i]
  x <- read.csv(sprintf("results/supp/S%d_%s_per_run.csv",i,metric),stringsAsFactors=FALSE)
  stopifnot(nrow(x)==2320L,!anyDuplicated(key(x)),setequal(key(x),key(main)))
  j <- match(key(x),key(main));stopifnot(all(abs(x$value-main[[metric]][j])<1e-12))
}
missing <- read.csv("results/supp/S5_missing_configurations.csv",stringsAsFactors=FALSE)
grid <- expand.grid(lapply(main[,key_cols],unique),stringsAsFactors=FALSE)
stopifnot(nrow(grid)==2340L,nrow(missing)==20L,
  setequal(key(missing),setdiff(key(grid),key(main))))
cat("PASS: table identities, provenance, all four corrected metrics and 20 missing configurations.\n")
