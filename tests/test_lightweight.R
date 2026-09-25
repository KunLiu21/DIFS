# These tests do not require SC3 and do not claim a full DIFS execution.
source("R/difs.R")
e <- difs_load_engine(full=FALSE)
fails <- function(expr) stopifnot(inherits(try(force(expr),silent=TRUE),"try-error"))
stopifnot(e$difs_gate_threshold(180,k=3,rule="submitted")==35,
          e$difs_gate_threshold(8569,k=14,rule="submitted")==0.05*8569,
          e$difs_gate_threshold(1000,k=2,rule="submitted")==e$difs_gate_threshold(1000,k=20,rule="submitted"))
d <- readRDS("demo/demo_data.rds")
r1 <- e$difs_stage1_ranking(d$logcounts,gate_rule="submitted",gate_k=3)
r2 <- e$difs_stage1_ranking(d$logcounts,gate_rule="submitted",gate_k=3)
stopifnot(identical(r1,r2),length(r1)==422L,!anyDuplicated(r1))
# Directly compare loaded bodies/formals with the canonical benchmark nodes.
for (file in c("Functions_mixture_auto.R","Functions_controlled.R","difs_sc3_complete.R")) {
  for (expr in parse(file.path("benchmark",file),keep.source=FALSE)) {
    if (is.call(expr) && as.character(expr[[1]]) %in% c("<-","=") && is.symbol(expr[[2]]) &&
        is.call(expr[[3]]) && identical(expr[[3]][[1]],as.name("function"))) {
      name <- as.character(expr[[2]])
      if (name %in% c("num_features_decider","sc3_clustering") && file=="Functions_mixture_auto.R") next
      if (exists(name,e,inherits=FALSE) && name!="difs_null_path") {
        target <- eval(expr[[3]],e)
        if (!identical(body(e[[name]]),body(target)) || !identical(formals(e[[name]]),formals(target)))
          stop("Loaded function differs from its canonical source: ",file," / ",name)
      }
    }
  }
}
difs_validate_counts(d$counts)
z <- d$counts; colnames(z)[2] <- colnames(z)[1];fails(difs_validate_counts(z))
z <- d$counts;z[1,1] <- NA;fails(difs_validate_counts(z))
z <- d$counts;z[,1] <- 0;fails(difs_validate_counts(z))
z <- d$counts;z[1,1] <- -1;fails(difs_validate_counts(z))
# Fractional estimates are valid benchmark inputs, both dense and sparse.
z <- d$counts;z[1,1] <- 0.25;before <- z;difs_validate_counts(z);stopifnot(identical(z,before))
zs <- methods::as(Matrix::Matrix(z,sparse=TRUE),"dgCMatrix")
before <- zs;difs_validate_counts(zs);stopifnot(identical(zs,before))
for (invalid in c(NA_real_,NaN,Inf,-Inf,-0.25)) {
  zd <- z;zd[1,1] <- invalid;fails(difs_validate_counts(zd))
  zz <- zs;zz[1,1] <- invalid;fails(difs_validate_counts(zz))
}
fails(difs_check_integer(2.5,"k"));fails(difs_require("DIFS_intentionally_missing_package"))
truth <- setNames(c("a","a","b","b"),letters[1:4])
pred <- setNames(c("x","x","x","y"),letters[1:4])
stopifnot(isTRUE(all.equal(unname(difs_metrics(truth,pred)),c(0,1/sqrt(6),0.25,0.75))))
stopifnot(identical(difs_metrics(truth,pred),difs_metrics(truth[4:1],pred)))
stopifnot(all(difs_metrics(truth,truth)==1))
bad <- pred;names(bad)[1] <- "missing";fails(difs_metrics(truth,bad))
bad <- pred;bad[1] <- NA;fails(difs_metrics(truth,bad))
source("benchmark/difs_signrank.R")
stopifnot(difs_p(c(-1,2,-3,4,5,6,7,8,0))==0.0546875)
fails(difs_p(c(1,NA)))
cat("PASS: gate policy, deterministic ranking, canonical function bodies, input/label guards, hand-calculated metrics and exact signed-rank.\n")
