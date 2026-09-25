difs_args <- function(args, env = parent.frame()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) stop("cannot parse argument: ", a)
    v <- kv[3]; num <- suppressWarnings(as.numeric(v))
    assign(kv[2], if (nzchar(v) && !is.na(num)) num else v, envir = env)
  }
}
set <- "ncurve"
difs_args(commandArgs(TRUE))
source("difs_signrank.R")   ## exact paired test; see that file for why not wilcox.test

dir_in  <- file.path("..", "results", set)
dir_out <- file.path(dir_in, "summary")
dir.create(dir_out, recursive = TRUE, showWarnings = FALSE)

fs <- list.files(dir_in, pattern = "\\.rds$", full.names = TRUE)
cat("reading", length(fs), "result files\n")
ls_ <- lapply(fs, readRDS)
nms <- unique(unlist(lapply(ls_, names)))
d <- do.call(rbind, lapply(ls_, function(x) {
  v <- setNames(rep(list(NA), length(nms)), nms)
  for (n in names(x)) v[[n]] <- if (is.null(x[[n]])) NA else x[[n]]
  as.data.frame(v, stringsAsFactors = FALSE)
}))
d$n_label <- ifelse(d$n_requested == 0, "gate_n", as.character(d$n_requested))
if (is.null(d$k_policy)) d$k_policy <- "unknown"
write.csv(d, file.path(dir_out, "ncurve_raw.csv"), row.names = FALSE)

cat("\nk policies present:", paste(sort(unique(d$k_policy)), collapse = ", "), "\n")

if (!is.null(d$k_hit)) {
  cat("\n========== N0  operating point: fraction of runs that hit the target k ==========\n")
  h <- aggregate(cbind(k_hit = as.numeric(k_hit)) ~ feature.selection.method +
                   clustering.method + k_policy, d, mean)
  hw <- reshape(h, idvar = c("feature.selection.method", "k_policy"),
                timevar = "clustering.method", direction = "wide")
  names(hw) <- sub("k_hit\\.", "", names(hw))
  num <- vapply(hw, is.numeric, logical(1)); hw[num] <- round(hw[num], 3)
  print(hw, row.names = FALSE)
  if (!is.null(d$n_cells_dropped)) {
    dr <- aggregate(n_cells_dropped ~ feature.selection.method + k_policy, d, mean)
    dr$n_cells_dropped <- round(dr$n_cells_dropped, 1)
    cat("\n  mean cells left unassigned (excluded from scoring):\n")
    print(dr[dr$n_cells_dropped > 0, ], row.names = FALSE)
    if (!any(dr$n_cells_dropped > 0)) cat("    none\n")
  }
}

for (kp in sort(unique(d$k_policy))) {
 dk <- d[d$k_policy == kp, ]
 cat("\n\n##################  k_policy = ", kp, "  ##################\n", sep = "")
 cat("\n========== N1  mean ARI by feature count (over datasets, clustering) ==========\n")
 a <- aggregate(ARI ~ feature.selection.method + n_requested, dk, mean)
 w <- reshape(a, idvar = "feature.selection.method", timevar = "n_requested",
              direction = "wide")
 names(w) <- sub("ARI\\.", "n=", names(w)); names(w) <- sub("n=0$", "n=gate_n", names(w))
 w[-1] <- round(w[-1], 3)
 print(w[order(-w[[ncol(w)]]), ], row.names = FALSE)

cat("\n========== N2  the feature count that maximises ARI, per dataset ==========\n")
best <- do.call(rbind, lapply(split(dk, list(dk$data.name, dk$feature.selection.method),
                                   drop = TRUE), function(z) {
  z2 <- aggregate(ARI ~ n_requested, z, mean)
  i <- which.max(z2$ARI)
  data.frame(data.name = z$data.name[1], method = z$feature.selection.method[1],
             n_cells = z$n_features[1] * 0 + NA, gate_n = z$gate_n[1],
             best_n = z2$n_requested[i], best_ARI = round(z2$ARI[i], 3),
             ARI_at_1000 = round(z2$ARI[z2$n_requested == 1000][1], 3),
             stringsAsFactors = FALSE)
}))
best$best_n_label <- ifelse(best$best_n == 0, "gate_n", best$best_n)
best <- best[order(best$method, best$gate_n), ]
print(best[, c("data.name", "method", "gate_n", "best_n_label", "best_ARI",
               "ARI_at_1000")], row.names = FALSE)
cat("\n  If best_n is well below 1000 on the small datasets, the fixed 1000 used\n")
cat("  throughout the submitted paper is the wrong operating point -- and the\n")
cat("  criticism applies to our own results first.\n")

cat("\n========== N3  DIFS minus GateOnly by feature count ==========\n")
cat("  (BOTH stages together, isolated from the gate -- NOT the ranking alone.\n",
    "   The ablation splits it: at 100 features the dip ranking and stage II\n",
    "   contribute about equally, and at 300 the ranking's share is negative.\n",
    "   See N8, and do not describe this number as the ranking's contribution.)\n",
    sep = "")
g <- aggregate(ARI ~ feature.selection.method + n_requested + data.name, dk, mean)
sp <- reshape(g, idvar = c("data.name", "n_requested"),
              timevar = "feature.selection.method", direction = "wide")
names(sp) <- sub("ARI\\.", "", names(sp))
if (all(c("DIFS", "GateOnly") %in% names(sp))) {
  sp$rank_gain <- round(sp$DIFS - sp$GateOnly, 3)
  gg <- aggregate(rank_gain ~ n_requested, sp, mean)
  gg$n_label <- ifelse(gg$n_requested == 0, "gate_n", gg$n_requested)
  print(gg[, c("n_label", "rank_gain")], row.names = FALSE)
}

cat("\n========== N4  DIFS margin over the best baseline, by feature count ==========\n")
base <- c("Variance", "Seurat", "FEAST")
if (all(c("DIFS", base) %in% names(sp))) {
  sp$best_base <- apply(sp[, base, drop = FALSE], 1, max, na.rm = TRUE)
  sp$margin <- round(sp$DIFS - sp$best_base, 3)
  mm <- aggregate(margin ~ n_requested, sp, mean)
  mm$n_label <- ifelse(mm$n_requested == 0, "gate_n", mm$n_requested)
  print(mm[, c("n_label", "margin")], row.names = FALSE)
  write.csv(sp, file.path(dir_out, paste0("ncurve_by_dataset_k-", kp, ".csv")),
            row.names = FALSE)
 }

 cat("\n========== N5  paired Wilcoxon, DIFS vs each method ==========\n")
 for (rng in list(all = unique(dk$n_requested), small = c(100, 200, 300))) {
   ss <- sp[sp$n_requested %in% rng, , drop = FALSE]
   cat("  --- n in {", paste(ifelse(rng == 0, "gate_n", rng), collapse = ", "), "} ---\n", sep = "")
   for (b in c("GateOnly", "Variance", "Seurat", "FEAST")) {
     if (!b %in% names(ss)) next
     keep <- !is.na(ss$DIFS) & !is.na(ss[[b]])
     if (sum(keep) < 5) next
     df <- ss$DIFS[keep] - ss[[b]][keep]
     pv <- difs_p(df)
     cat(sprintf("    DIFS vs %-9s  mean %+0.4f   win %2d/%2d   p = %.3g\n",
                 b, mean(df), sum(df > 0), length(df), pv))
   }
 }

 cat("\n========== N6  paired Wilcoxon at EACH feature count separately ==========\n")
 cat("  pairing unit = dataset (ARI first averaged over clustering methods)\n")
 cat("  read 'win w/N' before the p-value: a large mean with w ~ N/2 is two\n",
     "  datasets carrying the average, not an effect.\n", sep = "")
 bl <- c("FEAST", "Variance", "Seurat", "GateOnly", "DIFS_stage1only")
 bl <- bl[bl %in% names(sp)]
 cat(sprintf("\n%8s", "n"))
 for (b in bl) cat(sprintf("%28s", b))
 cat("\n")
 for (nn in sort(unique(sp$n_requested))) {
   ss <- sp[sp$n_requested == nn, , drop = FALSE]
   cat(sprintf("%8s", if (nn == 0) "gate_n" else nn))
   for (b in bl) {
     keep <- !is.na(ss$DIFS) & !is.na(ss[[b]])
     if (sum(keep) < 4) { cat(sprintf("%28s", "-")); next }
     df <- ss$DIFS[keep] - ss[[b]][keep]
     pv <- difs_p(df)
     cat(sprintf("   %+0.3f %2d/%2d p=%6.3f", mean(df), sum(df > 0), length(df), pv))
   }
   cat("\n")
 }

 cat("\n========== N7  DIFS minus each baseline, per dataset ==========\n")
 kmap <- tapply(dk$k_true, dk$data.name, function(z) z[1])
 gmap <- tapply(dk$gate_n, dk$data.name, function(z) z[1])
 show_n <- intersect(c(100, 200, 300, 500, 1000), unique(sp$n_requested))
 for (b in bl) {
   cat("\n--- DIFS - ", b, " ---\n", sep = "")
   cat(sprintf("%-16s%4s%8s", "dataset", "k", "gate_n"))
   for (nn in show_n) cat(sprintf("%9s", paste0("n=", nn)))
   cat("\n")
   for (x in sort(unique(sp$data.name))) {
     cat(sprintf("%-16s%4s%8s", x, kmap[[x]], gmap[[x]]))
     for (nn in show_n) {
       r <- sp[sp$data.name == x & sp$n_requested == nn, , drop = FALSE]
       v <- if (nrow(r) && !is.na(r$DIFS[1]) && !is.na(r[[b]][1]))
              r$DIFS[1] - r[[b]][1] else NA_real_
       cat(if (is.na(v)) sprintf("%9s", "--") else sprintf("%+9.3f", v))
     }
     cat("\n")
   }
 }

 if (all(c("DIFS", "DIFS_stage1only", "GateOnly") %in% names(sp))) {
   cat("\n========== N8  stage-by-stage decomposition ==========\n")
   ok3 <- !is.na(sp$DIFS) & !is.na(sp$DIFS_stage1only) & !is.na(sp$GateOnly)
   cat("  complete triples: ", sum(ok3), " of ", nrow(sp),
       " dataset x n cells\n", sep = "")
   if (sum(ok3) < nrow(sp))
     cat("  INCOMPLETE -- the arms below are not measured on the same datasets\n",
         "  as the totals elsewhere in this file.  Fill the gaps before using.\n",
         sep = "")
   sp <- sp[ok3, , drop = FALSE]
   if (!nrow(sp)) { cat("  nothing to decompose\n") } else {
   sp$rank_gain2  <- sp$DIFS_stage1only - sp$GateOnly
   sp$stage2_gain <- sp$DIFS - sp$DIFS_stage1only
   sp$total_gain  <- sp$DIFS - sp$GateOnly
   for (nn in sort(unique(sp$n_requested))) {
     ss <- sp[sp$n_requested == nn, , drop = FALSE]
     f <- function(v) {
       v <- v[!is.na(v)]
       if (length(v) < 4) return("      -")
       p <- difs_p(v)
       sprintf("%+0.4f %2d/%2d p=%6.3f", mean(v), sum(v > 0), length(v), p)
     }
     cat(sprintf("  n=%-7s ranking %s | stage II %s | total %s\n",
                 if (nn == 0) "gate_n" else nn,
                 f(ss$rank_gain2), f(ss$stage2_gain), f(ss$total_gain)))
   }
   cat("  (ranking + stage II = total, by construction, on these cells)\n")
   cat("\n  per dataset (mean over the n grid):\n")
   agg <- aggregate(cbind(ranking = rank_gain2, stageII = stage2_gain,
                          total = total_gain) ~ data.name,
                    sp, mean, na.action = na.omit)
   agg[-1] <- round(agg[-1], 4)
   print(agg[order(-agg$stageII), ], row.names = FALSE)
   cat("\n  If stage II's contribution is ~0, say so: the honest reading is that\n",
       "  the method is carried by the gate and the dip ranking, and stage II\n",
       "  should be presented as optional rather than defended.\n", sep = "")
   }
 }
}

cat("\nwritten to: ", normalizePath(dir_out), "\n", sep = "")
