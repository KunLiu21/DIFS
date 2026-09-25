
difs_signrank <- function(d, y = NULL, max_m = 20L) {
  if (!is.null(y)) {
    if (length(d) != length(y)) stop("difs_signrank: unequal lengths")
    d <- d - y
  }
  if (!is.numeric(d) || any(!is.finite(d)))
    stop("difs_signrank: all paired differences must be finite numeric values")
  n_all <- length(d)
  d <- d[d != 0]
  m <- length(d)
  out <- list(n = n_all, n_zero = n_all - m, m = m, statistic = NA_real_,
              p.value = NA_real_, method = NA_character_,
              n_pos = sum(d > 0), n_neg = sum(d < 0))
  if (m == 0L) { out$p.value <- 1; out$method <- "all differences zero"; return(out) }
  if (m > max_m) stop("difs_signrank: ", m, " nonzero pairs exceeds max_m = ", max_m,
                      "; this project never has that many, so something is wrong")
  r  <- rank(abs(d))                     # midranks: ties are handled, not avoided
  Wp <- sum(r[d > 0])
  mu <- sum(r) / 2                       # E[W+] under sign symmetry
  signs <- as.matrix(expand.grid(rep(list(c(0, 1)), m)))
  null  <- as.vector(signs %*% r)
  p <- mean(abs(null - mu) >= abs(Wp - mu) - 1e-12)
  out$statistic <- Wp
  out$p.value   <- p
  out$method    <- sprintf("exact sign-permutation signed-rank, m = %d (%d zero%s dropped)",
                           m, n_all - m, if (n_all - m == 1L) "" else "s")
  out
}

difs_p <- function(d, y = NULL) difs_signrank(d, y)$p.value

difs_ci <- function(d, y = NULL) {
  if (!is.null(y)) d <- d - y
  if (!is.numeric(d) || any(!is.finite(d)))
    stop("difs_ci: all paired differences must be finite numeric values")
  if (length(d) < 2L) return(c(NA_real_, NA_real_))
  as.numeric(stats::t.test(d)$conf.int)
}
