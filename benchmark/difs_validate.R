
DIFS_KEY <- c("data.name", "fmethod", "clust", "n", "k_policy", "gate_rule", "seed")

difs_key_string <- function(df, key = DIFS_KEY)
  do.call(paste, c(lapply(key, function(k) as.character(df[[k]])), sep = " | "))

expected_from_scenarios <- function(set) {
  f <- sprintf("../scenarios/scenarios_%s.Rdata", set)
  if (!file.exists(f)) stop("scenario table not found: ", f)
  e <- new.env(parent = emptyenv()); load(f, envir = e)
  s <- get("scenarios", envir = e)
  data.frame(
    data.name = gsub(".*/([^.]*).*", "\\1", as.character(s$data.path)),
    fmethod   = as.character(s$feature.selection.method),
    clust     = as.character(s$method),
    n         = as.integer(s$n_features),
    k_policy  = as.character(s$k_policy),
    gate_rule = if ("gate_rule" %in% names(s)) as.character(s$gate_rule) else "submitted",
    seed      = as.integer(s$seed),
    stringsAsFactors = FALSE)
}

difs_gate <- function(new, expected, old = NULL, require_full_labels = FALSE,
                      pair_key = DIFS_KEY,
                      required_new = c("ARI", "nfeat"),
                      required_old = c("ARI"),
                      expect_total = NULL) {
  fail <- character(0)
  note <- function(...) fail <<- c(fail, paste0(...))
  line <- function(ok, msg) cat(if (ok) "  [PASS] " else "  [FAIL] ", msg, "\n", sep = "")

  cat(strrep("=", 76), "\nGATE  completeness and identity\n", strrep("=", 76), "\n", sep = "")

  kn <- difs_key_string(new); ke <- difs_key_string(expected)
  dup     <- kn[duplicated(kn)]
  missing <- setdiff(ke, kn)
  extra   <- setdiff(kn, ke)
  ok <- length(dup) == 0L
  line(ok, sprintf("no duplicate configurations in the new results (%d duplicate key(s))", length(dup)))
  if (!ok) { note(length(dup), " duplicated configuration(s)"); print(utils::head(dup, 10)) }
  ok <- length(missing) == 0L
  line(ok, sprintf("every expected configuration present (%d of %d missing)", length(missing), length(ke)))
  if (!ok) { note(length(missing), " missing configuration(s)"); print(utils::head(missing, 15)) }
  ok <- length(extra) == 0L
  line(ok, sprintf("no unexpected configurations (%d extra)", length(extra)))
  if (!ok) { note(length(extra), " unexpected configuration(s)"); print(utils::head(extra, 10)) }

  chk_cols <- function(df, cols, what) {
    for (v in cols) {
      if (!v %in% names(df)) {
        line(FALSE, sprintf("%s: column '%s' is required and is absent", what, v))
        note(what, ": required column '", v, "' absent"); next
      }
      bad <- if (is.numeric(df[[v]])) !is.finite(df[[v]]) else is.na(df[[v]])
      ok <- !any(bad)
      line(ok, sprintf("%s: %s present and usable on every run (%d bad)", what, v, sum(bad)))
      if (!ok) note(what, ": ", v, " missing or non-finite on ", sum(bad), " run(s)")
    }
  }
  chk_cols(new, required_new, "new")
  if (!is.null(old)) chk_cols(old, required_old, "stored")

  if (require_full_labels) {
    na_drop <- sum(is.na(new$dropped))
    ok <- na_drop == 0L
    line(ok, sprintf("n_cells_dropped recorded for every run (%d missing -> cannot be verified)", na_drop))
    if (!ok) note(na_drop, " run(s) do not report n_cells_dropped")
    if (!is.null(expect_total)) {
      wrong <- sum(is.na(new$total) | new$total != expect_total)
      ok <- wrong == 0L
      line(ok, sprintf("every run covers all %d cells (%d run(s) report a different total)",
                       expect_total, wrong))
      if (!ok) { note(wrong, " run(s) do not cover all ", expect_total, " cells")
                 print(unique(new$total)) }
    } else {
      line(FALSE, "expect_total was not supplied -- cannot verify cell coverage")
      note("expect_total not supplied while require_full_labels = TRUE")
    }
    pos <- sum(!is.na(new$dropped) & new$dropped > 0)
    ok <- pos == 0L
    line(ok, sprintf("every run scored every cell (%d run(s) still drop cells)", pos))
    if (!ok) { note(pos, " run(s) still drop cells")
               print(new[!is.na(new$dropped) & new$dropped > 0,
                         intersect(c("data.name","fmethod","clust","n","total","dropped"), names(new))],
                     row.names = FALSE) }
  }

  if (!is.null(old)) {
    kn2 <- difs_key_string(new, pair_key)
    ko  <- difs_key_string(old, pair_key)
    dupo <- ko[duplicated(ko)]
    ok <- length(dupo) == 0L
    line(ok, sprintf("no duplicate configurations in the stored results (%d)", length(dupo)))
    if (!ok) { note(length(dupo), " duplicated stored configuration(s)"); print(utils::head(dupo, 10)) }
    unmatched <- setdiff(kn2, ko)
    ok <- length(unmatched) == 0L
    line(ok, sprintf("every new run has exactly one stored counterpart (%d unmatched)", length(unmatched)))
    if (!ok) { note(length(unmatched), " new run(s) with no stored counterpart"); print(utils::head(unmatched, 15)) }
  }

  cat(strrep("-", 76), "\n", sep = "")
  if (length(fail)) {
    cat("VERDICT: INVALID\n\n")
    for (f in fail) cat("  - ", f, "\n", sep = "")
    cat("\nNo decision is printed from an incomplete or ambiguous result set.\n",
        "Fix the listed problems -- resubmit the missing array indices, remove\n",
        "duplicates, or correct the scenario table -- and run this again.\n", sep = "")
    quit(save = "no", status = 2)
  }
  cat("VERDICT: complete and unambiguous -- proceeding\n\n")
  invisible(TRUE)
}
