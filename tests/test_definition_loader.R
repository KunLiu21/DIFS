# Import only requested definitions; never execute a script's other statements.
source("R/difs.R")
local({
  fixture <- tempfile(fileext=".R")
  on.exit(unlink(fixture))
  writeLines(c(
    'base::stop("qualified calls must not run")',
    '(function() stop("computed calls must not run"))()',
    'unwanted <- stop("unrequested assignments must not run")',
    'wanted <- function(x) x + 1',
    'constant = 7L'
  ),fixture)
  e <- new.env(parent=baseenv())
  difs_load_definitions(fixture,c("wanted","constant"),e)
  stopifnot(e$wanted(2)==3,e$constant==7L,
            setequal(ls(e),c("wanted","constant")))
  must_fail <- function(expr,pattern) {
    failure <- tryCatch({force(expr);NULL},error=conditionMessage)
    stopifnot(!is.null(failure),grepl(pattern,failure,fixed=TRUE))
  }
  must_fail(difs_load_definitions(fixture,"absent",new.env()),"Missing definitions")
  writeLines(c('wanted <- 1','wanted <- 2'),fixture)
  must_fail(difs_load_definitions(fixture,"wanted",new.env()),"Duplicate definition")
  # Exercise the exact public data-preparation source that exposed the bug.
  helpers <- new.env(parent=globalenv())
  wanted <- c("difs_pick_assay","difs_scale_report","build_seurat")
  difs_load_definitions("benchmark/prepare_datasets.R",wanted,helpers)
  stopifnot(setequal(ls(helpers),wanted),
            all(vapply(wanted,function(n) is.function(helpers[[n]]),logical(1))))
})
cat("PASS: qualified/computed calls skipped, definitions imported, missing/duplicate definitions rejected.\n")
