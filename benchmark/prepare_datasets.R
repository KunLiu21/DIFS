
difs_parse_args <- function(args, env = globalenv()) {
  for (a in args) {
    kv <- regmatches(a, regexec("^([A-Za-z._][A-Za-z0-9._]*)=(.*)$", a))[[1]]
    if (length(kv) != 3L) { message("ignoring argument: ", a); next }
    key <- kv[2L]; val <- kv[3L]
    val <- sub('^"(.*)"$', "\\1", val)
    val <- sub("^'(.*)'$", "\\1", val)
    v <- if (grepl("^(TRUE|FALSE|T|F)$", val)) as.logical(val)
         else if (grepl("^-?[0-9]+(\\.[0-9]+)?([eE][-+]?[0-9]+)?$", val)) as.numeric(val)
         else val
    assign(key, v, envir = env)
  }
  invisible(TRUE)
}
difs_parse_args(commandArgs(TRUE))

if (!exists("out_dir", inherits = FALSE)) out_dir <- "../source"
if (!exists("only",    inherits = FALSE)) only    <- NULL
if (!exists("cache",   inherits = FALSE)) cache   <- file.path(path.expand("~"), ".cache", "DIFS", "ExperimentHub")
if (!exists("hemberg_dir", inherits = FALSE)) hemberg_dir <- "../source/hemberg_raw"
if (!exists("norm",    inherits = FALSE)) norm    <- "lognorm"   # or "sct"
stopifnot(norm %in% c("lognorm", "sct"))

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

cat("libPaths:\n"); for (l in .libPaths()) cat("  ", l, "\n")
if (!nzchar(Sys.getenv("DIFS_LIB")))
  cat("NOTE: DIFS_LIB is not set.  If packages are reported missing, run\n",
      "      set DIFS_LIB only if a separate R library is required\n",
      "      set DIFS_EXCLUDE_LIB only to omit an incompatible library\n", sep = "")
for (p in c("DuoClustering2018", "scRNAseq", "Seurat")) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA)
  w <- tryCatch(dirname(find.package(p)),               error = function(e) "-")
  cat(sprintf("  %-20s %-10s %s\n", p, if (is.na(v)) "MISSING" else v, w))
}
cat("\n")

Sys.setenv(EXPERIMENT_HUB_CACHE = cache, ANNOTATION_HUB_CACHE = cache,
           BFC_CACHE = cache)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(SummarizedExperiment)
})

pick_labels <- function(sce, prefer = c("trueclass", "phenoid",
                                        "cell.type", "cell_type",
                                        "celltype", "cell type",
                                        "publishedClusters",
                                        "cell_type1", "cell_type2",
                                        "cluster_id", "label", "Group",
                                        "level1class", "broad_type",
                                        "CLUSTER", "cluster", "free_annotation"),
                        force = NULL, guess = FALSE) {
  cd <- SummarizedExperiment::colData(sce)
  if (!is.null(force)) {
    if (!force %in% colnames(cd))
      stop("labels= asked for column '", force, "', which does not exist. ",
           "Run inspect= on this dataset to see the real column names.")
    return(list(col = force, y = cd[[force]]))
  }
  for (p in prefer) if (p %in% colnames(cd)) return(list(col = p, y = cd[[p]]))
  if (!isTRUE(guess)) {
    stop("no known label column. Run  inspect=<dataset>  to see what colData ",
         "holds, then pass  labels=<dataset>=<column>.  Columns present: ",
         paste(colnames(cd), collapse = ", "))
  }
  for (p in colnames(cd)) {
    v <- cd[[p]]
    if (is.factor(v) || is.character(v)) {
      k <- length(unique(v))
      if (k >= 2 && k <= 30) {
        warning("guessing label column '", p, "' -- VERIFY THIS", call. = FALSE)
        return(list(col = p, y = v))
      }
    }
  }
  stop("no label column found; colData has: ", paste(colnames(cd), collapse = ", "))
}

difs_inspect <- function(sce, nm) {
  cat("\n########## ", nm, " ##########\n", sep = "")
  cat("class  : ", class(sce)[1], "\n", sep = "")
  cat("dim    : ", nrow(sce), " genes x ", ncol(sce), " cells\n", sep = "")
  cat("assays : ", paste(SummarizedExperiment::assayNames(sce), collapse = ", "),
      "\n", sep = "")
  ae <- tryCatch(SingleCellExperiment::altExpNames(sce), error = function(e) NULL)
  if (length(ae)) cat("altExps: ", paste(ae, collapse = ", "), "\n", sep = "")
  cd <- SummarizedExperiment::colData(sce)
  cat("colData columns (", ncol(cd), "):\n", sep = "")
  for (p in colnames(cd)) {
    v <- cd[[p]]
    if (!(is.factor(v) || is.character(v) || is.logical(v) || is.integer(v))) {
      cat(sprintf("  %-28s %-10s (not a label candidate)\n", p, class(v)[1]))
      next
    }
    u <- unique(as.character(v))
    cat(sprintf("  %-28s %-10s %4d distinct | %s\n", p, class(v)[1], length(u),
                paste(utils::head(u, 6), collapse = " / ")))
  }
  invisible(NULL)
}

difs_pick_assay <- function(sce, prefer = c("counts", "normcounts", "exprs",
                                            "logcounts", "tpm", "fpkm", "rpkm")) {
  an <- SummarizedExperiment::assayNames(sce)
  for (p in prefer) if (p %in% an) return(p)
  if (!length(an)) stop("object has no assays")
  an[1L]
}

difs_scale_report <- function(cnt, binary.bound = log(5), name = "",
                              block = 2000L) {
  lib <- colSums(cnt); lib[lib == 0] <- 1
  nc <- ncol(cnt); hits <- 0
  for (a in seq(1L, nc, by = block)) {
    b  <- min(a + block - 1L, nc)
    z  <- log1p(sweep(cnt[, a:b, drop = FALSE], 2, lib[a:b], "/") * 1e4)
    hits <- hits + sum(z > binary.bound)
  }
  pct_gate  <- 100 * hits / (as.numeric(nrow(cnt)) * nc)
  tot       <- sum(lib)
  top_share <- 100 * max(rowSums(cnt)) / tot
  nz        <- cnt[cnt != 0]
  nz        <- nz[seq_len(min(10000L, length(nz)))]
  is_int    <- all(abs(nz - round(nz)) < 1e-8)
  message(sprintf(
    "  scale check: %.3f%% of entries > log(5) | top gene = %.3f%% of library | integer = %s",
    pct_gate, top_share, is_int))
  if (pct_gate < 0.05)
    stop("REJECTED: only ", sprintf("%.4f", pct_gate), "% of entries exceed ",
         "log(5) after CP10K normalisation, so the stage I gate would pass no ",
         "genes.  This matrix is almost certainly rank- or quantile-transformed ",
         "rather than counts (cf. FEAST::Yan).  Do not use it.")
  if (top_share < 0.1)
    warning("dataset ", name, ": the largest gene is only ",
            sprintf("%.3f", top_share), "% of the library -- unusually flat ",
            "for counts.  Verify the source before using it.", call. = FALSE)
  list(pct_gate = pct_gate, top_share = top_share, integer = is_int)
}

build_seurat <- function(sce, labels, min.cells = 3, min.features = 200,
                         norm = "lognorm", assay_name = NULL) {
  if (is.null(assay_name)) assay_name <- difs_pick_assay(sce)
  if (!identical(assay_name, "counts"))
    message("  NOTE: no raw counts; using assay '", assay_name,
            "'.  Record this in Table 1 -- the log(5) gate is defined on ",
            "library-size-normalised counts.")
  attr(labels, "assay_used") <- assay_name
  cnt <- SummarizedExperiment::assay(sce, assay_name)
  cnt <- as.matrix(cnt)
  storage.mode(cnt) <- "numeric"
  if (is.null(colnames(cnt))) colnames(cnt) <- paste0("cell", seq_len(ncol(cnt)))
  if (is.null(rownames(cnt)))
    stop("this matrix has no gene names.  Refusing to invent them: a dataset ",
         "with synthetic identifiers (gene1, gene2, ...) is silently unusable ",
         "for any comparison of gene sets across datasets.  The identifiers ",
         "are usually in the object's rowData -- fix the loader, not this.")
  keep <- !duplicated(rownames(cnt))
  cnt <- cnt[keep, , drop = FALSE]
  sc <- difs_scale_report(cnt, name = assay_name)

  keep_lab <- !is.na(labels) & nzchar(trimws(as.character(labels)))
  if (any(!keep_lab)) {
    message("  ", sum(!keep_lab), " of ", length(labels),
            " cells have no label; dropped before building the object")
    cnt <- cnt[, keep_lab, drop = FALSE]
    labels <- labels[keep_lab]
  }

  seu <- Seurat::CreateSeuratObject(cnt, min.cells = min.cells,
                                    min.features = min.features)
  lab <- labels[match(colnames(seu), colnames(cnt))]
  seu$trueclass <- as.character(lab)
  seu@misc$assay_used   <- assay_name
  seu@misc$scale_report <- sc
  seu <- if (norm == "sct") {
    Seurat::SCTransform(seu, verbose = FALSE)
  } else {
    Seurat::NormalizeData(seu, normalization.method = "LogNormalize",
                          scale.factor = 1e4, verbose = FALSE)
  }
  seu
}

save_one <- function(seu, name, out_dir) {
  f <- file.path(out_dir, paste0(name, ".rds"))
  saveRDS(seu, f)
  data.frame(dataset = name, status = "ok",
             assay = paste0(Seurat::DefaultAssay(seu), " <- ",
                            if (is.null(seu@misc$assay_used)) "counts" else
                              seu@misc$assay_used),
             cells = ncol(seu), genes = nrow(seu),
             classes = length(unique(seu$trueclass)),
             pct_gate  = if (is.null(seu@misc$scale_report)) NA else
                           round(seu@misc$scale_report$pct_gate, 3),
             top_share = if (is.null(seu@misc$scale_report)) NA else
                           round(seu@misc$scale_report$top_share, 3),
             file = f, stringsAsFactors = FALSE)
}

need_pkg <- function(pkg) {
  if (pkg %in% loadedNamespaces()) return(invisible(TRUE))
  inst <- tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) NA)
  if (is.na(inst)) stop("package ", pkg, " is not installed")
  tryCatch(loadNamespace(pkg), error = function(e)
    stop("package ", pkg, " ", inst, " is installed but fails to load: ",
         conditionMessage(e)))
  invisible(TRUE)
}

fail <- function(name, msg)
  data.frame(dataset = name, status = "FAILED", assay = NA, cells = NA, genes = NA,
             classes = NA, pct_gate = NA, top_share = NA,
             file = substr(msg, 1, 200), stringsAsFactors = FALSE)

duo <- c(Koh = "sce_full_Koh", KohTCC = "sce_full_KohTCC",
         Kumar = "sce_full_Kumar", KumarTCC = "sce_full_KumarTCC",
         Trapnell = "sce_full_Trapnell", TrapnellTCC = "sce_full_TrapnellTCC",
         SimKumar4easy = "sce_full_SimKumar4easy",
         SimKumar4hard = "sce_full_SimKumar4hard",
         Zhengmix4eq   = "sce_full_Zhengmix4eq",
         Zhengmix4uneq = "sce_full_Zhengmix4uneq",
         Zhengmix8eq   = "sce_full_Zhengmix8eq")

scrnaseq <- c(
  Darmanis     = "DarmanisBrainData",        #   466 cells, human brain, ~9 types
  Fletcher     = "FletcherOlfactoryData",    #   849 cells, mouse olfactory
  Baron        = "BaronPancreasData",        #  8569 cells, inDrop,     GSE84133
  Muraro       = "MuraroPancreasData",       #  2126 cells, CEL-Seq2,   GSE85241
  Segerstolpe  = "SegerstolpePancreasData",  #  3514 cells, Smart-seq2, E-MTAB-5061
  Lawlor       = "LawlorPancreasData",       #   638 cells, Fluidigm C1, GSE86469
  Zeisel       = "ZeiselBrainData",          #  3005 cells, STRT-seq, 7 / 47 types
  Tasic        = "TasicBrainData",           #  1809 cells, SMARTer,  17 / 49 types
  Shekhar      = "ShekharRetinaData",        # 27499 cells, Drop-seq, 19 types -- LARGE
  Pollen       = "PollenGliaData",           # SC3 gold standard,   301 cells
  Kolodziejczyk= "KolodziejczykESCData",     # SC3 gold standard,   704 cells
  Usoskin      = "UsoskinBrainData",         # SC3 silver standard, 622 cells
  Romanov      = "RomanovBrainData",         # FEAST,              2881 cells
  Nestorowa    = "NestorowaHSCData"          # FEAST, 1656 cells -- CONTINUUM,
)

HEMBERG_BASE <- "https://scrnaseq-public-datasets.s3.amazonaws.com/scater-objects/"
hemberg <- c(
  Yan           = "yan",            #  90 cells, human preimplantation embryo
  Biase         = "biase",          #  56 cells, mouse embryo
  Deng          = "deng-reads",     # 268 cells, mouse embryo (READS, not rpkm)
  Goolam        = "goolam",         # 124 cells, mouse embryo
  PollenH       = "pollen",         # human, SC3 gold standard
  KolodziejczykH= "kolodziejczyk"   # mouse ESC, SC3 gold standard
)

CONQUER_BASE <- "http://imlspenticton.uzh.ch/robinson_lab/conquer/data-mae/"
conquer <- c(
  DengC          = "GSE45719",   # Deng          SC3 gold standard / FEAST
  BiaseC         = "GSE57249",   # Biase         SC3 gold standard
  KolodziejczykC = "EMTAB2600"   # Kolodziejczyk SC3 gold standard
)

aliasreg <- c(Zeisel2 = "Zeisel", Tasic2 = "Tasic")
alias_label <- c(Zeisel2 = "level2class", Tasic2 = "broad_type")

difs_conquer_sce <- function(acc) {
  need_pkg("MultiAssayExperiment")
  dir.create(hemberg_dir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(hemberg_dir, paste0(acc, "_conquer.rds"))
  if (!file.exists(f) || file.size(f) < 1000) {
    u <- paste0(CONQUER_BASE, acc, ".rds")
    message("  downloading ", u)
    utils::download.file(u, f, mode = "wb", quiet = FALSE)
  } else message("  using cached ", f)
  mae <- readRDS(f)
  ex  <- MultiAssayExperiment::experiments(mae)
  nmx <- if ("gene" %in% names(ex)) "gene" else names(ex)[1]
  sce <- ex[[nmx]]
  an  <- SummarizedExperiment::assayNames(sce)
  if ("count" %in% an && !"counts" %in% an)
    SummarizedExperiment::assay(sce, "counts") <-
      round(SummarizedExperiment::assay(sce, "count"))
  cd <- as.data.frame(MultiAssayExperiment::colData(mae))
  cd <- cd[match(colnames(sce), rownames(cd)), , drop = FALSE]
  if ("characteristics_ch1" %in% colnames(cd)) {
    parts <- strsplit(as.character(cd$characteristics_ch1), "\\s*;\\s*")
    keys  <- unique(unlist(lapply(parts, function(z)
               sub(":.*$", "", z[grepl(":", z)]))))
    for (k in keys) {
      v <- vapply(parts, function(z) {
             h <- z[startsWith(z, paste0(k, ":"))]
             if (length(h)) trimws(sub("^[^:]*:", "", h[1])) else NA_character_
           }, character(1))
      cn <- make.names(paste0("ch1_", k))
      cd[[cn]] <- v
    }
  }
  SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(cd)
  sce
}

pkgdata <- c(YanF = "FEAST::Yan")

difs_data_sce <- function(spec) {
  z <- strsplit(spec, "::", fixed = TRUE)[[1]]
  pkg <- z[1]; dname <- z[2]
  need_pkg(pkg)
  e <- new.env()
  utils::data(list = dname, package = pkg, envir = e)
  objs <- mget(ls(e), envir = e)
  if (!length(objs)) stop("data(", dname, ", package=", pkg, ") loaded nothing")
  message("  data() produced: ", paste(names(objs), collapse = ", "))
  for (o in objs) if (methods::is(o, "SummarizedExperiment")) return(o)
  is2d <- vapply(objs, function(o)
            is.matrix(o) || is.data.frame(o) || methods::is(o, "Matrix"),
            logical(1))
  if (!any(is2d))
    stop("no matrix among: ", paste(names(objs), collapse = ", "))
  big <- objs[is2d]
  cnt <- big[[which.max(vapply(big, function(o) prod(dim(o)), numeric(1)))]]
  cnt <- as.matrix(cnt)
  if (is.null(colnames(cnt))) colnames(cnt) <- paste0("cell", seq_len(ncol(cnt)))
  cand <- Filter(function(o) (is.character(o) || is.factor(o)) &&
                   length(o) == ncol(cnt), objs)
  cd <- if (length(cand)) as.data.frame(lapply(cand, as.character),
                                        stringsAsFactors = FALSE) else
          data.frame(row.names = colnames(cnt))
  cd$from_colnames <- sub("[._-]*[0-9]+[._-]*$", "",
                          sub("\\.(RPKM|FPKM|TPM|counts?)$", "", colnames(cnt),
                              ignore.case = TRUE))
  rownames(cd) <- colnames(cnt)
  SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = cnt), colData = S4Vectors::DataFrame(cd))
}

eh_prefix <- c(
  Segerstolpe   = "Segerstolpe pancreas",   # EH2575 / EH2576 / EH2577
  Zeisel        = "Zeisel brain",           # EH2580 / EH2581 / EH2582
  Tasic         = "Tasic brain",            # EH2578 / EH2579 (no rowData)
  Kolodziejczyk = "Kolodziejczyk ESC",      # EH3107 only -- NO colData
  Nestorowa     = "Nestorowa HSC",
  Pollen        = "Pollen Glia",
  Usoskin       = "Usoskin brain")

difs_eh_records <- function(prefix) {
  need_pkg("ExperimentHub")
  eh <- ExperimentHub::ExperimentHub()
  ti <- eh$title
  hit <- which(grepl(tolower(prefix), tolower(ti), fixed = TRUE))
  data.frame(id = names(eh)[hit], title = ti[hit], stringsAsFactors = FALSE)
}

difs_eh_sce <- function(prefix) {
  r <- difs_eh_records(prefix)
  if (!nrow(r))
    stop("no ExperimentHub record whose title contains '", prefix, "'. ",
         "Run  eh=<dataset>  to list what the hub actually has.")
  message("  ExperimentHub records matching '", prefix, "':")
  for (i in seq_len(nrow(r))) message("    ", r$id[i], "  ", r$title[i])
  need_pkg("ExperimentHub")
  eh  <- ExperimentHub::ExperimentHub()
  get1 <- function(pat) {
    j <- grep(pat, r$title, ignore.case = TRUE)
    if (!length(j)) return(NULL)
    j <- j[which.min(nchar(r$title[j]))]
    message("  using ", r$id[j], "  (", r$title[j], ")")
    eh[[r$id[j]]]
  }
  cnt <- get1("counts?$")
  if (is.null(cnt)) cnt <- get1("counts")
  if (is.null(cnt)) stop("no counts record among: ",
                         paste(r$title, collapse = " | "))
  cd <- get1("colData|coldata")
  rd <- get1("rowData|rowdata")
  if (is.null(colnames(cnt)) && !is.null(cd)) colnames(cnt) <- rownames(cd)
  sce <- SingleCellExperiment::SingleCellExperiment(
           assays = list(counts = cnt))
  if (!is.null(cd)) {
    cd <- methods::as(cd, "DataFrame")
    if (nrow(cd) == ncol(sce)) SummarizedExperiment::colData(sce) <- cd else
      warning("colData has ", nrow(cd), " rows but the matrix has ",
              ncol(sce), " columns -- colData dropped", call. = FALSE)
  }
  if (!is.null(rd)) {
    rd <- methods::as(rd, "DataFrame")
    if (nrow(rd) == nrow(sce)) {
      SummarizedExperiment::rowData(sce) <- rd
      if (is.null(rownames(sce))) {
        cand <- grep("^(symbol|gene|feature|name|id)", tolower(colnames(rd)),
                     value = TRUE)
        if (!length(cand))
          cand <- colnames(rd)[vapply(as.list(rd), function(z)
                    is.character(z) || is.factor(z), logical(1))]
        if (length(cand)) {
          nmv <- as.character(rd[[cand[1]]])
          message("  rownames taken from rowData$", cand[1],
                  "  (e.g. ", paste(utils::head(nmv, 3), collapse = ", "), ")")
          rownames(sce) <- make.unique(nmv)
        } else {
          warning("no usable identifier column in rowData; columns are: ",
                  paste(colnames(rd), collapse = ", "), call. = FALSE)
        }
      }
    } else {
      warning("rowData has ", nrow(rd), " rows but the matrix has ", nrow(sce),
              " -- rowData dropped", call. = FALSE)
    }
  }
  sce
}

difs_fetch_sce <- function(nm) {
  if (nm %in% names(aliasreg)) nm <- aliasreg[[nm]]   # same source, finer labels
  if (nm %in% names(duo)) {
    need_pkg("DuoClustering2018")
    return(getExportedValue("DuoClustering2018", duo[[nm]])())
  }
  if (nm %in% names(scrnaseq)) {
    need_pkg("scRNAseq")
    fn <- getExportedValue("scRNAseq", scrnaseq[[nm]])
    return(tryCatch(fn(), error = function(e1) {
      message("  gypsum backend failed (", conditionMessage(e1), ")")
      message("  try 2: same loader with legacy = TRUE ...")
      tryCatch(fn(legacy = TRUE), error = function(e2) {
        message("  legacy path failed too (", conditionMessage(e2), ")")
        pf <- eh_prefix[[nm]]
        if (is.null(pf))
          stop("both scRNAseq paths failed and no ExperimentHub title prefix ",
               "is registered for ", nm, ".  Add one to eh_prefix after ",
               "running  eh=", nm, ".  Last error: ", conditionMessage(e2))
        message("  try 3: assembling from ExperimentHub directly, no altExp ...")
        difs_eh_sce(pf)
      })
    }))
  }
  if (nm %in% names(pkgdata)) return(difs_data_sce(pkgdata[[nm]]))
  if (nm %in% names(conquer)) return(difs_conquer_sce(conquer[[nm]]))
  if (nm %in% names(hemberg)) {
    f <- file.path(hemberg_dir, paste0(hemberg[[nm]], ".rds"))
    if (!file.exists(f))
      stop("the Hemberg S3 bucket returns HTTP 403 and is no longer public. ",
           "Put the object at ", f, " by hand, or use the conquer entry ",
           "(DengC / BiaseC / KolodziejczykC) instead.")
    return(readRDS(f))
  }
  stop("unknown dataset: ", nm)
}

wanted <- c(names(duo), names(scrnaseq), names(pkgdata),
            names(conquer), names(hemberg), names(aliasreg))
if (!is.null(only)) wanted <- intersect(wanted, strsplit(only, ",")[[1]])

label_override <- list()
if (exists("labels", inherits = FALSE) && is.character(labels)) {
  for (kv in strsplit(labels, ",")[[1]]) {
    z <- strsplit(kv, "=", fixed = TRUE)[[1]]
    if (length(z) == 2L) label_override[[trimws(z[1])]] <- trimws(z[2])
  }
}

for (a in names(alias_label))
  if (is.null(label_override[[a]])) label_override[[a]] <- alias_label[[a]]

if (exists("eh", inherits = FALSE) && is.character(eh)) {
  for (nm in strsplit(eh, ",")[[1]]) {
    nm <- trimws(nm)
    pf <- if (nm %in% names(eh_prefix)) eh_prefix[[nm]] else nm
    cat("\n########## ", nm, "  (prefix '", pf, "') ##########\n", sep = "")
    r <- tryCatch(difs_eh_records(pf),
                  error = function(e) { cat("  FAILED: ", conditionMessage(e),
                                            "\n", sep = ""); NULL })
    if (is.null(r)) next
    if (!nrow(r)) cat("  no matching record\n") else
      for (i in seq_len(nrow(r))) cat(sprintf("  %-12s %s\n", r$id[i], r$title[i]))
  }
  cat("\neh mode: nothing built.\n")
  quit(save = "no")
}

if (exists("inspect", inherits = FALSE) && is.character(inspect)) {
  for (nm in strsplit(inspect, ",")[[1]]) {
    nm <- trimws(nm)
    message("=== fetching ", nm, " for inspection ===")
    tryCatch(difs_inspect(difs_fetch_sce(nm), nm),
             error = function(e) cat("\n########## ", nm,
                                     " ##########\n  FAILED: ",
                                     conditionMessage(e), "\n", sep = ""))
  }
  cat("\ninspect mode: nothing built. Choose a column, then rerun with\n",
      "  only=<Name> labels=<Name>=<Column>\n", sep = "")
  quit(save = "no")
}

if (isTRUE(get0("list", ifnotfound = FALSE))) {
  for (pk in c("DuoClustering2018", "scRNAseq")) {
    cat("\n==== ", pk, " ====\n", sep = "")
    ok <- tryCatch({ loadNamespace(pk); TRUE }, error = function(e) {
      cat("  not loadable: ", conditionMessage(e), "\n", sep = ""); FALSE })
    if (!ok) next
    ex0 <- sort(getNamespaceExports(pk)); ex <- ex0
    ex <- if (pk == "DuoClustering2018") grep("^sce_full_", ex, value = TRUE) else
            grep("Data$", ex, value = TRUE)
    cat(paste0("  ", ex, collapse = "\n"), "\n")
    if (pk == "scRNAseq" && "listDatasets" %in% ex0) {
      cat("\n  --- listDatasets() ---\n")
      print(utils::head(as.data.frame(
        getExportedValue(pk, "listDatasets")()), 100))
    }
    known <- if (pk == "DuoClustering2018") duo else scrnaseq
    miss <- setdiff(known, ex)
    if (length(miss)) cat("  !! referenced here but NOT exported: ",
                          paste(miss, collapse = ", "), "\n", sep = "")
  }
  for (pk in c("FEAST", "SC3", "TSCAN", "monocle")) {
    d <- tryCatch(utils::data(package = pk)$results, error = function(e) NULL)
    if (is.null(d) || !nrow(d)) next
    cat("\n==== data(package = \"", pk, "\") ====\n", sep = "")
    print(d[, c("Item", "Title"), drop = FALSE])
  }
  cat("\nlist mode: nothing built.\n")
  quit(save = "no")
}

res <- list()

for (nm in intersect(wanted, names(duo))) {
  message("=== ", nm, " (DuoClustering2018) ===")
  res[[nm]] <- tryCatch({
    need_pkg("DuoClustering2018")
    sce <- difs_fetch_sce(nm)
    lb  <- pick_labels(sce, force = label_override[[nm]])
    message("  labels from colData$", lb$col, "  (", length(unique(lb$y)), " classes)")
    save_one(build_seurat(sce, lb$y, norm = norm), nm, out_dir)
  }, error = function(e) { message("  FAILED: ", conditionMessage(e))
                           fail(nm, conditionMessage(e)) })
}

for (nm in intersect(wanted, names(scrnaseq))) {
  message("=== ", nm, " (scRNAseq) ===")
  res[[nm]] <- tryCatch({
    sce <- difs_fetch_sce(nm)
    lb  <- pick_labels(sce, force = label_override[[nm]])
    message("  labels from colData$", lb$col, "  (", length(unique(lb$y)), " classes)")
    save_one(build_seurat(sce, lb$y, norm = norm), nm, out_dir)
  }, error = function(e) { message("  FAILED: ", conditionMessage(e))
                           fail(nm, conditionMessage(e)) })
}

for (nm in intersect(wanted, names(aliasreg))) {
  message("=== ", nm, " (alias of ", aliasreg[[nm]], ", labels = ",
          label_override[[nm]], ") ===")
  res[[nm]] <- tryCatch({
    sce <- difs_fetch_sce(nm)
    lb  <- pick_labels(sce, force = label_override[[nm]])
    message("  labels from colData$", lb$col, "  (", length(unique(lb$y)), " classes)")
    save_one(build_seurat(sce, lb$y, norm = norm), nm, out_dir)
  }, error = function(e) { message("  FAILED: ", conditionMessage(e))
                           fail(nm, conditionMessage(e)) })
}

for (nm in intersect(wanted, names(pkgdata))) {
  message("=== ", nm, " (", pkgdata[[nm]], ") ===")
  res[[nm]] <- tryCatch({
    sce <- difs_fetch_sce(nm)
    lb  <- pick_labels(sce, force = label_override[[nm]])
    message("  labels from colData$", lb$col, "  (", length(unique(lb$y)), " classes)")
    save_one(build_seurat(sce, lb$y, norm = norm), nm, out_dir)
  }, error = function(e) { message("  FAILED: ", conditionMessage(e))
                           fail(nm, conditionMessage(e)) })
}

for (nm in intersect(wanted, names(conquer))) {
  message("=== ", nm, " (conquer ", conquer[[nm]], ") ===")
  res[[nm]] <- tryCatch({
    sce <- difs_fetch_sce(nm)
    lb  <- pick_labels(sce, force = label_override[[nm]])
    message("  labels from colData$", lb$col, "  (", length(unique(lb$y)), " classes)")
    save_one(build_seurat(sce, lb$y, norm = norm), nm, out_dir)
  }, error = function(e) { message("  FAILED: ", conditionMessage(e))
                           fail(nm, conditionMessage(e)) })
}

for (nm in intersect(wanted, names(hemberg))) {
  message("=== ", nm, " (Hemberg lab S3) ===")
  res[[nm]] <- tryCatch({
    sce <- difs_fetch_sce(nm)
    if (!methods::is(sce, "SummarizedExperiment"))
      stop("legacy ", class(sce)[1], " object: install the matching scater ",
           "release or convert it once with scater::toSingleCellExperiment()")
    message("  assays: ", paste(SummarizedExperiment::assayNames(sce),
                                collapse = ", "))
    lb <- pick_labels(sce, force = label_override[[nm]])
    message("  labels from colData$", lb$col, "  (",
            length(unique(lb$y)), " classes)")
    save_one(build_seurat(sce, lb$y, norm = norm), nm, out_dir)
  }, error = function(e) { message("  FAILED: ", conditionMessage(e))
                           fail(nm, conditionMessage(e)) })
}

out <- do.call(rbind, res)
cat("\n================ summary ================\n")
print(out, row.names = FALSE)
utils::write.csv(out, file.path(out_dir, "dataset_build_report.csv"),
                 row.names = FALSE)

cat("\nNot attempted here:\n")
cat("  Hemberg lab S3 (scrnaseq-public-datasets) returns HTTP 403 since\n")
cat("             2026-09-08 and is no longer public, so Yan / Goolam have no\n")
cat("             working automatic source.  Deng, Biase and Kolodziejczyk are\n")
cat("             available through conquer instead (DengC/BiaseC/KolodziejczykC).\n")
cat("  Lung_Human -- provenance not recorded in the manuscript (Reviewer 2,\n")
cat("             minor 3).  Either recover the accession or drop it and say so.\n")
cat("\nnormalisation used: ", norm,
    if (norm == "lognorm") "  (matches the original objects)" else
      "  (NOTE: seurat_k_cluster() expects RNA_snn_res.*)", "\n", sep = "")
cat("files written to: ", normalizePath(out_dir), "\n")
