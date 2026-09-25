# Dependencies and execution record

The original cluster ran R 4.4 with Bioconductor 3.19-compatible libraries.
The dependency list below is a setup recipe, not a verified lockfile. An exact
reproduction claim requires recording the versions in the environment that
actually runs the full example. Do not mix incompatible Bioconductor releases.

Full example: Seurat, SeuratObject, diptest, SC3, SingleCellExperiment,
SummarizedExperiment, S4Vectors, matrixStats, FEAST, magrittr, and
DuoClustering2018 for downloading Kumar. FEAST supplies the MSE criterion used
by the benchmark even when no comparator arm is run.

From an R 4.4 environment, a candidate installation recipe is:

```r
install.packages("BiocManager")
BiocManager::install(version="3.19", ask=FALSE)
BiocManager::install(c("SC3", "FEAST", "SingleCellExperiment",
                      "DuoClustering2018"), ask=FALSE, update=FALSE)
install.packages(c("Seurat", "diptest", "matrixStats", "magrittr"))
```

Use a coherent environment supported by the packages; the latest CRAN package
set may differ from the original cluster. SeuratObject v5 access is adapted by
the example API; use the original compatible environment when testing the
historical benchmark scripts, which retain older API calls.

The full benchmark's legacy imports additionally require monocle, dplyr,
scRNAseq and comparator-specific dependencies. Its optional scran/glmGamPoi
handling is retained. The API avoids loading unused comparator code at startup.

Run `Rscript scripts/check_environment.R` to save installed versions. The full
example also saves `sessionInfo()`. No claim of clean-environment execution is
made until its actual output and log have been checked. Do not download an
unknown container: an author's existing R environment can run the prepared
example, and its exact environment must then be documented.
