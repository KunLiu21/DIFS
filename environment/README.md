# Dependencies and execution record

The successful Kumar run used Ubuntu 22.04.4, R 4.4.1, Seurat 4.4.0,
SeuratObject 5.0.2, SC3 1.32.0, FEAST 1.12.0, Matrix 1.7-0,
diptest 0.77-2 and matrixStats 1.4.1. The complete recorded namespace inventory
is [validated-session-packages.csv](validated-session-packages.csv); the original
session report is under `example/validated/kumar-20260925/sessionInfo.txt`.
These are observed versions from an existing environment, not a tested lockfile
or proof of a fresh installation. In particular, the data-download package was
not needed by this prepared-input run.

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

Run `Rscript scripts/check_environment.R > environment-check.txt` to record
installed versions. The full example also saves `sessionInfo()`. The validated
run used an existing author environment; clean-environment installation remains
untested. Use a trusted existing container or prepare a compatible R environment;
no private container path is sufficient to establish a portable installation.
