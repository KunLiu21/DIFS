# Fresh installation and data reproduction

This recipe is under validation on branch `codex/clean-kumar-reproduction`.
A workflow starting is not a successful reproduction. See its actual logs and
artifacts before citing an outcome; the released cluster example remains a
separate, already-verified result.

The workflow `.github/workflows/clean-kumar.yml` uses a new Ubuntu 22.04 runner,
R 4.4.1, an empty user package library and empty ExperimentHub/AnnotationHub
caches. It does not restore a package cache, mount an author library, or consume
the author's prepared Kumar input. R's freshly installed base/recommended
library is available. Every required non-base dependency is installed from
public repositories.

CRAN dependencies initially come from Posit's 2024-09-30 Jammy snapshot; key
versions are then pinned with CRAN archives (Matrix 1.7-0, SeuratObject 5.0.2,
Seurat 4.4.0, diptest 0.77-2, matrixStats 1.4.1 and igraph 2.0.3).
Bioconductor 3.19 supplies SC3 1.32.0, FEAST 1.12.0 and the public data download
packages. Installed versions and OS packages are recorded. This is a recipe,
not a claim of a byte-identical reconstruction of the author's whole system.

The steps run `environment/install_clean.R`, then the existing documented
`example/prepare_kumar.R` with explicitly empty data caches, then
`example/run_kumar.R`, the saved-output verifier and the old/fresh comparison.
Success requires four metric agreement, complete labels, the same selected
gene set and equivalent preliminary/final cell partitions. An RDS checksum
alone is not required to be identical across serializations/environments.

The uploaded evidence excludes the downloaded expression matrix and R library.
It includes dependency installation logs, download/preparation provenance,
input checksums, selected genes, cell labels, full output and comparison logs.
Failure at any step fails the workflow; partial artifacts remain available.
