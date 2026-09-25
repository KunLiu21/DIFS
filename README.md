# DIFS: discriminative feature selection for single-cell clustering

DIFS combines a dip-based stage-I ranking with a preliminary-cluster-dependent
stage-II feature selection and mixing procedure. This revision provides the
complete two-stage entry point, a real-data worked example, the benchmark
sources, and corrected supplementary results.

**Validation status:** the Kumar full example is prepared but is not yet
certified as executed in a clean environment. See [VALIDATION.md](VALIDATION.md)
for the checks actually completed. The old synthetic stage-I example is a quick
diagnostic, not evidence of full-pipeline reproducibility.

## Complete worked example: Kumar

From the repository root, in an R environment with the dependencies listed in
[environment/README.md](environment/README.md):

```sh
Rscript example/prepare_kumar.R
Rscript example/run_kumar.R
```

Alternatively, use the Kumar Seurat object already prepared for the benchmark:

```sh
Rscript example/run_kumar.R input=/your/data/Kumar.rds out=example/output-local
```

Kumar contains 246 cells and three annotated classes after benchmark preparation.
The example requests 100 features with seed 1 and **true-k=3**. This is the
paper's oracle configuration: k uses the number of annotated classes. The
individual labels are kept out of feature selection and used for evaluation.
It is a worked example, not a new comparison demonstrating superiority.

The runner executes stage-I scoring, the intermediate feature-count search,
repaired SC3 preliminary clustering, stage-II scoring, mixing-ratio search and
refined Louvain final clustering. Required dependencies are never silently
skipped. It saves the selected genes, preliminary/final labels, metrics, input
checksums, configuration, R session and a log. See [example/README.md](example/README.md).

## Use the complete method

```r
source("R/difs.R")
# counts: prepared raw counts, unique genes in rows and cells in columns
# k: supplied number of clusters; state how it was obtained
fit <- difs_fit(counts, k=3, n_features=100, seed=1)
fit$features
fit$labels
```

The API loads the scientific function definitions directly from `benchmark/`;
there is no second independently edited implementation of the selection logic.
The gate used by the revision is **submitted**: expression above ln(5) in at
least max(0.05*N,35) cells after log(1+10,000*count/library-size) normalization.
The bundled Gaussian null lookup is deterministic but numerically approximate.
Stage-I and stage-II scores do not provide general p-value/FDR guarantees.
Actual returned feature counts and the nominal mixing ratio are recorded;
they need not equal the requested count or realized stage membership.

The historical SC3 wrapper defaults to seed 1 and two cores. The API preserves
that behavior and records it; the outer seed is not claimed to control every
internal algorithm's random state. Environment differences may change results.

## Reproduce the paper's analyses

- [benchmark/README.md](benchmark/README.md): preparation, configuration, runs,
  repair provenance, summary tables and figures.
- [results/README.md](results/README.md): corrected S0–S5 and diagnostic tables.
- [simulation/](simulation/): the simulation analysis sources.
- [provenance/source_manifest.json](provenance/source_manifest.json): source
  identities and packaging changes.
- [CHANGELOG.md](CHANGELOG.md): implementation and reporting corrections.

The main grid has **2,320 observed configurations out of 2,340 planned**, plus
117 ablation runs and 52 Hartigan-reference runs. It is not 2,340 successful
runs. Raw matrices and per-run RDS files are not redistributed here; download
instructions and the corrected tabular results are provided. A single worked
example is not a clean-environment rerun of all benchmarks.

To regenerate corrected figures without rerunning clustering:

```sh
Rscript scripts/rebuild_figures.R
```

Only the published comparator scope is used in the main figures. Historical
internal method arms remain in the data for traceability and are not new
required comparisons.

## Lightweight diagnostic

```sh
Rscript demo/run_demo.R
Rscript tests/test_lightweight.R
```

These require only R and `diptest`; they do not replace the full example.
The earlier implementation under `code/` is retained for history and is not the
revision entry point.

## Availability and citation

Code is licensed under MIT; dependencies and source datasets retain their own
terms. See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff). Cite the exact
commit/release you use. A revision DOI is not claimed until one is issued.
