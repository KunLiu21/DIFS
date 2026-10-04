# DIFS: discriminative feature selection for single-cell clustering

DIFS combines a dip-based stage-I ranking with a preliminary-cluster-dependent
stage-II feature selection and mixing procedure. This revision provides the
complete two-stage entry point, a real-data worked example, the benchmark
sources, and corrected supplementary results.

**Validation status:** the complete Kumar example passed both in the author's
existing environment (2026-09-25) and in a fresh Ubuntu 22.04 / R 4.4.1 installation
with a public-data download (verified 2026-10-04). The fresh run reproduced the
100 selected genes, their order, both partitions and all four evaluation metrics.
See [fresh-run evidence](example/validated/kumar-clean-20261004),
[earlier run](example/validated/kumar-20260925) and [VALIDATION.md](VALIDATION.md).
This is one tested configuration; it is not a full-grid or all-platform guarantee.

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
# counts: prepared counts/count estimates, genes in rows and cells in columns
# k: supplied number of clusters; state how it was obtained
fit <- difs_fit(counts, k=3, n_features=100, seed=1)
fit$features
fit$labels
```

Inputs must be finite and nonnegative with positive cell library sizes.
Fractional count estimates are retained without rounding, matching benchmark
preparation. This does not establish the scale of an arbitrary matrix: supply
the count-scale input, not already log-normalized expression or scaled residuals.

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

**Second-round materials (pre-release, 2026-10-04):**
[analysis guide](benchmark/REVISION2.md) and
[organized result snapshot](results/revision2/README.md) add component controls,
repeated runs, k sensitivity, markers, runtime and full-method simulations.
The component and marker summaries include the 12-decimal signed-rank zero/tie
correction, independently checked from their observations. This material update
has not yet been tagged as `v2.1-revision`.
The validated Kumar API and historical v2.0 result archive are unchanged.

- [benchmark/README.md](benchmark/README.md): preparation, configuration, runs,
  repair provenance, summary tables and figures.
- [results/README.md](results/README.md): corrected S0–S5 and diagnostic tables.
- [simulation/](simulation/): the simulation analysis sources.
- [provenance/source_manifest.json](provenance/source_manifest.json): source
  identities and packaging changes.
- [CHANGELOG.md](CHANGELOG.md): implementation and reporting corrections.

The historical v2.0 main grid has **2,320 observed configurations out of 2,340 planned**, plus
117 ablation runs and 52 Hartigan-reference runs. It is not 2,340 successful
runs. Raw matrices and per-run RDS files are not redistributed here; download
instructions and the corrected tabular results are provided. A single worked
example is not a clean-environment rerun of all benchmarks.

The second-round manuscript excludes the historical Variance control and uses
1852/1872 main-grid configurations. Its workbook also includes S0b/S5b and
S6–S9b; see the snapshot guide for exact filtering and missing/skipped counts.

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

These use R, `diptest` and `Matrix`; they do not replace the full example.
The earlier implementation under `code/` is retained for history and is not the
revision entry point.

## Availability and citation

Code is licensed under MIT; dependencies and source datasets retain their own
terms. See [LICENSE](LICENSE) and [CITATION.cff](CITATION.cff). Cite the exact
commit/release you use. A revision DOI is not claimed until one is issued.
