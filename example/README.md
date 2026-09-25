# Full two-stage worked example

The prepared-input example was executed and verified on 2026-09-25. See
[validated/kumar-20260925](validated/kumar-20260925) for its actual outputs,
logs and configuration. To recheck those saved outputs without rerunning the
method, run from the repository root:

```sh
Rscript scripts/verify_kumar_output.R example/validated/kumar-20260925
```

The data download/preparation entry point below was not rerun as part of that
validation; the run used the author's existing prepared Kumar object.

The input is Kumar from `DuoClustering2018::sce_full_Kumar()` (GSE60749), with
`phenoid` annotations. The preparation script uses the same helper function as
the packaged benchmark: minimum 3 cells per gene, minimum 200 detected genes per
cell, removal of duplicate gene names, and log-normalization with scale 10,000.
The processed object must have 246 cells and three classes. The package version,
label field and prepared input checksum are written to `example/input/`.

Run from the repository root. Install the listed dependencies, then execute
`Rscript example/prepare_kumar.R` and `Rscript example/run_kumar.R`.
For a prepared benchmark object use `input=/path/to/Kumar.rds`.
Choose a new `out=` directory for each execution; existing outputs are not
silently overwritten. A prepared object must retain counts and `trueclass`.

The output comprises:

- `input_diagnostics.csv`, `sessionInfo_before_fit.txt`: input validity,
  fractional-value counts and the environment captured before feature selection;
- `selected_genes.csv`, `cell_labels.csv`: explicit gene/cell identities;
- `summary.csv`: four metrics, stage sizes, gate size and actual feature count;
- `difs_fit.rds`: method output including the repaired SC3 audit;
- `sessionInfo.txt`, `input_checksums.csv`, `run.log`: execution provenance;
- `historical_comparison.csv`: differences from the recorded benchmark row.

`reference/benchmark_Kumar_n100_true_seed1.csv` is an existing benchmark result,
not an output invented for this example and not proof that the example has run.
Its ARI is approximately 0.988771469. In the historical environment, numerical
reproduction requires agreement of metrics (target absolute tolerance 1e-8),
feature count and cell coverage. A larger difference should be investigated
against input checksums, package versions, thread settings and random behavior;
it is not a reason to tune the example until it matches.

Execution success requires both stages to finish and all cells to have labels.
Numerical reproduction is a separate check. The reference does not contain the
historical selected genes, so identical gene membership cannot be verified from
that record. Newly generated genes and labels will be available for later tests.

True-k=3 is an explicitly annotation-informed example. The `difs_fit` API can
accept a k estimated independently from expression, but that is a different
configuration and should not be compared against this oracle reference row.
