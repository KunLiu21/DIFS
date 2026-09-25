# Benchmark reproduction

The code here comes from the corrected revision working version v0.5. Long
historical discussion comments were removed; scientific function definitions
are preserved. The repaired SC3 wrapper is loaded by the revision runner.
The public API reuses the same function definitions without executing unrelated
top-level package imports. `provenance/source_manifest.json` records the inputs.

## New benchmark runs

Use a compatible environment (environment/README.md). Run from `benchmark/`:

```sh
Rscript prepare_datasets.R only=Kumar
Rscript scenario_generator_ncurve.R
Rscript Method_comparison_ncurve.R k=1 set=ncurve
```

This generates a scenario grid and runs ONE scenario by its row number; inspect
the generated `../scenarios/scenarios_ncurve.Rdata` to choose the desired row.
Prepare all specified datasets before scheduling the entire grid. Scenario
generators for ablation and reference-distribution comparisons are included.
`scenario_generator_controlled.R` retains the separate Array A design.
Cluster scheduling is infrastructure-specific; these scripts do not require
the author's personal Slurm account or scratch directory.

The historical full grid includes internal controls as well as formal
comparators. Do not turn multiple budgets from the same dataset into independent
replicates, or treat single-seed ncurve/ablation as three-seed runs.

## Corrected historical results

The bundled CSVs contain corrected values. Their raw RDS inputs are not bundled.
If those inputs are available in `results/ncurve`, `results/abl`, `results/hart`
and `results/baronfix`, together with the actual corresponding scenario tables,
run from `benchmark/`:

```sh
Rscript difs_baronfix_merge.R
Rscript difs_hart_compare.R corrected=../results/diag/baron_hart_repaired.csv
Rscript difs_tables_export.R
```

Do not overwrite the corrected bundled tables by blindly summarizing old raw
Baron results. The repair-aware exporters and figures use corrected sources.
Actual cluster-submission scenario files must be retained separately; generated
scenarios or test fixtures are not proof of the original submission history.

From the repository root, `Rscript scripts/rebuild_figures.R` uses the bundled
corrected tables and does not require raw RDS files or new clustering runs.
