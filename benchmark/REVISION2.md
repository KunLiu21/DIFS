# Second-round revision materials (2026-10-04 source snapshot)

These additions cover the extended component analysis, repeated runs, cluster
number sensitivity, runtime/memory, marker recovery and additional simulations.
Sources are mapped to the author's organized `DIFS_revision_package` in
[`provenance/revision2_sources.json`](../provenance/revision2_sources.json).
The validated Kumar implementation and existing benchmark functions are unchanged.

**Status: untagged material update.** The component and marker summaries include
the 2026-10-04 correction for numerical zero/tie handling, with a documented
12-decimal difference convention. Read [the result audit status](../results/revision2/README.md)
before quoting inferential results. No main-grid clustering rerun was performed
as part of this packaging. The manuscript and response letter are not published here.

## What each script does

| Analysis | Entry point | Archived evidence |
|---|---|---|
| Additional ablation, 5 seeds, 5 stratified subsamples, k offsets | `scenario_generator_rev2.R`, `Method_comparison_rev2.R`, `difs_rev2_overlay.R` | S0b, S5b, `rev2_*` tables |
| Analyze those arrays | `difs_rev2_analysis.R` | Comparisons with corrected signed-rank helper |
| Runtime / scheduler memory | `difs_runtime_memory.R` | `06_runtime/runtime_scaling.csv` |
| Hurdle, level, detection and Wilcoxon marker references | `difs_marker_coverage.R` | `08_markers/`, S9/S9b |
| Full DIFS on simulated panels | `../simulation/difs_sim_complete.R` | `07_simulation/sim_complete_*.csv` |
| Stage-I assumption experiment | `../simulation/difs_sim_E8_assumption.R` | `07_simulation/E8_assumption.csv` |
| Gaussian vs uniform null illustration | `../simulation/difs_fig_dip_null.R` | `07_simulation/fig_dip_null_pvalues.csv` |
| Repaired feature exports and cross-dataset overlap | `difs_export_features_v2.R`, `difs_feature_overlap.R` | `05_feature_sets_gate_ranking/` |
| Ranking diagnostics and RBM39 | `difs_ranking_check.R`, `difs_key_check.R`, `difs_gene_rank.R` | Same directory |
| Earlier SC3 repair sensitivity | `scenario_generator_sc3fix.R`, `Method_comparison_sc3fix.R`, `difs_sc3fix_compare.R` | `04_repairs/`; historical diagnostic, not the final main-grid runner |
| S0b/S5b from raw runs | `difs_tables_export_rev2.R` | `supplementary_tables_S0-S5/` |
| Recompute corrected comparisons from public CSVs | `../scripts/regen_tables_20261004.R` | 72 component and 768 marker comparisons |
| Supplementary workbook | `../scripts/build_additional_file1.py` | S0 through S9b, including companion sheets |

The original Monte Carlo illustration is not the manuscript Figure 4 script.
Figure 4 uses the bundled Gaussian lookup and `difs_fig_dip_null.R`.

## Prepare a separate run directory

From the repository root, Python 3.9+:

```sh
python scripts/prepare_revision2_workspace.py
```

This creates `work/revision2/`, copies the analysis sources, the bundled lookup
table and the corrected baseline CSVs, and refuses conflicting existing files.
It does not download expression matrices or run an experiment. Set `--workspace`
to another directory if desired; then set `DIFS_WORKSPACE` to that same absolute
directory when using the shell runner. This keeps generated results separate
from archived evidence and from v2.0 tables.

Install the dependencies in [environment/README.md](../environment/README.md),
then prepare the datasets using `prepare_datasets.R` from the workspace's
`benchmark/` directory. That directory must contain all 13 prepared inputs under
`../source/` before generating the complete second-round scenario tables.

```sh
cd work/revision2/benchmark
Rscript prepare_datasets.R only=Kumar     # one dataset; not preparation of all 13
# Prepare each remaining dataset before the following two commands.
Rscript scenario_generator_rev2.R
Rscript difs_rev2_smoke.R
```

Inspect the generated task table and smoke results before submitting expensive
arrays. The generated scenarios describe a new run, not historical job receipts.
From the repository root, these examples run only the explicitly selected task:

```sh
bash scripts/run_revision2.sh rev2 set=abl2 k=1
bash scripts/run_revision2.sh markers task=list
bash scripts/run_revision2.sh simulation task=list
# E8 runs all of its panels; submit with appropriate resources rather than on a login node.
# bash scripts/run_revision2.sh e8
```

The marker task list also checks that all prepared dataset files exist; it
deliberately fails with their names until dataset preparation is complete.

Optional Slurm wrappers are `benchmark/run_rev2.sh`,
`benchmark/run_marker_coverage.sh`, and `simulation/run_sim_complete.sh`. Supply
your own account, partition, resources, output paths and array indices on the
`sbatch` command. Create log directories before submission. If using a compatible
Singularity image, export `DIFS_SIF`; otherwise the runner uses Rscript directly.
An optional `DIFS_LIB` is your own R library path. No author's account or private
container is required by these wrappers.

## Rebuild tables without pretending to rerun clustering

The raw per-run RDS and prepared input matrices are not included. Analyses that
read them require a corresponding local run or the original author's files and
scenario tables. The archived CSVs can rebuild the supplementary workbook:

```sh
python scripts/build_additional_file1.py --out work/Additional_file_1_snapshot.xlsx
```

Requires pandas and openpyxl. The builder refuses an existing output path, and
the workbook README records the corrected numerical convention and release status.
To regenerate the two corrected comparison tables into a separate empty directory,
and independently verify means, denominators, win counts and exact p-values:

```sh
Rscript scripts/regen_tables_20261004.R work/recomputed-tables
python scripts/verify_revision2_statistics.py
python scripts/verify_revision2_snapshot.py
```

The R regeneration uses the same comparison definitions as the author scripts;
the Python check independently uses integer doubled ranks and subset-sum counts.
Neither command reruns feature selection or clustering. Confidence intervals
remain paired t intervals; multiplicity adjustments are not added by this update.

## Interpretation boundaries

Repeated seeds/subsamples are nested within datasets, not independent datasets.
The nine-dataset sensitivity analysis excludes the two Kumar-derived simulations
and two Zheng mixtures. Unadjusted p-values are not evidence of simultaneous
significance across all comparisons. Nonsignificance does not establish equivalence.

E8 matches variance on the latent log-expression scale; this does not imply equal
VST scores after count generation and normalization. Its scale factor is the
median library size. E2 continuous nulls start on the log-expression scale, while
count nulls are normalized with background genes. Neither experiment establishes
universal calibration. Stage II requires at least 10% of cells in some preliminary
cluster to exceed the expression threshold; it does not test literally every gene.
