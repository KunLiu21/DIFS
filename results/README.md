# Corrected tabular evidence

`supp/S0_per_run.csv` has 2601 rows: 2320 main-grid, 117 ablation, 52 Hartigan
and 112 raw Baron repair-source rows. Select `is_analysis_row=TRUE`, then the
appropriate `result_set`; do not pool experiments or double-count repair rows.
`baron_repaired` identifies substitutions, `source_result_set` their source, and
`ARI_before_repair` preserves the previous metric. All four metrics and the
saved stage metadata were replaced together.

S1–S4 contain 2320 main-grid rows each, for ARI, FM, cell-pair Jaccard and purity.
S5 contains 20 missing configurations. Main figures use DIFS/FEAST/Seurat; other
stored arms remain for provenance and internal analyses. Array A and simulation
outputs are separate and are not disguised as part of these six tables.

`diag/grid_corrected.csv`, `diag/abl_corrected.csv`, and the corrected Hartigan
tables support summaries and figures. This package does not include raw count
matrices or the original per-run RDS files. S0 stores recorded metrics, not the
full cell-level predictions needed to recalculate each metric from scratch.
