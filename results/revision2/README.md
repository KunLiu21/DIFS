# Organized second-round result snapshot

These 55 CSV files are copied byte-for-byte from the author's organized
2026-10-04 revision package. SHA256 identities are in
[`provenance/revision2_sources.json`](../../provenance/revision2_sources.json).
They do not overwrite the v2.0 result archive under the other `results/` paths.

**Corrected statistics / untagged material update:**
`02_ablation/rev2_abl2_comparisons.csv` and `08_markers/marker_coverage_tests.csv`
include the author's 2026-10-04 precision correction. Paired differences are
rounded to 12 decimal places before zero removal and midranking. This prevents
CSV/RDS serialization noise from becoming spurious differences or broken ties.
The hurdle-marker, 100-feature, nine-dataset comparison against Seurat now has
p=0.4140625. All 72 component and 768 marker comparisons were independently
checked from observations by `scripts/verify_revision2_statistics.py`.
The per-run observations were not altered or rerun. This update alone is not a
release tag, manuscript acceptance check, or independent full-grid rerun.

## Scope and denominators

- The source CSVs preserve the historical arms: S0 has **2601 rows**, and its
  main grid has **2320 observed of 2340 planned**, 20 missing. The workbook
  builder removes `feature_method == "Variance"` for the manuscript scope.
- After that explicit filter: **2115 S0 rows** = 1852 main-grid + 117 ablation +
  52 Hartigan + 94 Baron repair provenance rows; the main grid is **1852/1872**.
  Also filter `is_analysis_row` before aggregation to avoid counting repair
  provenance as extra experiments. No new comparator experiment is introduced.
- S0b: **2097 records** = 234 extended-ablation + 1553 repeats + 310 k-offset
  records. The latter include 12 `skipped_k_below_2` records, leaving 298 completed
  k-offset results. S5b lists 9 missing tasks (7 repeats, 2 k offsets).
- E8: 300 rows (60 panels times 5 rankings).
- Full-method simulation: 1263 rows, including the completed rerun of task 25.
  Missing FEAST outputs remain missing; do not treat each row as an independent panel.
- Marker per-dataset table: 3120 rows; marker tests: 768 source-summary rows.

## Locate an analysis

| Directory | Contents |
|---|---|
| `01_main_grid` | Corrected main grid, n-curves and ranks |
| `02_ablation` | Original component decomposition and extended component comparisons |
| `03_reps_krange` | Repeated-run means/spread and cluster-number comparisons |
| `04_repairs` | Baron, SC3 and Hartigan comparison provenance |
| `05_feature_sets_gate_ranking` | Repaired overlap, gate diagnostics, ranking and RBM39 |
| `06_runtime` | Runtime and recorded scheduler memory summary |
| `07_simulation` | E1–E8, complete-method panels and null-distribution illustration |
| `08_markers` | Reference markers, per-dataset coverage, summary and paired tests |
| `supplementary_tables_S0-S5` | S0–S5 plus S0b/S5b |

The workbook builder maps overlap to S6, runtime to S7, simulations to S8a/S8b,
and marker coverage/tests to S9/S9b. It leaves non-Baron memory values empty as
in the manuscript workbook. Absence of a reliable peak is not zero memory usage.

Raw expression matrices, private cluster logs, the response letter and manuscript
drafts are not part of this snapshot. The older null-table file is not substituted:
all public runs use `inst/extdata/dip_null_tables.rds` (MD5
`29406eeaabe9806321450b034074c868`).
