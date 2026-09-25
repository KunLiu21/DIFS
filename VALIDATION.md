# Validation status

This file distinguishes code preparation from execution evidence.

| Item | Status |
|---|---|
| Complete two-stage Kumar entry point | Executed and verified; job 36842924, commit 1c64256 |
| Kumar historical reference | Extracted from the verified main-grid S0 row |
| Original full-grid clean rerun | Not performed for this packaging task |
| Corrected S0–S5 | Previously validated against independent RDS decoding |
| SC3 repair outputs | Previously validated for 112 configurations / 8569 cells |
| Local lightweight/API/figure checks | See provenance/local_validation.txt after execution |
| Release scope | Validated single-example revision; no DOI claimed |

Local R has diptest and an installed Seurat package, but loading Seurat is blocked
by Windows Application Control on a data.table DLL. SC3, FEAST and
DuoClustering2018 are also unavailable. This restriction was not bypassed.
The full example subsequently ran in the author's compatible cluster environment.
Missing required packages cause an error; passing the quick demo alone does not
certify the full example.

The recorded historical benchmark row is explicitly labeled as historical.
Actual selected genes, cell labels, metrics, fitted output, logs and session
information from the successful run are included under
`example/validated/kumar-20260925/`.

## Cluster attempt on 2026-09-25

Job 36842830 at commit 1576ca4 failed at the new API input validation before
feature selection. The error combined several checks and did not identify the
trigger. Inspection found an integer-only requirement absent from benchmark
preparation. This restriction is removed; dense and sparse fractional estimates
are preserved, while negative and nonfinite values still fail regression tests.
The runner now records input diagnostics and exposes errors in both run and
Slurm logs. The retry outcome is recorded below.

## Verified retry on 2026-09-25

Job 36842924 executed commit `1c64256c0c395d42c9f289cfe8a6758c37524651`.
The uploaded checkout has no tracked source changes relative to that commit.
The full run produced 100 genes, preliminary and final labels for 246 cells,
stage-I size 453, stage-II candidate size 616, and mixing ratio 1:1. Its SC3 audit
records zero missing labels. ARI=0.988771468816226, FM=0.992545322773649,
Jaccard=0.985200431245712 and Purity=0.995934959349594 agree with the historical
row to within 5e-16.

The saved-output R verifier passed on the cluster and again locally. Independent
Python contingency-table calculations from `cell_labels.csv` also reproduced all
four metrics; these checks did not rely on the printed success statement.
Input diagnostics record 6,174,993 fractional stored values with no negative or
nonfinite values, confirming why the initial integer-only wrapper rejected this
input. The source expression matrix was not rounded or altered.

The existing prepared input's MD5 and null-table MD5 are recorded. The source
matrix is not redistributed. A fresh data download/preparation, a fresh package
installation, and a complete rerun of all paper experiments have not been
performed as part of this release validation. A session package inventory is an
execution record, not a tested dependency lockfile.
