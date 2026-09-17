# DIFS — Discriminative Feature Selection for single-cell RNA-seq clustering

Code for the manuscript *DIFS: Discriminative Feature Selection for Cell
Clustering Based on Single-Cell RNA Sequencing Data* (BMC Bioinformatics, under
revision).

DIFS selects genes in two stages. **Stage I** keeps genes whose expression,
*among the cells that express them*, departs from unimodality, scored with the
dip statistic of Hartigan & Hartigan calibrated against a Gaussian reference at
each gene's own number of expressing cells. **Stage II** adds cell-type-specific
genes — which are unimodal among expressing cells and therefore invisible to
stage I — using a Fisher exact test against a preliminary clustering.

---

## Quick start — 10 seconds, no Bioconductor, no network

```bash
git clone https://github.com/KunLiu21/DIFS.git
cd DIFS
Rscript demo/run_demo.R
```

Requires only **R ≥ 4.1** and **diptest** (`install.packages("diptest")`).

The demo runs stage I on a small simulated panel and checks three properties of
the method. Its output is committed at [`demo/expected_output.txt`](demo/expected_output.txt)
so you can diff against it:

```
gene class                      gated  in top 30  median rank
BIMODAL (stage I target)       30/30          30         15.5
ONOFF (stage II target)        30/30           0          248
background                    368/540          0        227.5
```

Read that table before anything else — it is the method's design, measured.
Stage I recovers **every** bimodal gene in the top 30, and ranks **on/off
markers no better than background**. That is not a failure: an on/off marker is
unimodal among the cells that express it, so stage I cannot see it, and
recovering those genes is exactly what stage II is for. If you only remember one
thing about DIFS, remember that the two stages target disjoint regimes.

If **Seurat** is installed the demo also clusters on the selected features and
reports an adjusted Rand index; if it is not, the demo says so and exits 0.

---

## What is here

| Path | What it is | Needs |
|---|---|---|
| `R/difs_core.R` | The method: expression gate + stage I ranking. Self-contained. | R, diptest |
| `inst/extdata/dip_null_tables.rds` | Precomputed null quantiles of the dip statistic (446 KB) | — |
| `demo/` | The 10-second example above, its data, and its expected output | R, diptest |
| `benchmark/` | The harness that produced the paper's tables | see below |

`R/difs_core.R` is copied verbatim from the benchmark harness, so the two cannot
drift. If you only want to *use* DIFS, that one file and the null table are
everything you need.

### Using DIFS on your own data

```r
source("R/difs_core.R")

# logmat: genes x cells, log-normalised as log(1 + 1e4 * count / library size)
ranked <- difs_stage1_ranking(logmat,
                              min.expression = log(5),   # a cell "expresses" a
                                                         # gene above 4 per 10k
                              gate_rule = "combined",
                              gate_k    = k)             # expected cluster count
features <- head(ranked, 300)
```

`gate_rule = "combined"` is the rule used in the revision:
a gene must exceed `log(5)` in at least `max(min(0.05n, n/(4k)), 30)` cells.
`gate_rule = "submitted"` reproduces the originally submitted rule,
`max(0.05n, 35)`.

---

## Reproducing the paper

**This part needs a cluster. We are explicit about that rather than implying a
laptop will do.** The benchmark is 13 datasets × 5 feature-selection methods ×
2 clustering methods × 9 feature budgets × 2 policies for the number of
clusters = **2,340 runs**, plus a 117-run ablation.

| | |
|---|---|
| Longest single run | ~2.5 h (SC3 clustering on Baron, 8,569 cells) |
| Peak memory | see `benchmark/` notes; 64 GB requested per task |
| Scheduler | SLURM array jobs |
| Extra dependencies | Seurat, SingleCellExperiment, SC3, FEAST, monocle, TSCAN, mclust |
| Datasets | fetched by `benchmark/prepare_datasets.R` from `DuoClustering2018` and `scRNAseq`; **not redistributed here** |

Order:

```bash
Rscript benchmark/prepare_datasets.R          # builds ../source/*.rds
Rscript benchmark/scenario_generator_ncurve.R # writes the scenario table
sbatch  benchmark/run_ncurve.sh               # the array
Rscript benchmark/results_summary_ncurve.R    # tables
Rscript benchmark/difs_figures.R              # figures
```

Raw per-run `.rds` files (2,340 of them) are **not** in this repository. The
summary tables the figures are built from are, under `results/`.

---

## Verification status

We say what was actually checked, and where.

| Check | Status |
|---|---|
| `demo/run_demo.R` tier 1 (stage I, diptest only) | **verified**, R 4.3.3 / diptest 0.77.0, output committed |
| `demo/run_demo.R` tier 2 (Seurat clustering) | **not yet verified in a clean environment** |
| Full benchmark from a clean clone | **not yet verified**; the paper's results were produced on the cluster with the code in `benchmark/` |

A "runs end to end" claim is only worth something if someone actually ran it, so
the tiers that have not been re-run from a clean clone are marked as such.

---

## Corrections to the originally submitted code

The version tagged `v1.0-submitted` is what accompanied the first submission.
The following defects were found during the revision and are fixed here. They
are listed because at least one of them changes results.

*(This section is completed as each fix lands; see CHANGELOG.md.)*

---

## Citation

If you use this code, please cite the manuscript. See `CITATION.cff`.

## Licence

MIT — see [LICENSE](LICENSE). The packages DIFS calls carry their own licences.
