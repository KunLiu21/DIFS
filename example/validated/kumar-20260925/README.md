# Executed Kumar example, 2026-09-25

These are actual outputs from Slurm job 36842924, using source commit
`1c64256c0c395d42c9f289cfe8a6758c37524651`. Later release changes archive the
evidence and update documentation/citation; they do not modify this executed
algorithm. The uploaded execution checkout has no tracked changes relative to
the recorded commit.

Configuration: prepared Kumar input (246 cells, 44,508 genes), true-k=3,
100 requested features, outer seed=1, repaired SC3 seed=1 / two cores,
Gaussian lookup, submitted gate, and final Refined Louvain.

Observed output: 100 selected genes, stage-I size 453, stage-II candidate size
616, ratio 1:1, all 246 cells labeled. Four metrics agree with the historical
reference within 5e-16. See `summary.csv`, `historical_comparison.csv` and
`independent-check.json`.

`source-manifest.json` gives SHA256 identities and sizes of files copied from
the execution, before adding this README and the independent check. Original
file bytes are retained in Git, including original log line endings.
`difs_fit.rds` retains preliminary/final labels, selected genes and SC3 audit.
`cell_labels.csv` permits metric recalculation without R package dependencies.

From the repository root, verify these saved outputs with:

```sh
Rscript scripts/verify_kumar_output.R example/validated/kumar-20260925
```

The cluster verification was repeated locally and the four metrics were also
recomputed independently from contingency counts in Python. Verifying saved
outputs is not a new execution of feature selection/clustering.

The input MD5 is `d0dbd41cb5e691d7ffa2661004b1238d`; the expression matrix is not
redistributed. The example documentation supplies its DuoClustering2018/Kumar
source and preparation route. This run used an existing prepared object and
existing Linux/R environment; fresh download/preparation and fresh installation
have not been validated. The observed environment inventory is
`environment/validated-session-packages.csv` in the repository.
