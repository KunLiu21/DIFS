# Fresh-installation Kumar reproduction, verified 2026-10-04

GitHub Actions run [36806258451](https://github.com/KunLiu21/DIFS/actions/runs/36806258451)
executed commit `368959ca80ab70b19a28b54ad47deb7d5aae26d5` on Ubuntu 22.04 / R 4.4.1.
The workflow installed packages into an empty library, downloaded public data
with empty data caches, prepared Kumar, and ran the full two-stage example.
No author-prepared input, package cache, private library or container was used.

All 246 cells were labeled and 100 genes were selected. Gene identities and order
match the earlier prepared-input run; preliminary and final partitions each have
ARI 1 against that run. Independent contingency calculations give:
ARI 0.9887714688162261, FM 0.9925453227736489, Jaccard 0.985200431245712,
Purity 0.9959349593495935. The saved-output verifier also passed.

`evidence/` retains installation, preparation and verification logs, package
inventories and empty-library/cache checks. `preparation/` records public-data
provenance; the input expression matrix is not redistributed. Original artifact
bytes are identified by `source-manifest.json`; independent checks are in
`independent-check.json`. The prepared-object MD5 differs from the older object
because object serialization/provenance is not asserted byte-identical.

From the repository root:

```sh
Rscript scripts/verify_kumar_output.R example/validated/kumar-clean-20261004
```

This verifies one configuration, not the complete benchmark or every platform.
An earlier clean run (36805242810) failed in Seurat GroupSingletons. The next
commit added failure diagnostics, without changing the scientific algorithm;
the successful run skipped those diagnostics. Its success does not establish the
cause of the earlier failure or guarantee that it cannot recur. Earlier setup
failures and their fixes remain in Git history and the Actions logs.
