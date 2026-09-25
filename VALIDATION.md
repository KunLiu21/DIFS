# Validation status

This file distinguishes code preparation from execution evidence.

| Item | Status |
|---|---|
| Complete two-stage Kumar entry point | Implemented; full execution pending |
| Kumar historical reference | Extracted from the verified main-grid S0 row |
| Original full-grid clean rerun | Not performed for this packaging task |
| Corrected S0–S5 | Previously validated against independent RDS decoding |
| SC3 repair outputs | Previously validated for 112 configurations / 8569 cells |
| Local lightweight/API/figure checks | See provenance/local_validation.txt after execution |
| Public final release / DOI | Not yet issued for this working branch |

Local R has diptest and an installed Seurat package, but loading Seurat is blocked
by Windows Application Control on a data.table DLL. SC3, FEAST and
DuoClustering2018 are also unavailable. This restriction was not bypassed.
The full example must run in a suitable environment. Missing required packages
cause an error; passing the quick demo does not certify the full example.

No fabricated selected-gene list, success log or expected full-example output
is included. The recorded historical benchmark row is explicitly labeled as
historical. Send the generated full-example output and session information for
verification before changing this status to completed.

## Cluster attempt on 2026-09-25

Job 36842830 at commit 1576ca4 failed at the new API input validation before
feature selection. The error combined several checks and did not identify the
trigger. Inspection found an integer-only requirement absent from benchmark
preparation. This restriction is removed; dense and sparse fractional estimates
are preserved, while negative and nonfinite values still fail regression tests.
The runner now records input diagnostics and exposes errors in both run and
Slurm logs. A cluster retry is required; full-example success is not established.
