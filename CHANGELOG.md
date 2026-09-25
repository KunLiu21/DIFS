# Revision working branch

- Validate the full Kumar prepared-input example on the cluster (job 36842924),
  then recheck labels/metrics locally and with independent contingency counts.
  Archive actual outputs, original logs, checksums and observed package versions.
  All four metrics match the historical row within 5e-16. This is a single
  configuration in an existing environment, not a clean full-grid rerun.
- Remove the example API's extra integer-only input restriction; retain
  fractional estimates unchanged as benchmark preparation does. Distinguish
  negative, nonnumeric and nonfinite inputs. Add dense/sparse regression checks,
  input diagnostics, visible error logs and a Slurm retry with a fresh output
  directory per job. The first full-example attempt (36842830) failed validation;
  the subsequent retry is now validated as recorded above.
- Add a complete two-stage API and a Kumar worked-example runner using shared
  benchmark definitions; retain the old code for historical reference.
- Correct the public documentation's gate from combined to submitted, matching
  the manuscript and corrected benchmark. Identify the old demo as stage I only.
- Include the SC3 hybrid-label completion used for the Baron repair and preserve
  its recorded defaults. Do not silently discard unassigned cells.
- Add corrected supplementary and diagnostic CSVs, explicit missing records,
  original/repaired provenance, and fixed figure-C dataset labels (13/9).
- Include repair-aware exporters, paired Hartigan checks and exact signed-rank
  code; provide scripts to regenerate figures from the bundled summaries.
- Add environment guidance, real-data preparation, output and identity checks,
  and a transparent validation status. Full Kumar execution is verified for the
  documented configuration; no revision release DOI is claimed.

This is not a claim that all original pipeline limitations are repaired. Feature
count exceptions, post-selection inference, random-start coverage and resource
measurements remain as documented in the revision.
