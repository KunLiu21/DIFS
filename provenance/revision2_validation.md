# Revision-2 packaging validation (2026-10-04)

The source package was refreshed during the sync to include the author's
2026-10-04 12-decimal signed-rank corrections. The earlier uncorrected snapshot
was not committed or published as revision-2 material.

Local checks completed:

- All 79 imported/updated source and result hashes match the provenance manifest.
- Main-grid, supplementary, missing/skipped, simulation and marker denominators
  match `results/revision2/README.md` (including the explicit manuscript filter).
- Independent Python subset-sum checks reproduce all 72 component and 768 marker
  comparison means, pair counts, win/loss counts and exact p-values.
- The portable R regeneration script reproduces both corrected comparison CSVs
  including paired-t confidence interval columns (tolerance 1e-12).
- The workbook builder reconstructs all 15 data sheets of the updated manuscript
  workbook within 1e-12; its README additionally identifies release/audit status.
  Attempting to overwrite the existing output is rejected without changing it.
- All 45 benchmark/simulation R sources parse. All newly imported function
  definitions agree with the author's organized sources. Of 207 shared function
  comparisons, four pre-existing public packaging differences remain: finite-input
  rejection in `difs_signrank` / `difs_ci`, and checks in `read_runs` / `xcheck`
  in the existing repair-aware exporter. These were not replaced with weaker
  author-workspace checks. The only shared helper update adds numerical rounding.
- The existing lightweight method/metric tests and definition-loader tests pass;
  the new zero/tie/nonfinite regression tests pass.
- The full-simulation task list enumerates 85 panels without loading heavy
  clustering dependencies. The marker task list correctly refuses to proceed
  until all 13 prepared datasets exist.
- The null-figure script runs locally and reproduces all 1000 archived numerical
  rows within 1e-12. No visual-layout certification is claimed here.
- Generic shell entry points pass syntax checks; private account/path defaults
  were removed from imported runnable material.

These checks validate packaging, summary arithmetic and the stated example.
They are not a new run of all feature selection/clustering analyses. Raw matrices
and per-run RDS inputs are not redistributed. Repeated-run independence,
selection bias and method scope remain scientific limitations, not issues a
successful packaging check can remove. See the fresh Kumar evidence separately.
