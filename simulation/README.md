# Simulation analysis sources

Run from this directory with a compatible R and diptest installation:

```sh
Rscript difs_sim_run.R preset=full gate=implemented out_dir=../results/simulation-new
Rscript difs_sim_report.R out_dir=../results/simulation-new
```

The named implemented gate uses ln(5), while documented selects the historical
function-signature value 0.5. Specify the gate explicitly. The full preset can
be expensive; it is not part of the lightweight validation or the single-data
worked example. Use a new output directory to avoid confusing cached null
tables or previous runs with fresh simulations.

The simulator uses separate simulation-specific functions and its own cached
null tables. Those are not the production stage-I lookup bundled in inst/.
Its entire execution has not been rerun during this packaging task. Scenario
parameters, seeds and input table identity must accompany a reproduction claim.
