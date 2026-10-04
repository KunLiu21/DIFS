#!/usr/bin/env bash
# Run one explicitly requested task. Choose Slurm resources at submission time.
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
workspace="${DIFS_WORKSPACE:-$repo/work/revision2}"
if [[ ! -f "$workspace/benchmark/Functions_controlled.R" ]]; then
  printf '%s\n' 'Prepare an isolated workspace with scripts/prepare_revision2_workspace.py first.' >&2
  exit 2
fi
mode="${1:?Usage: run_revision2.sh <rev2|markers|simulation|e8> [R key=value arguments]}"
shift
case "$mode" in
  rev2) script='Method_comparison_rev2.R' ;;
  markers) script='difs_marker_coverage.R' ;;
  simulation) script='../simulation/difs_sim_complete.R' ;;
  e8) script='../simulation/difs_sim_E8_assumption.R' ;;
  *) printf 'Unknown mode: %s\n' "$mode" >&2; exit 2 ;;
esac
cd -- "$workspace/benchmark"
export DIFS_NULL_TABLES="$(cd ../inst/extdata && pwd)/dip_null_tables.rds"
export DIFS_TABLES="$DIFS_NULL_TABLES"
runner=(Rscript)
if [[ -n "${DIFS_SIF:-}" ]]; then
  runner=(singularity exec --bind "$workspace:$workspace" "$DIFS_SIF" Rscript)
fi
exec "${runner[@]}" "$script" "$@"
