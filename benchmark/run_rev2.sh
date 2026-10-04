#!/usr/bin/env bash
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
set_name="${1:?Specify abl2, reps, or krange}"
case "$set_name" in abl2|reps|krange) ;; *) exit 2 ;; esac
exec bash "$repo/scripts/run_revision2.sh" rev2 "set=$set_name" "k=${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array}"
