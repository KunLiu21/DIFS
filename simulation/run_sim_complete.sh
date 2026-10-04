#!/usr/bin/env bash
set -euo pipefail
repo="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
exec bash "$repo/scripts/run_revision2.sh" simulation "k=${SLURM_ARRAY_TASK_ID:?Submit as a Slurm array}"
