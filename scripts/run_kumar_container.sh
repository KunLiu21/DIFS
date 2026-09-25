#!/usr/bin/env bash
# Invoke explicitly using an existing, trusted container image.
# bash scripts/run_kumar_container.sh /path/to/image.sif /path/to/Kumar.rds
set -euo pipefail
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
  echo "Usage: bash scripts/run_kumar_container.sh /path/to/image.sif /path/to/Kumar.rds [new-output-directory]" >&2
  exit 2
fi
repo_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
image_path=$1
input_path=$2
output_path=${3:-example/output-container}
if [ ! -f "$image_path" ] || [ ! -f "$input_path" ]; then
  echo "The container and Kumar input must both exist." >&2; exit 2
fi
input_dir=$(cd "$(dirname "$input_path")" && pwd)
input_name=$(basename "$input_path")
runner=$(command -v singularity || command -v apptainer) || { echo "Singularity/Apptainer not found" >&2; exit 2; }
cd "$repo_root"
library_args=()
if [ -n "${DIFS_LIB:-}" ]; then
  if [ ! -d "$DIFS_LIB" ]; then echo "DIFS_LIB is not a directory" >&2; exit 2; fi
  library_args=(--bind "$DIFS_LIB:/difs-lib:ro" --env R_LIBS=/difs-lib --env R_LIBS_USER=/difs-lib)
fi
"$runner" exec "${library_args[@]}" --bind "$repo_root:/difs" --bind "$input_dir:/difs-input:ro" \
  --pwd /difs "$image_path" Rscript example/run_kumar.R \
  "input=/difs-input/$input_name" "out=$output_path"
