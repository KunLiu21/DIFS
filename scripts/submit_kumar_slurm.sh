#!/usr/bin/env bash
# Source your existing DIFS environment first; invoke on the login node.
set -euo pipefail
: "${DIFS_ROOT:?Source your DIFS environment first}" "${DIFS_SIF:?}" "${DIFS_ACCOUNT:?}"
repo_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
input_file=${1:-${DIFS_SOURCE:-$DIFS_ROOT/source}/Kumar.rds}
if [ ! -f "$input_file" ] || [ ! -f "$DIFS_SIF" ]; then
  echo "Missing input or container: $input_file ; $DIFS_SIF" >&2; exit 2
fi
export DIFS_EXAMPLE_REPO="$repo_dir"
export DIFS_EXAMPLE_INPUT=$(realpath "$input_file")
sbatch --parsable --job-name=difs_kumar_demo --partition=normal \
  --account="$DIFS_ACCOUNT" --time=02:00:00 --nodes=1 --ntasks=1 \
  --cpus-per-task=2 --mem=16G --export=ALL --chdir="$repo_dir" \
  --output="$repo_dir/kumar-slurm-%j.log" <<'WORKER'
#!/usr/bin/env bash
set -euo pipefail
cd "$DIFS_EXAMPLE_REPO"
out="example/output-slurm-${SLURM_JOB_ID:?}"
echo "Git revision: $(git rev-parse HEAD)"
echo "Output: $out"
bash scripts/run_kumar_container.sh "$DIFS_SIF" "$DIFS_EXAMPLE_INPUT" "$out" \
  2>&1 | tee "kumar-console-${SLURM_JOB_ID}.log"
runner=$(command -v singularity || command -v apptainer)
"$runner" exec --bind "$PWD:/difs" --pwd /difs "$DIFS_SIF" \
  Rscript scripts/verify_kumar_output.R "$out" \
  2>&1 | tee "kumar-verification-${SLURM_JOB_ID}.log"
echo "PASS: full example execution and output verification completed."
WORKER
