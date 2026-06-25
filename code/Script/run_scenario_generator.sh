#!/bin/bash

### set the number of cores
#SBATCH -n 1
### set the expected runtime
#SBATCH --time=40:00:0
### set the node id
#SBATCH -p normal
### group account
#SBATCH --account=coa_lji226_uksr 

### set the output file name
#SBATCH --output=../Rout/slurm-%j.out

### job array this can be used to automatically generate an array of task id with value 1, 2 and 3
#SBATCH --array=1 # job array index

# Run this script as follow:
# sbatch runR.sh
# Then the values 1 or 2 or 3  will be assigned to the n argument in each job.

export TMPDIR=~/back

# Load the Singularity module
singularity exec "../../method_comparison/code/hotpot_rstudio-archR.sif" Rscript --verbose scenario_generator.R k=${SLURM_ARRAY_TASK_ID} >../Rout/z${SLURM_ARRAY_TASK_ID}.Rout