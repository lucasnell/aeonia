#!/bin/bash -l
#SBATCH --array=1-6
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=1-20:00:00
#SBATCH --job-name=sobol-sens
#SBATCH --output=sobol-sens-%a.out
#SBATCH --error=sobol-sens-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL


./sobol-sens.R "${SLURM_CPUS_PER_TASK}" "${SLURM_ARRAY_TASK_ID}"

