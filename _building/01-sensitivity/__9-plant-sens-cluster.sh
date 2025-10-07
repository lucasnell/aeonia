#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=80G
#SBATCH --time=1-20:00:00
#SBATCH --job-name=9plant-sens
#SBATCH --output=9plant-sens.out
#SBATCH --error=9plant-sens.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

# For an interactive job, you could use the following to replace above:
# srun -n 1 -c 50 -N 1 --mem=80G --time=1-20:00:00 \
#     --job-name="9plant-sens" --pty bash -l


chmod +x "9plant-sens-cluster.R"
Rscript --vanilla "9plant-sens-cluster.R" "$SLURM_CPUS_PER_TASK"

