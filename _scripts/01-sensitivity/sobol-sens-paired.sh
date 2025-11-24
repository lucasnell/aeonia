#!/bin/bash -l

#SBATCH --array=1-8
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=40G
#SBATCH --time=1-20:00:00
#SBATCH --job-name=sobol-sens-paired
#SBATCH --output=logs/sobol-sens-paired-%a.out
#SBATCH --error=logs/sobol-sens-paired-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

#'
#' Small-scale sensitivity via Sobol indices
#'
#'

#' I first moved the sensitivity preamble and this script over to bioHPC
#' using the following:
#'
#' cd ~/GitHub/Cornell/aeonia/_scripts/01-sensitivity \
#'     && scp sobol-preamble.R sobol-sens-paired.sh \
#'     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/sobol-paired/
#'
#' This was then run on BioHPC in a non-interactive job started with the following:
#'
#' cd /home2/lan68/sobol-paired \
#'     && sbatch sobol-sens-paired.sh
#'





Rscript - << EOF

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.4")

n_threads <- suppressWarnings(as.integer("${SLURM_CPUS_PER_TASK}"))
if (is.na(n_threads)) stop("SLURM_CPUS_PER_TASK must be an integer")

curr_idx <- suppressWarnings(as.integer("${SLURM_ARRAY_TASK_ID}"))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")

max_idx <- suppressWarnings(as.integer("${SLURM_ARRAY_TASK_MAX}"))
if (is.na(max_idx)) stop("SLURM_ARRAY_TASK_MAX must be an integer")

# Set threads for simulations:
options("mc.cores" = n_threads)

source("sobol-preamble.R")

# Verify that indices align with 'sobol_mat' object created in 'sobol-preamble.R':
stopifnot(nrow(sobol_mat) %% max_idx == 0L)

# Get seed for this set of simulations:
.seed <- c(834704589, 1938725610, 1927014495, 1236399649, 741338666, 820277141,
           482712696, 587689356)[curr_idx]

# number of rows to simulate per job:
n_rows_pj <- nrow(sobol_mat) %/% max_idx

# Start and stop for rows to simulate this job:
curr_start <- (curr_idx - 1L) * n_rows_pj + 1L
curr_stop <- curr_idx * n_rows_pj

# Takes ~6.5 min per job (of 8 total jobs) with 9 parameters, N=2^12, and 20 threads:
cat("Starting simulations...\n")
t0 <- Sys.time()
set.seed(.seed)
sobol_sims <- some_sobol_sims(curr_start:curr_stop)
t1 <- Sys.time()
difftime(t1, t0, units = "min")

cat("Writing output...\n")
t0 <- Sys.time()
write_rds(sobol_sims, sprintf("sobol-sims-paired-%i.rds", curr_idx), compress = "gz")
t1 <- Sys.time()
difftime(t1, t0, units = "hours")

EOF
