#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-6
#SBATCH --cpus-per-task=11
#SBATCH --mem=50G
#SBATCH --time=04:00:00
#SBATCH --job-name=large-r_fly
#SBATCH --output=logs/large-r_fly-%a.out
#SBATCH --error=logs/large-r_fly-%a.err
#SBATCH --mail-type=END,FAIL

#
# Larger landscape simulations that vary r_fly
#

# I first moved this script and the preamble over to bioHPC using the following:
#
# cd ~/GitHub/Cornell/aeonia/_scripts
# ssh biohpc 'mkdir -p /home2/$USER/03-large-r_fly && mkdir -p /home2/$USER/03-large-r_fly/logs'
# scp ./02-large-hpc/03-large-r_fly.sh biohpc:/home2/$USER/03-large-r_fly/
#
# This was then run on BioHPC in a non-interactive batch job started with the following:
#
# cd /home2/$USER/03-large-r_fly/
# sbatch 03-large-r_fly.sh
#
# Then, when the jobs are done (assuming you're back on your desktop in
# directory `~/GitHub/Cornell/aeonia/_scripts`):
#
# scp biohpc:/home2/$USER/03-large-r_fly/large-r_fly*.rds ./interm-data/large-r_fly/
#
#





Rscript - << EOF

source("../00-large-hpc-preamble.R")

# --------------*
# Create dataframe of possible landscapes and parameter values:

sim_df <- sim_df |>
        # add row for without *Pseudomonas*
        add_row(n_pseudo = 0L, wt_vp = NA, wt_pp = NA,
                landscape = list(array(c(1L, rep(0L, 9999L)), c(100L, 100L, 100L)))) |>
        # add other parameters:
        crossing(crossing(wasp_resp = c("strong", "weak"),
                          radius = an_environ[["radius"]] * c(0.5, 1, 2),
                          p_load = 0.5,
                          sd_N = 0,
                          virus_attract = 1,
                          pseudo_repel = 1)) |>
        # placeholder for simulation output:
        mutate(sim = rep(list(NA), n()))



# --------------*
# Read inputs from job manager:

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")

max_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX")))
if (is.na(max_idx)) stop("SLURM_ARRAY_TASK_MAX must be an integer")
if (curr_idx > max_idx) stop("SLURM_ARRAY_TASK_ID cannot exceed SLURM_ARRAY_TASK_MAX")

# Get seed for this set of simulations:
.seed <- c(1866878966, 60212338, 1296054450, 372395175, 2130242721, 609579743)[curr_idx]

# Verify that indices align with 'sim_df':
stopifnot(nrow(sim_df) %% max_idx == 0L)

# number of rows to simulate per job:
n_rows_pj <- nrow(sim_df) %/% max_idx

# Start and stop for rows to simulate this job:
curr_start <- (curr_idx - 1L) * n_rows_pj + 1L
curr_stop <- curr_idx * n_rows_pj


sim_df <- sim_df[curr_start:curr_stop,]



# --------------*
# Main simulator function:



large_simmer <- function(sim_df_row) {

    radius <- sim_df[["radius"]][[sim_df_row]]

    landscape <- sim_df[["landscape"]][[sim_df_row]]
    wasp_resp <- sim_df[["wasp_resp"]][[sim_df_row]]
    p_load <- sim_df[["p_load"]][[sim_df_row]]
    sd_N <- sim_df[["sd_N"]][[sim_df_row]]
    virus_attract <- sim_df[["virus_attract"]][[sim_df_row]]
    pseudo_repel <- sim_df[["pseudo_repel"]][[sim_df_row]]

    args <- list(radius = radius,
                 landscape = landscape,
                 sd_N = sd_N,
                 virus_attract = virus_attract,
                 pseudo_repel = pseudo_repel,
                 Y0 = 250,
                 N0 = 55,
                 zeta = ifelse(wasp_resp == "weak", 0.1, 0.9),
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = 100L,
                 summ = "all")

    return(do.call(big_plantscape, args))

}


# Took ~30 min per job, where each job processed 21 rows using 11 threads

cat("Starting simulations...\n")
t0 <- Sys.time()

set.seed(.seed)
sim_df[["sim"]] <- map(1:nrow(sim_df), large_simmer,
                       .progress = list(clear = FALSE,
                                        format = paste("{cli::pb_bar}",
                                                       "{cli::pb_percent}",
                                                       "[{cli::pb_elapsed}] |",
                                                       "ETA: {cli::pb_eta}")))
t1 <- Sys.time()
difftime(t1, t0, units = "min")

# --------------*


cat("Writing output...\n")

# Takes about 1 sec:
write_rds(sim_df, sprintf("large-r_fly-%02i.rds", curr_idx), compress = "gz")



cat("\nFINISHED!!\n\n")

EOF
