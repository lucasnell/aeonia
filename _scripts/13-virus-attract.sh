#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-12
#SBATCH --cpus-per-task=20
#SBATCH --mem=60G
#SBATCH --time=04:00:00
#SBATCH --job-name=virus-attract
#SBATCH --output=logs/virus-attract-%a.out
#SBATCH --error=logs/virus-attract-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

#
# Larger landscape simulations for virus attraction
#

# I first moved this script and the preamble over to bioHPC using the following:
#
# cd ~/GitHub/Cornell/aeonia/_scripts
# scp 13-virus-attract.sh \
#     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/13-virus-attract/
#
# This was then run on BioHPC in a non-interactive batch job started with the following:
#
# cd /home2/lan68/13-virus-attract/
# mkdir -p logs
# sbatch 13-virus-attract.sh
#
# Then, when the jobs are done (assuming you're back on your desktop in
# directory `~/GitHub/Cornell/aeonia/_scripts`):
#
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/13-virus-attract/virus-attract*.rds \
#     ./interm-data/virus-attract/
#
#





Rscript - << EOF

source("../03-large-preamble.R")

# --------------*
# Create dataframe of possible landscapes and parameter values:

sim_df <- sim_df |>
        filter(wt_pp == 1) |>
        # add row for without *Pseudomonas*
        add_row(n_pseudo = 0L, wt_vp = NA, wt_pp = NA,
                landscape = list(array(c(1L, rep(0L, 9999L)), c(100L, 100L, 100L)))) |>
        # add other parameters:
        crossing(crossing(wasp_resp = c("strong", "weak"),
                          p_load = 1,
                          fly_p = c(0.1, 0.5),
                          virus_attract = c(1, 5, 100))) |>
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
.seed <- c(421624097,  2128854588, 620380686,  343417722,  1731760016,
           1362204649, 1553761294, 1692916644, 1739955373, 1364476424,
           709020402,  1766198254)[curr_idx]

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

    landscape <- sim_df[["landscape"]][[sim_df_row]]
    wasp_resp <- sim_df[["wasp_resp"]][[sim_df_row]]
    p_load <- sim_df[["p_load"]][[sim_df_row]]
    fly_p <- sim_df[["fly_p"]][[sim_df_row]]
    virus_attract <- sim_df[["virus_attract"]][[sim_df_row]]

    args <- list(landscape = landscape,
                 virus_attract = virus_attract,
                 fly_p = fly_p,
                 zeta = ifelse(wasp_resp == "weak", 0.1, 0.9),
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 sd_N = 0,
                 pseudo_repel = 1,
                 Y0 = 250,
                 N0 = 55,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = 100L,
                 summ = "all")

    return(do.call(big_plantscape, args))

}


# Took ~30 min per job, where each job processed 11 rows using 20 threads

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
write_rds(sim_df, sprintf("virus-attract-%02i.rds", curr_idx), compress = "gz")



cat("\nFINISHED!!\n\n")

EOF
