#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-6
#SBATCH --cpus-per-task=25
#SBATCH --mem=40G
#SBATCH --time=04:00:00
#SBATCH --job-name=large-w1
#SBATCH --output=logs/large-w1-%a.out
#SBATCH --error=logs/large-w1-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL


#
# Larger landscape simulations with w = 1
#

# I first moved this script over to bioHPC using the following:
#
# cd ~/GitHub/Cornell/aeonia/_scripts
# scp 07-large-w1.sh \
#     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/07-large-w1/
#
# This was run on BioHPC in a non-interactive batch job started with the following:
#
# cd /home2/lan68/07-large-w1/
# mkdir -p logs
# sbatch 07-large-w1.sh
#
# Then, when the jobs are done (assuming you're back on your desktop in
# directory `~/GitHub/Cornell/aeonia/_scripts`):
#
# export RDS_DIR="/home2/lan68/07-large-w1"
# scp lan68@cbsugreischar.biohpc.cornell.edu:${RDS_DIR}/large-w1*.rds \
#     ./interm-data/large-w1
#
#


Rscript - << EOF

source("../03-large-preamble.R")


# --------------*
# Create dataframe of possible landscapes and parameter values:

sim_df <- sim_df |>
        # add row for without *Pseudomonas*
        add_row(n_pseudo = 0, wt_vp = NA, wt_pp = NA,
                landscape = list(array(0L, c(100L, 100L, 100L)))
        ) |>
        # add other parameters:
        crossing(crossing(wasp_resp = c("strong", "weak"),
                          virus_attract = c(1, 5, 100))) |>
        # placeholder for simulation output:
        mutate(sim = rep(list(NA), n()))


# --------------*
# Adjust landscapes to start with 5 randomly located virus-infected plants:
#
# this whole section differs from w = 0.2 simulations

# Number of virus-infected plants to start.
n_inf <- 5L

set.seed(1928347678)
sim_df[["landscape"]] <- map2(sim_df[["landscape"]], sim_df[["wt_vp"]],
                                function(land, wt_vp) {
                                    # reset to have no virus-infected plants:
                                    land[land == 1L] <- 0L
                                    land[land == 3L] <- 2L
                                    stopifnot(all(land %in% c(0L, 2L)))
                                    for (j in 1:dim(land)[3]) {
                                        tmp <- land[,,j]
                                        if (is.na(wt_vp) || wt_vp == 1) {
                                            # Totally random
                                            idx <- sample.int(length(tmp), n_inf)
                                        } else if (wt_vp < 1) {
                                            # NO overlap with Pseudomonas
                                            idx <- sample(which(tmp == 0L), n_inf)
                                        } else {
                                            # ALL overlap with Pseudomonas
                                            idx <- sample(which(tmp == 2L), n_inf)
                                        }
                                        tmp[idx] <- tmp[idx] + 1L
                                        land[,,j] <- tmp
                                    }
                                    return(land)
                                })



# --------------*
# Read inputs from job manager:

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")

max_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX")))
if (is.na(max_idx)) stop("SLURM_ARRAY_TASK_MAX must be an integer")
if (curr_idx > max_idx) stop("SLURM_ARRAY_TASK_ID cannot exceed SLURM_ARRAY_TASK_MAX")

# Get seed for this set of simulations:
.seed <- c(921782975, 1753992136, 506798594, 1245256270, 1635039610, 50641579)[curr_idx]

# Verify that indices align with 'sim_df' object:
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
    virus_attract <- sim_df[["virus_attract"]][[sim_df_row]]

    args <- list(landscape = landscape,
                 sd_N = 0,  # << different from w = 0.2 sims
                 virus_attract = virus_attract,
                 pseudo_repel = 1,  # << different from w = 0.2 sims
                 Y0 = 250,
                 N0 = 55,
                 zeta = ifelse(wasp_resp == "weak", 0.1, 0.9),
                 p_load_alate = 1,  # << different from w = 0.2 sims
                 p_load_plant = 1,  # << different from w = 0.2 sims
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = 100L,
                 summ = "all",
                 fly_p = 0.5,  # << different from w = 0.2 sims
                 w = 1)  # << different from w = 0.2 sims

    return(do.call(big_plantscape, args))

}






# Took ~50 min per job, where each job processed 42 rows using 25 threads

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
write_rds(sim_df, sprintf("large-w1-%02i.rds", curr_idx), compress = "gz")



cat("\nFINISHED!!\n\n")

EOF
