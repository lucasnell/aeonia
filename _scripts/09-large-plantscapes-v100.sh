#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-4
#SBATCH --cpus-per-task=25
#SBATCH --mem=25G
#SBATCH --time=20:00:00
#SBATCH --job-name=large-plantscapes-v100
#SBATCH --output=logs/large-plantscapes-v100-%a.out
#SBATCH --error=logs/large-plantscapes-v100-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL


#
# Larger landscape simulations with w = 1
#

# I first moved this script over to bioHPC using the following:
#
# cd ~/GitHub/Cornell/aeonia/_scripts
# scp 09-large-plantscapes-v100.sh \
#     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/09-large-plantscapes-v100/
#
# This was run on BioHPC in a non-interactive batch job started with the following:
#
# cd /home2/lan68/09-large-plantscapes-v100/
# mkdir -p logs
# sbatch 09-large-plantscapes-v100.sh
#
# Then, when the jobs are done (assuming you're back on your desktop in
# directory `~/GitHub/Cornell/aeonia/_scripts`):
#
# export RDS_DIR="/home2/lan68/09-large-plantscapes-v100"
# scp lan68@cbsugreischar.biohpc.cornell.edu:${RDS_DIR}/large-plantscapes-v100*.rds \
#     ./interm-data/
#
#


Rscript - << EOF

source("../03-large-preamble.R")

# --------------*
# Create dataframe of possible landscapes and parameter values:


sim_df <- sim_df |>
        # add row for without *Pseudomonas*
        add_row(n_pseudo = 0, wt_vp = NA, wt_pp = NA,
                landscape = list(array(0L, c(100L, 100L, 100L))) # << different from w = 0.2 sims
        ) |>
        # add other parameters:
        crossing(crossing(type = c("low", "high"),
                          # outbreaks = c("small", "big"),  # << different from w = 0.2 sims
                          # sd_N = c(0, 50),  # << different from w = 0.2 sims
                          virus_attract = c(1, 100),
                          pseudo_repel = c(1, 5))) |>
        # placeholder for simulation output:
        mutate(sim = rep(list(NA), n()))



# --------------*
# Adjust landscapes to start with 5 randomly located virus-infected plants:
#
# this whole section differs from w = 0.2 simulations

set.seed(1284259666)
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
                                            idx <- sample.int(length(tmp), 5)
                                        } else if (wt_vp < 1) {
                                            # NO overlap with Pseudomonas
                                            idx <- sample(which(tmp == 0L), 5)
                                        } else {
                                            # ALL overlap with Pseudomonas
                                            idx <- sample(which(tmp == 2L), 5)
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

# Get seed for this set of simulations:
.seed <- c(2074960264, 604827250, 777549395, 1752795709)[curr_idx]

# Verify that indices align with 'sobol_mat' object created in 'sobol-preamble.R':
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
    type <- sim_df[["type"]][[sim_df_row]]
    # outbreaks <- sim_df[["outbreaks"]][[sim_df_row]]  # << different from w = 0.2 sims
    # sd_N <- sim_df[["sd_N"]][[sim_df_row]]  # << different from w = 0.2 sims
    virus_attract <- sim_df[["virus_attract"]][[sim_df_row]]
    pseudo_repel <- sim_df[["pseudo_repel"]][[sim_df_row]]

    p_load <- 1  # << different from w = 0.2 sims
    Y0 <- 150
    N0 <- ifelse(type == "low", 150, 20)
    # >>>>>>>>>>>>>>>>>>>>>>>>
    # this chunk is uncommented and used instead of chunk above in  w = 0.2 sims
    # if (outbreaks == "big") {
    #     p_load <- 0.5
    #     Y0 <- 150
    #     N0 <- ifelse(type == "low", 150, 20)
    # } else {
    #     p_load <- 0.05
    #     Y0 <- ifelse(type == "low", 220, 300)
    #     N0 <- ifelse(type == "low", 60, 35)
    # }
    # <<<<<<<<<<<<<<<<<<<<<<<<

    zeta <- ifelse(type == "low", 1, 0.1)

    args <- list(landscape = landscape,
                 sd_N = 0,  # << different from w = 0.2 sims
                 virus_attract = virus_attract,
                 pseudo_repel = pseudo_repel,
                 Y0 = Y0,
                 N0 = N0,
                 zeta = zeta,
                 p_load_alate = p_load,
                 p_load_plant = p_load,
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
write_rds(sim_df, sprintf("large-plantscapes-v100-%02i.rds", curr_idx), compress = "gz")



cat("\nFINISHED!!\n\n")

EOF
