#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-4
#SBATCH --cpus-per-task=25
#SBATCH --mem=25G
#SBATCH --time=20:00:00
#SBATCH --job-name=large-plantscapes
#SBATCH --output=logs/large-plantscapes-%a.out
#SBATCH --error=logs/large-plantscapes-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

#
# Larger landscape simulations
#

# I first moved this script and the preamble over to bioHPC using the following:
#
# cd ~/GitHub/Cornell/aeonia/_scripts
# scp 03-large-preamble.R \
#     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/
# scp 06-large-plantscapes.sh \
#     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/06-large-plantscapes/
#
# This was then run on BioHPC in a non-interactive batch job started with the following:
#
# cd /home2/lan68/06-large-plantscapes/
# mkdir -p logs
# sbatch 06-large-plantscapes.sh
#
# Then, when the jobs are done (assuming you're back on your desktop in
# directory `~/GitHub/Cornell/aeonia/_scripts`):
#
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/06-large-plantscapes/large-plantscapes*.rds \
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
                landscape = list(array(c(1L, rep(0L, 99999L)), c(100L, 100L, 1L)))) |>
        # add other parameters:
        crossing(crossing(type = c("low", "high"),
                          outbreaks = c("small", "big"),
                          sd_N = c(0, 50),
                          virus_attract = c(1, 5),
                          pseudo_repel = c(1, 5))) |>
        # placeholder for simulation output:
        mutate(sim = rep(list(NA), n()))



# --------------*
# Read inputs from job manager:

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")

max_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX")))
if (is.na(max_idx)) stop("SLURM_ARRAY_TASK_MAX must be an integer")

# Get seed for this set of simulations:
.seed <- c(1409962309, 969620847, 693008744, 828602012)[curr_idx]

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
    outbreaks <- sim_df[["outbreaks"]][[sim_df_row]]
    sd_N <- sim_df[["sd_N"]][[sim_df_row]]
    virus_attract <- sim_df[["virus_attract"]][[sim_df_row]]
    pseudo_repel <- sim_df[["pseudo_repel"]][[sim_df_row]]

    if (outbreaks == "big") {
        p_load <- 0.5
        Y0 <- 150
        N0 <- ifelse(type == "low", 150, 20)
    } else {
        p_load <- 0.05
        Y0 <- ifelse(type == "low", 220, 300)
        N0 <- ifelse(type == "low", 60, 35)
    }

    zeta <- ifelse(type == "low", 1, 0.1)

    args <- list(landscape = landscape,
                 sd_N = sd_N,
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
                 summ = "all")

    return(do.call(big_plantscape, args))

}





# Took ~3 hours per job, where each job processed 168 rows using 25 threads

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
write_rds(sim_df, sprintf("large-plantscapes-%02i.rds", curr_idx), compress = "gz")



cat("\nFINISHED!!\n\n")

EOF
