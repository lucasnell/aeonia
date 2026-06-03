#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-32
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=04:00:00
#SBATCH --job-name=large-interior
#SBATCH --output=logs/large-interior-%a.out
#SBATCH --error=logs/large-interior-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

#
# Larger landscape simulations with interior virus starting locations
#

# I first moved this script and the preamble over to bioHPC using the following:
#
# cd ~/GitHub/Cornell/aeonia/_scripts
# scp 07b-large-interior.sh \
#     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/07b-large-interior/
#
# This was then run on BioHPC in a non-interactive batch job started with the following:
#
# cd /home2/lan68/07b-large-interior/
# mkdir -p logs
# sbatch 07b-large-interior.sh
#
# Then, when the jobs are done (assuming you're back on your desktop in
# directory `~/GitHub/Cornell/aeonia/_scripts`):
#
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/07b-large-interior/large-interior*.rds \
#     ./interm-data/large-interior/
#
#





Rscript - << EOF

source("../03-large-preamble.R")

# --------------*
# Create new dataframe of landscapes with  viruses starting at x=50, y=50

set.seed(48353817)
sim_df <- crossing(n_pseudo = as.integer(c(1e3, 3e3, 5e3, 7e3, 9e3)),
                   wt_vp = c(1e-6, 100),
                   wt_pp = c(1, 3)) |>
    mutate(landscape = pmap(
        across(everything()),
        \(n_pseudo, wt_vp, wt_pp) {
            .virus_starts <- cbind(50, 50)
            if (wt_vp < 1) {
                .pseudo_starts <- cbind(100, 100)
            } else .pseudo_starts <- .virus_starts
            sim_plant_types(n_x = 100,
                            n_y = 100,
                            n_lands = 100,
                            n_virus = 1,
                            n_pseudo = n_pseudo,
                            wt_vp = wt_vp,
                            wt_pp = wt_pp,
                            virus_starts = .virus_starts,
                            pseudo_starts = .pseudo_starts)
        }))


# --------------*
# Create dataframe of possible landscapes and parameter values:

sim_df <- sim_df |>
        # add row for without *Pseudomonas* (note landscape diff. from corner sims)
        add_row(n_pseudo = 0L, wt_vp = NA, wt_pp = NA,
                landscape = list(array(c(rep(0L, 4949L), 1L, rep(0L, 5050L)), rep(100, 3)))) |>
        # add other parameters:
        crossing(crossing(wasp_resp = c("strong", "weak"),
                          p_load = c(0.5, 0.05),
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
if (curr_idx > max_idx) stop("SLURM_ARRAY_TASK_ID cannot exceed SLURM_ARRAY_TASK_MAX")

# Get seed for this set of simulations:
.seed <- c(1930128157, 1662648893, 1607036552, 1072583135, 505412883,
           631177519,  716115502,  1996740403, 439275693,  1788179576,
           923580245,  2054256536, 835671563,  644843138,  1880075451,
           1118415247, 1772031408, 1011779535, 1762018110, 73177731,
           186331081,  1202052869, 715261268,  267763020,  1812005114,
           1590545119, 439528765,  2088582417, 2011625946, 1527599829,
           47903130, 1100864230)[curr_idx]

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
    sd_N <- sim_df[["sd_N"]][[sim_df_row]]
    virus_attract <- sim_df[["virus_attract"]][[sim_df_row]]
    pseudo_repel <- sim_df[["pseudo_repel"]][[sim_df_row]]

    args <- list(landscape = landscape,
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





# Took ~1.5 hours per job, where each job processed 63 rows using 8 threads

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
write_rds(sim_df, sprintf("large-interior-%02i.rds", curr_idx), compress = "gz")



cat("\nFINISHED!!\n\n")

EOF
