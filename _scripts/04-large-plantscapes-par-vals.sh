#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --array=1-6
#SBATCH --cpus-per-task=25
#SBATCH --mem=25G
#SBATCH --time=20:00:00
#SBATCH --job-name=large-plantscapes-par-vals
#SBATCH --output=logs/large-plantscapes-par-vals-%a.out
#SBATCH --error=logs/large-plantscapes-par-vals-%a.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL



#
# Find par values for larger landscape simulations
#

# I first moved this script over to bioHPC using the following:
#
# cd ~/GitHub/Cornell/aeonia/_scripts
# scp 03-large-preamble.R lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/
# scp 04-large-plantscapes-par-vals.sh \
#     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/04-large-plantscapes-par-vals/
#
# This was then run on BioHPC in a non-interactive job started with the following:
#
# cd /home2/lan68/04-large-plantscapes-par-vals/
# mkdir -p logs
# sbatch 04-large-plantscapes-par-vals.sh
#
# Then, when the job is done (assuming you're still in `~/GitHub/Cornell/aeonia/_scripts`):
#
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/04-large-plantscapes-par-vals/large-plantscapes-par-vals*.rds \
#     ./interm-data/
#
# Then run script 05-large-plantscapes-par-vals.R for analysis
#



Rscript - << EOF

source("../03-large-preamble.R")


large_simmer <- function(landscape, Y0, N0, zeta, p_load) {

    args <- list(landscape = landscape,
                 sd_N = 0,
                 virus_attract = 5,  # << because this maximizes p_emerge
                 pseudo_repel = 1,
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


.n_pseudo <- 7000L  # << maximizes effect of Pseudomonas when it promotes outbreaks

landscape1 <- sim_df |>
    filter(n_pseudo == .n_pseudo, wt_vp == 1e-6, wt_pp == 1) |>
    getElement("landscape") |> getElement(1)
landscape0 <- array(c(1L, rep(0L, 99999L)), c(100L, 100L, 1L))




one_test <- function(p_load, Y0, N0, zeta) {

    n_inf1 <- large_simmer(landscape = landscape1,
                          Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
        getElement("n_infected")
    n_inf0 <- large_simmer(landscape = landscape0,
                          Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
        getElement("n_infected")

    tibble(n_pseudo = c(.n_pseudo, 0L),
           p_load = .env\$p_load, Y0 = .env\$Y0, N0 = .env\$N0, zeta = .env\$zeta,
           outbreak_size = c(mean(n_inf1[n_inf1 > 1]), mean(n_inf0[n_inf0 > 1])),
           p_emerge = c(mean(n_inf1 > 1), mean(n_inf0 > 1)),
           n_infected = c(mean(n_inf1), mean(n_inf0)))

}



test_sim_df <- crossing(p_load = c(0.01, 0.05, 0.5),
                        Y0 = round(seq(50, 500, 50)),
                        N0 = 55,
                        zeta = c(0.1, 0.9))

# --------------*
# Read inputs from job manager:

curr_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")))
if (is.na(curr_idx)) stop("SLURM_ARRAY_TASK_ID must be an integer")

max_idx <- suppressWarnings(as.integer(Sys.getenv("SLURM_ARRAY_TASK_MAX")))
if (is.na(max_idx)) stop("SLURM_ARRAY_TASK_MAX must be an integer")

# Get seed for this set of simulations:
.seed <- c(769895830, 1692138374, 1774036889, 2036965360, 1974372476, 1538998712)[curr_idx]

# Verify that indices align with 'test_sim_df' object:
stopifnot(nrow(test_sim_df) %% max_idx == 0L)

# number of rows to simulate per job:
n_rows_pj <- nrow(test_sim_df) %/% max_idx

# Start and stop for rows to simulate this job:
curr_start <- (curr_idx - 1L) * n_rows_pj + 1L
curr_stop <- curr_idx * n_rows_pj


test_sim_df <- test_sim_df[curr_start:curr_stop,]

# --------------*
# Run simulations:

# Each job takes ~2 min per row using 25 threads
t0 <- Sys.time()
test_sim_df <- test_sim_df |>
    pmap(one_test, .progress = list(clear = FALSE,
                                    format = paste("{cli::pb_bar}",
                                                   "{cli::pb_percent}",
                                                   "[{cli::pb_elapsed}] |",
                                                   "ETA: {cli::pb_eta}"))) |>
    list_rbind()
t1 <- Sys.time()
difftime(t1, t0, units = "min")

write_rds(test_sim_df, sprintf("large-plantscapes-par-vals-%02i.rds", curr_idx), compress = "gz")

cat("\nFINISHED!!\n\n")

EOF
