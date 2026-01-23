#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --time=20:00:00
#SBATCH --job-name=large-plantscapes
#SBATCH --output=large-plantscapes.out
#SBATCH --error=large-plantscapes.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL

#'
#' Larger landscape simulations
#'

#' I first moved the sensitivity preamble and this script over to bioHPC
#' using the following:
#'
#' cd ~/GitHub/Cornell/aeonia/_scripts \
#'     && scp 03-large-plantscapes.sh lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/
#'
#' This was then run on BioHPC in a non-interactive job started with the following:
#'
#' cd /home2/lan68/ \
#'     && sbatch 03-large-plantscapes.sh
#'
#' Then, when the job is done (assuming you're still in `~/GitHub/Cornell/aeonia/_scripts`):
#'
#' scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/large-plantscapes.rds ./interm-data/
#'
#'





Rscript - << EOF

.libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.5")

##> NOTE: Running multiple simulations at a time (and multithreading
##> per simulation) via the future package took longer than just using
##> 50 threads per simulation.

suppressPackageStartupMessages({
    library(tidyverse)
    library(aeonia)
})

n_threads <- suppressWarnings(as.integer("${SLURM_CPUS_PER_TASK}"))
if (is.na(n_threads)) stop("SLURM_CPUS_PER_TASK must be an integer")
options("mc.cores" = n_threads)


# --------------*
# Create dataframe of possible landscapes and parameter values:

cat("Creating simulation inputs...\n")

# Note: creating landscapes first like this bc I want to use the exact same 13
# landscapes over and over, to reduce dependency on exact landscape
# configurations.

set.seed(727577311)
sim_df <- crossing(n_pseudo = c(1e3, 3e3, 5e3, 7e3, 9e3),
                   wt_vp = c(1e-6, 100),
                   wt_pp = c(1, 3)) |>
    mutate(landscape = pmap(
        across(everything()),
        \(n_pseudo, wt_vp, wt_pp) {
            if (wt_vp < 1) {
                .pseudo_starts <- cbind(100, 100)
            } else .pseudo_starts <- cbind(1, 1)
            sim_plant_types(n_x = 100,
                            n_y = 100,
                            n_lands = 1,
                            n_virus = 1,
                            n_pseudo = n_pseudo,
                            wt_vp = wt_vp,
                            wt_pp = wt_pp,
                            virus_starts = cbind(1, 1),
                            pseudo_starts = .pseudo_starts)
        })) |>
    # add row for without *Pseudomonas*
    add_row(n_pseudo = 0, wt_vp = NA, wt_pp = NA,
            landscape = list(array(c(1L, rep(0L, 99999L)), c(100L, 100L, 1L)))) |>
    # add other parameters:
    crossing(crossing(type = c("low", "high"),
                      sd_N = c(0, 50),
                      virus_attract = c(1, 5),
                      pseudo_repel = c(1, 5))) |>
    # placeholder for simulation output:
    mutate(sim = rep(list(NA), n()))



large_simmer <- function(sim_df_row) {

    landscape <- sim_df[["landscape"]][[sim_df_row]]
    type <- sim_df[["type"]][[sim_df_row]]
    sd_N <- sim_df[["sd_N"]][[sim_df_row]]
    virus_attract <- sim_df[["virus_attract"]][[sim_df_row]]
    pseudo_repel <- sim_df[["pseudo_repel"]][[sim_df_row]]

    args <- list(landscape = landscape,
                 sd_N = sd_N,
                 virus_attract = virus_attract,
                 pseudo_repel = pseudo_repel,
                 Y0 = 100,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = 100,
                 summ = "all")

    if (type == "low") {
        args <- list_assign(args,
                            N0 = 110,
                            zeta = 1)
    } else {
        args <- list_assign(args,
                            N0 = 10,
                            zeta = 0.1)
    }

    return(do.call(big_plantscape, args))

}




# Took ~ min with 50 threads and no nested parallelism

cat("Starting simulations...\n")
t0 <- Sys.time()


set.seed(1772344300)
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
write_rds(sim_df, "large-plantscapes.rds", compress = "gz")


cat("\nFINISHED!!\n\n")

EOF
