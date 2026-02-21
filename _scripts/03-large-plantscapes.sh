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

#' I first moved this script over to bioHPC using the following:
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

source("03-large-preamble.R")


# --------------*
# Create dataframe of possible landscapes and parameter values:


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
