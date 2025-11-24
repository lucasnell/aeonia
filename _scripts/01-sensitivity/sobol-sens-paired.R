
#'
#' Small-scale sensitivity via Sobol indices
#'
#'

#' I first moved the sensitivity preamble over to bioHPC using the following:
#'
#' scp ~/GitHub/Cornell/aeonia/_scripts/01-sensitivity/sobol-preamble.R \
#'     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/sobol-paired/
#'
#' This was then run on BioHPC in a non-interactive job started with the following:
#'
#' cd /home2/lan68/sobol-paired
#'
#'
#' # srun -n 1 -c 50 -N 1 --mem=80G --time=1-20:00:00 \
#' #   --job-name="sobol-paired" --pty bash -l
#'
#'
#' Then I started R using:
#' R --vanilla
#'
#' Then I ran the following R code interactively:
#'




slurm_nt <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (slurm_nt == "") stop("Variable SLURM_CPUS_PER_TASK not found!")
n_threads <- suppressWarnings(as.integer(slurm_nt))
if (is.na(n_threads)) stop("Argument must be an integer")
# For purrr progress bar:
.prog_args <- list(clear = FALSE,
                   format = paste("{cli::pb_bar}",
                                  "{cli::pb_percent}",
                                  "[{cli::pb_elapsed}] |",
                                  "ETA: {cli::pb_eta}"))


# Set threads for simulations:
options("mc.cores" = n_threads)



suppressPackageStartupMessages({
    library(tidyverse)
    library(aeonia)
})


source("sobol-preamble.R")


# Takes ~12 min with 10 parameters, N=2^10, and 50 threads:
t0 <- Sys.time()
sobol_sims <- all_sobol_sims(.prog_args)
write_rds(sobol_sims, "sobol-sims-paired.rds", compress = "gz")
t1 <- Sys.time()
t1 - t0


