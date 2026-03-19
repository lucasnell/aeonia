#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=100G
#SBATCH --time=20:00:00
#SBATCH --job-name=large-plantscapes-par-vals
#SBATCH --output=large-plantscapes-par-vals.out
#SBATCH --error=large-plantscapes-par-vals.err
#SBATCH --mail-user=lan68@cornell.edu
#SBATCH --mail-type=END,FAIL


#'
#' Find par values for larger landscape simulations
#'

#' I first moved this script over to bioHPC using the following:
#'
#' cd ~/GitHub/Cornell/aeonia/_scripts
#' scp 04-large-plantscapes-par-vals.sh \
#'     lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/
#'
#' This was then run on BioHPC in a non-interactive job started with the following:
#'
#' cd /home2/lan68/
#' sbatch 04-large-plantscapes-par-vals.sh
#'
#' Then, when the job is done (assuming you're still in `~/GitHub/Cornell/aeonia/_scripts`):
#'
#' scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/large-plantscapes-par-vals.rds \
#'     ./interm-data/
#'
#'



Rscript - << EOF

source("03-large-preamble.R")


large_simmer <- function(landscape, sd_N, virus_attract, pseudo_repel,
                         Y0, N0, zeta, p_load) {

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



landscape1 <- sim_df |>
    filter(n_pseudo == 5000, wt_vp == 1e-6, wt_pp == 1,
           # These do not affect landscape:
           type == "low", sd_N == 0, virus_attract == 1,
           pseudo_repel == 1) |>
    getElement("landscape") |> getElement(1)
landscape0 <- array(c(1L, rep(0L, 99999L)), c(100L, 100L, 1L))



one_test <- function(p_load, Y0, N0, zeta) {

    n_inf1 <- large_simmer(landscape = landscape1,
                          sd_N = 0, virus_attract = 1, pseudo_repel = 1,
                          Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
        getElement("outbreak_size")
    n_inf0 <- large_simmer(landscape = landscape0,
                          sd_N = 0, virus_attract = 1, pseudo_repel = 1,
                          Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
        getElement("outbreak_size")

    tibble(n_pseudo = c(5000L, 0L),
           p_load = .env\$p_load, Y0 = .env\$Y0, N0 = .env\$N0, zeta = .env\$zeta,
           outbreak_size = c(mean(n_inf1[n_inf1 > 1]), mean(n_inf0[n_inf0 > 1])),
           p_emerge = c(mean(n_inf1 > 1), mean(n_inf0 > 1)))

}

# Takes ~3 hrs 50 min
test_sims <- crossing(p_load = round(seq(0.05, 0.2, 0.025), 3),
                      Y0 = round(seq(90, 110, 5)),
                      N0 = c(9:11, 100, 110, 120),
                      zeta = c(0.08, 0.1, 0.12, 0.9, 0.95, 1)) |>
    # Filter for two scenarios: Pseudomonas promotes or inhibits viruses, resp.:
    filter((zeta < 0.2 & N0 < 20) | (zeta > 0.8 & N0 > 90)) |>
    pmap(one_test, .progress = list(clear = FALSE,
                                    format = paste("{cli::pb_bar}",
                                                   "{cli::pb_percent}",
                                                   "[{cli::pb_elapsed}] |",
                                                   "ETA: {cli::pb_eta}"))) |>
    list_rbind()

write_rds(test_sims, "large-plantscapes-par-vals.rds", compress = "gz")

# Then run
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/large-plantscapes-par-vals.rds ~/GitHub/Cornell/aeonia/_scripts/interm-data/


# Then run script 05-large-plantscapes-par-vals.R for analysis


EOF
