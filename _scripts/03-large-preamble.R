#'
#' Preamble for larger landscape simulations
#'

#' I moved this script over to bioHPC using the following:
#'
#' cd ~/GitHub/Cornell/aeonia/_scripts \
#'     && scp 03-large-preamble.R lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/
#'


# Check if running locally or on cluster:
if (grepl("biohpc.cornell.edu$", Sys.info()[["nodename"]]) &&
    Sys.info()[["user"]] == "lan68") {
    .libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.5")
    n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
} else {
    n_threads <- max(1L, parallel::detectCores()-2L)
}
if (is.na(n_threads)) stop("n_threads cannot be NA")
options("mc.cores" = n_threads)

suppressPackageStartupMessages({
    library(tidyverse)
    library(aeonia)
})



# Takes ~40 sec with 6 threads
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
                            n_lands = 100,
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

