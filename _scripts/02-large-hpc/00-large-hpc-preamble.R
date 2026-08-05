#'
#' Preamble for large landscape simulations
#'

#' I moved this script over to bioHPC using the following:
#'
#' cd ~/GitHub/Cornell/aeonia/_scripts \
#'     && scp ./02-large-hpc/00-large-hpc-preamble.R lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/
#'


# Check if running locally or on cluster:
if (grepl("biohpc.cornell.edu$", Sys.info()[["nodename"]]) &&
    Sys.info()[["user"]] == "lan68") {
    .libPaths("/home/lan68/R/x86_64-pc-linux-gnu-library/4.6")
    n_threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
    file_dir <- "/home2/lan68"
} else {
    n_threads <- max(1L, parallel::detectCores()-2L)
    file_dir <- "_scripts/interm-data"
}
if (is.na(n_threads)) stop("n_threads cannot be NA")
options("mc.cores" = n_threads)
options("readr.num_threads" = n_threads)

suppressPackageStartupMessages({
    library(tidyverse)
    library(aeonia)
})

landscape_rds <- paste0(file_dir, "/large-landscapes.rds")



# Takes ~40 sec with 6 threads
if (!file.exists(landscape_rds)) {

    set.seed(727577311)
    sim_df <- crossing(n_pseudo = as.integer(c(1e3, 3e3, 5e3, 7e3, 9e3)),
                       wt_vp = c(1e-6, 100),
                       wt_pp = c(1, 3)) |>
        mutate(landscape = pmap(
            across(everything()),
            \(n_pseudo, wt_vp, wt_pp) {
                .virus_starts <- cbind(1, 1)
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

    write_rds(sim_df, landscape_rds, compress = "gz")

} else {

    sim_df <- read_rds(landscape_rds)

}




