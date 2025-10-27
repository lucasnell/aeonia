
#'
#' Testing for whether parameter combos that seemed to result in a negative
#' effect of Pseudomonas continue to show this pattern with greater numbers
#' of simulations.
#'
#' Start job on cluster:
#'
#' cd /home2/lan68/sobol-paired
#' srun -n 1 -c 100 -N 1 --mem=80G --time=20:00:00 --job-name="sobol-paired-test" --pty bash -l
#'
#' Then:
#'
#' R --vanilla


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

#' Output from `sobol-sens-paired.R`:
sobol_sims <- read_rds("sobol-sims-paired.rds")[[2L]]



# Summarize each set of simulations for differences between with and
# without Pseudo:
# Takes ~13 sec  (multithreading doesn't help)
diff_sobol_summs <- map(sobol_sims, \(sim_set) {
    out_df <- sim_set[[1]][1, c(names(struct_pars), names(vary_pars))]
    out_df[["n_pseudo"]] <- NULL
    for (y in yvars) {
        y_p <- sim_set$pseudo[[y]]
        y_np <- sim_set$no_pseudo[[y]]
        if (y != "infect_time" && any(is.na(c(y_p, y_np)))) {
            stop(y, " has NA values")
        }
        # using rounding bc weird, very small values show up here (~1e-16)
        # if I don't. Plus, I know that the difference can't have more than
        # 2 decimal digits bc n_sims = 100.
        out_df[[y]] <- round(mean(y_p, na.rm = TRUE) - mean(y_np, na.rm = TRUE), 2)
    }
    # now do prob. outbreak happened:
    out_df[["p_outbreak"]] <- round(mean(y_p > 1) - mean(y_np > 1), 2)
    # and outbreak size when there was one:
    out_df[["outbreak_size2"]] <- round(mean(y_p[y_p > 1]) - mean(y_np[y_np > 1]), 2)

    return(out_df)

}) |>
    list_rbind()



sum(diff_sobol_summs$outbreak_size > 0)
# [1] 1959


# Test the top 100 parameter combinations:
diff_test_df <- diff_sobol_summs |>
    arrange(desc(outbreak_size)) |>
    select(all_of(names(vary_pars)), outbreak_size) |>
    slice_head(n = 100)



all_diff_test_sims <- function(.prog_args = FALSE) {

    args0 <- struct_pars[2,] |> as.list()

    sim_outs <- map(1:nrow(diff_test_df), \(j) {

        args_j <- args0
        for (n in names(vary_pars)) {
            args_j[[n]] <- unname(diff_test_df[[n]][[j]])
        }
        args_j[["n_sims"]] <- 1000L

        sim <- do.call(one_combo, args_j)

        args_j[["n_pseudo"]] <- 0L
        sim0 <- do.call(one_combo, args_j)

        out <- list(pseudo = sim, no_pseudo = sim0)

        return(out)

    }, .progress = .prog_args)

    return(sim_outs)
}



# Takes ~ min with 100 combos of 9 parameters and 100 threads:
t0 <- Sys.time()
sobol_test_sims <- all_diff_test_sims(.prog_args)
write_rds(sobol_test_sims, "sobol-test-sims-paired.rds", compress = "gz")
t1 <- Sys.time()
t1 - t0

