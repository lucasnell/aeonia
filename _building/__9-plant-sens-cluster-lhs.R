#!/usr/bin/env -S Rscript --vanilla


if (!interactive()) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) != 1) stop("Must have exactly 1 argument")
    n_threads <- suppressWarnings(as.integer(args[[1]]))
    if (is.na(n_threads)) stop("Argument must be an integer")
    # For purrr progress bar:
    .prog_args <- FALSE
} else {
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
}



suppressPackageStartupMessages({
    # make sure these come before tidyverse to avoid masking dplyr::select:
    library(lhs)
    library(ppcor)
    # ^^
    library(tidyverse)
    library(aeonia)
    library(future.apply)
})


# Set threads for simulations:
options("mc.cores" = n_threads)

# And for using future.apply functions:
plan(multisession, workers = options()[["mc.cores"]])





#'
#' Parameters that aren't especially interesting or relevant:
#'   - fly_p (set to 0.05 below)
#'   - epsilon (set to 1 below)
#'   - sigma_Y (set to 1 below)
#'   - sigma_N (set to 1 below)


#'
#' Parameters that I need to structure simulations within:
#'
struct_pars <- crossing(mu_Y = c(-Inf, -2, 0, 2),
                        mu_N = c(0, 2, 4),
                        n_pseudo = c(0, 2, 4) |> as.integer(),
                        K = 12500 * c(0.75, 1, 1.5))

#'
#' Min and max values for each parameter to vary:
#'
vary_pars <- list(B = c(0, 0.15), # >= 0.1655339 results in carrying capacity of ~0
                  alpha = c(0, 5),
                  beta = c(-5, 0),
                  wasp_disp_m0 = 0.3 * c(0, 5),
                  wasp_disp_m1 = 0.349 * c(0, 5))

par_names <- c(names(struct_pars), names(vary_pars))


# Takes ~ 11 sec with 6 threads:
lhs_df <- struct_pars |>
    mutate(other_pars = future_lapply(1:n(), \(i) {
        maximinLHS(1000, length(vary_pars)) |>
            as.data.frame() |>
            set_names(names(vary_pars)) |>
            as_tibble() |>
            mutate(across(everything(), \(x) {
                p <- vary_pars[[cur_column()]]
                p[1] + x * (p[2] - p[1])
            })) |>
            mutate(n_sims = 100)
    }, future.seed = 1980974943, future.packages = c("tidyverse", "lhs")))

# Note: on cluster, you can ignore warnings of the form
# "'package:XXXX' may not be available when loading"


one_combo <- function(mu_Y, mu_N, n_pseudo, K,
                      B, alpha, beta, wasp_disp_m0, wasp_disp_m1,
                      n_sims, ...) {

    fly_p <- 0.05
    epsilon <- 1
    sigma_Y <- 1
    sigma_N <- 1

    n_x <- 3L
    n_y <- 3L
    n_plants <- n_x * n_y

    land <- array(0L, c(n_x, n_y, n_sims))
    N0 <- array(0.0, c(n_x, n_y, n_sims))
    Y0 <- array(0.0, c(n_x, n_y, n_sims))
    for (i in 1:n_sims) {
        land[1,1,i] <- 1L
        if (n_pseudo > 0) {
            k <- sample.int(n_plants - 1L, n_pseudo)
            x <- k - n_x * (k %/% n_x) + 1L
            y <- k %/% n_x + 1L
            land[cbind(x,y,i)] <- 2L
        }
        N0[,,i] <- rlnorm(n_plants, mu_N, sigma_N)
        Y0[,,i] <- rlnorm(n_plants, mu_Y, sigma_Y)
    }

    insect_args <- list(K = K, B = B, fly_p = fly_p,
                        wasp_disp_m0 = wasp_disp_m0,
                        wasp_disp_m1 = wasp_disp_m1)
    plant_args <- list(landscapes = land,
                       max_t = 100,
                       N0 = N0,
                       W0 = array(0.0, c(n_x, n_y, n_sims)),
                       Y0 = Y0,
                       alpha = alpha,
                       beta = beta,
                       epsilon = epsilon,
                       infect_time_n = 5,
                       delta_a = 0.5,
                       delta_p = 0.5,
                       radius = 1,
                       infect_stop = FALSE,
                       summ = "all")
    other_args <- list(...)
    if (length(other_args) > 0) {
        stopifnot(!is.null(names(other_args)) && !any(names(other_args) == ""))
        stopifnot(all(names(other_args) %in% c(names(formals(sim_plantscape)),
                                               names(formals(make_insect_ptr)))))
        not_allowed <- c("landscapes", "N0", "W0", "Y0")
        if (any(names(other_args) %in% not_allowed)) {
            not_allowed <- names(other_args)[names(other_args) %in% not_allowed]
            stop("The following are not allowed in `make_arg_list`: ",
                 paste(not_allowed, collapse = ", "))
        }

        nm_insect_args <- names(other_args)[names(other_args) %in%
                                                names(formals(make_insect_ptr))]
        for (n in nm_insect_args) insect_args[[n]] <- other_args[[n]]

        nm_plant_args <- names(other_args)[names(other_args) %in%
                                               names(formals(sim_plantscape))]
        for (n in nm_plant_args) {
            if (n %in% c("N0", "W0", "Y0", "landscapes")) {
                plant_args[[n]] <- array(other_args[[n]], c(n_x, n_y, n_sims))
            } else plant_args[[n]] <- other_args[[n]]
        }
    }

    plant_args[["insect_ptr"]] <- do.call(make_insect_ptr, insect_args)

    out <- do.call(sim_plantscape, plant_args) |>
        mutate(mu_N = mu_N, mu_Y = mu_Y, n_pseudo = n_pseudo,
               B = B, K = K, alpha = alpha, beta = beta,
               wasp_disp_m0 = wasp_disp_m0, wasp_disp_m1 = wasp_disp_m1)

    return(out)

}




# Takes ~24 min with 108 combos of 1000 reps and 50 threads:
set.seed(891393509)
lhs_sims <- lhs_df |>
    pmap(\(mu_Y, mu_N, n_pseudo, K, other_pars) {
        other_pars |>
            pmap(one_combo, mu_Y = mu_Y, mu_N = mu_N,
                 n_pseudo = n_pseudo, K = K) |>
            list_rbind() |>
            select(all_of(par_names), everything())
    }, .progress = .prog_args)
write_rds(lhs_sims, paste0(gsub("home/", "home2/", getwd()), "/lhs-sims.rds"),
          compress = "gz")




