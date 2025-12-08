
if (! "tidyverse" %in% .packages()) suppressPackageStartupMessages(library(tidyverse))
if (! "aeonia" %in% .packages()) suppressPackageStartupMessages(library(aeonia))

suppressPackageStartupMessages({
    library(sensobol)
    library(lhs)
})



#'
#'
#'
#'
#' `spat_config` indicates the spatial configuration of *Pseudomonas* and
#' the values indicate the following (v = virus, p = *Pseudomonas*, B = both):
#'
#' - `-1`: random configuration (including on initial virus plant)
#' - `0`: random configuration (never starting on initial virus plant)
#' - `1`:
#'     |---|---|---|
#'     | v | - | p |
#'     | - | p | - |
#'     | p | - | - |
#'     |---|---|---|
#' - `2`:
#'     |---|---|---|
#'     | v | p | - |
#'     | p | p | - |
#'     | - | - | - |
#'     |---|---|---|
#' - `3`:
#'     |---|---|---|
#'     | v | - | - |
#'     | - | - | p |
#'     | - | p | p |
#'     |---|---|---|
#' - `4`:
#'     |---|---|---|
#'     | B | p | - |
#'     | p | - | - |
#'     | - | - | - |
#'     |---|---|---|
#'
#'


vary_pars <- list(Y0 = c(1, 9),
                  mean_N = c(10, 100),
                  sd_N = c(0, 50),
                  K = 12500 * c(0.5, 2),
                  virus_attract = c(1, 10),
                  pseudo_repel = c(1, 10),
                  pseudo_surv = c(0.85, 1), # <= 0.8344661 results in carrying capacity of ~0
                  zeta = c(0, 1))

# Response variables
yvars <- col_namer("all", FALSE, FALSE) |> (\(x) x[x != "rep"])()


N <- 2^12

# Takes just a second or two
set.seed(641272456)
sobol_mat <- (\(i) {
    mat <- sobol_matrices(N = N, params = names(vary_pars))
    for (n in names(vary_pars)) {
        p1 <- min(vary_pars[[n]])
        p2 <- max(vary_pars[[n]])
        if (! n %in% c("spat_config", "alate_dens")) {
            mat[,n] <- p1 + mat[,n] * (p2 - p1)
        } else {
            mat[,n] <- qinteger(mat[,n], p1, p2)
        }
    }
    return(mat)
})()


#' Should be TRUE
identical(sort(names(vary_pars)), sort(colnames(sobol_mat)))



# Note: `.empir_N` forces N0 to follow `mean_N` and `sd_N` exactly
#        (or nearly so, given numerical inexactness in R)
one_combo <- function(n_pseudo, alate_dens,
                      Y0, mean_N, sd_N, K, virus_attract, pseudo_repel, pseudo_surv,
                      zeta, n_sims,
                      spat_config = -1L,
                      .empir_N = TRUE,
                      ...) {

    stopifnot(alate_dens %in% 0:1)


    fly_p <- 0.05
    epsilon <- 1
    p_load_alate <- 0.5
    p_load_plant <- 0.5

    n_x <- 3L
    n_y <- 3L
    n_plants <- n_x * n_y

    # Convert from mean and sd of lognormal distribution to parameters
    # to use for lognormal (mean and sd of underlying normal distribution):
    mu_N <- log(mean_N^2 / sqrt(mean_N^2 + sd_N^2))
    sigma_N <- sqrt(log(1 + sd_N^2 / mean_N^2))

    if (alate_dens == 0L) {
        alate_slope <- 0
        alate_max  <- 0.075
    } else {
        alate_slope <- pop_info$alate_slope
        alate_max <- 1
    }


    land <- array(0L, c(n_x, n_y, n_sims))
    N0 <- array(0.0, c(n_x, n_y, n_sims))
    for (i in 1:n_sims) {
        land[1,1,i] <- 1L
        if (n_pseudo > 0) {
            if (spat_config == -1L) {
                k <- sample.int(n_plants, n_pseudo)
                x <- k - n_x * ((k-1L) %/% n_x)
                y <- (k-1L) %/% n_x + 1L
                # Deal with situation where Pseudomonas on same plant as virus:
                if (any(x == 1 & y == 1)) {
                    idx <- which(x == 1 & y == 1)
                    x <- x[-idx]
                    y <- y[-idx]
                    land[1,1,i] <- 3L
                }
            } else if (spat_config == 0L) {
                k <- sample.int(n_plants - 1L, n_pseudo)
                x <- k - n_x * (k %/% n_x) + 1L
                y <- k %/% n_x + 1L
            } else if (spat_config == 1L) {
                x <- 3:1
                y <- 1:3
            } else if (spat_config == 2L) {
                x <- c(1, 2, 2)
                y <- c(2, 1, 2)
            } else if (spat_config == 3L) {
                x <- c(2, 3, 3)
                y <- c(3, 2, 3)
            } else if (spat_config == 4L) {
                land[1,1,i] <- 3L
                x <- c(1, 2)
                y <- c(2, 1)
            }
            land[cbind(x,y,i)] <- 2L
        }
        if (.empir_N) {
            N0[,,i] <- rnorm(9, 0, 1) |>
                (\(x) ((x - mean(x)) / sd(x)) * sigma_N + mu_N)() |>
                exp()
        } else {
            N0[,,i] <- rlnorm(n_plants, mu_N, sigma_N)
        }

    }

    insect_args <- list(K = K, pseudo_surv = pseudo_surv, fly_p = fly_p,
                        zeta = zeta, alate_slope = alate_slope,
                        alate_max = alate_max)
    plant_args <- list(landscapes = land,
                       max_t = 100,
                       N0 = N0,
                       W0 = array(0.0, c(n_x, n_y, n_sims)),
                       Y0 = Y0,
                       virus_attract = virus_attract,
                       pseudo_repel = pseudo_repel,
                       epsilon = epsilon,
                       infect_time_n = 5,
                       p_load_alate = p_load_alate,
                       p_load_plant = p_load_plant,
                       radius = 1,
                       infect_stop = FALSE,
                       summ = "all")
    other_args <- list(...)
    if (length(other_args) > 0) {
        stopifnot(!is.null(names(other_args)) && !any(names(other_args) == ""))
        if (!all(names(other_args) %in% c(names(formals(sim_plantscape)),
                                          names(formals(make_insect_ptr))))) {
            all_names <- c(names(formals(sim_plantscape)),
                           names(formals(make_insect_ptr)))
            print(names(other_args)[!names(other_args) %in% all_names])
            stop("\nNot all names in ... match args in sim_plantscape ",
                 "or make_insect_ptr")
        }
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

    out <- do.call(sim_plantscape, plant_args)

    # Now add arguments to the output dataframe:
    in_args <- as.list(match.call(expand.dots = TRUE))[-1]  # -1 removes fxn name
    # Removing boring parameters:
    boring_args <- c("n_sims", "summ", "max_t", "infect_stop",
                     "out_pseudo", "show_progress", "n_threads")
    in_args <- in_args[-which(names(in_args) %in% boring_args)]
    # Now add to output:
    for (n in names(in_args)) {
        if (length(in_args[[n]]) == 1) out[[n]] <- in_args[[n]]
        else out[[n]] <- list(in_args[[n]])
    }
    # Make these args come first:
    out <- out |> select(all_of(names(in_args)), everything())

    return(out)

}


