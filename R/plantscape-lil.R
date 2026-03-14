


# make a little landscape
make_lil_lands <- function(n_x, n_y, n_sims, n_pseudo, spat_config) {

    n_plants <- n_x * n_y

    land <- array(0L, c(n_x, n_y, n_sims))

    for (i in 1:n_sims) {
        land[1,1,i] <- 1L
        if (n_pseudo > 0) {
            if (spat_config == "random") {
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
            } else if (spat_config == "no virus") {
                stopifnot(n_pseudo < n_plants)
                k <- sample.int(n_plants - 1L, n_pseudo)
                x <- k - n_x * (k %/% n_x) + 1L
                y <- k %/% n_x + 1L
            } else if (spat_config == "diagonal") {
                stopifnot(n_pseudo == 3)
                x <- 3:1
                y <- 1:3
            } else if (spat_config == "near virus") {
                stopifnot(n_pseudo == 3)
                x <- c(1, 2, 2)
                y <- c(2, 1, 2)
            } else if (spat_config == "far virus") {
                stopifnot(n_pseudo == 3)
                x <- c(2, 3, 3)
                y <- c(3, 2, 3)
            } else if (spat_config == "over virus") {
                stopifnot(n_pseudo == 3)
                land[1,1,i] <- 3L
                x <- c(1, 2)
                y <- c(2, 1)
            } else stop("improper spat_config")
            land[cbind(x,y,i)] <- 2L
        }
    }

    return(land)

}








#' Simulate a little (3 x 3) plantscape.
#'
#' @details # Spatial configuration
#'
#' In the below descriptions, `v` = virus-infected plant,
#' `p` = *Pseudomonas*-inhabited plant, `B` = plant with both
#'
#' `"random"`: random configuration (including on initial virus plant)
#'
#' `"no virus`: random configuration (never starting on initial virus plant)
#'
#' `"diagonal"`:
#' ```
#' |---|---|---|
#' | v | - | p |
#' | - | p | - |
#' | p | - | - |
#' |---|---|---|
#' ```
#'
#' `"near virus"`:
#' ```
#' |---|---|---|
#' | v | p | - |
#' | p | p | - |
#' | - | - | - |
#' |---|---|---|
#' ```
#'
#' `"far virus"`:
#' ```
#' |---|---|---|
#' | v | - | - |
#' | - | - | p |
#' | - | p | p |
#' |---|---|---|
#' ```
#'
#' `"over virus"`:
#' ```
#' |---|---|---|
#' | B | p | - |
#' | p | - | - |
#' | - | - | - |
#' |---|---|---|
#' ```
#'
#'
#' @param n_sims Single integer giving the number of simulations to run.
#'     Must be `>= 1`.
#' @param n_pseudo Single integer giving the number of patches containing
#'     *Pseudomonas*. Must be `>= 0`.
#' @param N0 If a single numeric, it gives the mean for the lognormal
#'     distribution used to generate non-winged aphid densities per patch.
#'     If a numeric matrix, it must be 3x3 and gives the abundance of non-winged
#'     aphids per plant that is the same across repetitions.
#'     If a 3D array, it must be 3x3x`n_sims` and gives the abundance of
#'     non-winged aphids per plant across repetitions.
#' @param sd_N Single numeric giving the standard deviation for the lognormal
#'     distribution used to generate non-winged aphid densities per patch.
#'     Must be `>= 0`.
#' @param Y0 If a single numeric, it gives the starting wasp density that will
#'     be the same across repetitions.
#'     If a length-`n_sims` vector, it gives the starting wasp density that can
#'     differ across repetitions.
#' @param W0 If a single numeric, it gives the mean for the lognormal
#'     distribution used to generate winged aphid densities per patch.
#'     If a numeric matrix, it must be 3x3 and gives the abundance of winged
#'     aphids per plant that is the same across repetitions.
#'     If a 3D array, it must be 3x3x`n_sims` and gives the abundance of
#'     winged aphids per plant across repetitions.
#'     Defaults to `0`.
#' @param sd_W Single numeric giving the standard deviation for the lognormal
#'     distribution used to generate winged aphid densities per patch.
#'     Must be `>= 0`.
#'     Defaults to `0`.
#' @param spat_config A single string indicating the spatial configuration of
#'     *Pseudomonas* inhabited patches.
#'     See "Spatial configuration" section below for details.
#'     Defaults to `"random"`.
#' @param force_N_distr Single logical for whether to force starting aphid
#'     densities to exactly match a given mean and standard deviation.
#'     Otherwise, they are drawn from distributions with those parameters, but
#'     the exact values of the summary stats will differ stochastically.
#'     Defaults to `TRUE`.
#' @param ... Other parameters for functions [sim_plantscape()],
#'     [make_disease_ptr()], or [make_insects_ptr()].
#' @inheritParams sim_plantscape
#' @inheritParams make_disease_ptr
#' @inheritParams make_insects_ptr
#'
#' @returns A tibble with columns following the description in the
#'     "Summarizing" section in the docs for function [sim_plantscape()].
#'
#' @export
#'
#'
#'
lil_plantscape <- function(n_sims,
                           n_pseudo,
                           pseudo_repel,
                           pseudo_surv,
                           zeta,
                           N0,
                           sd_N,
                           Y0,
                           W0 = 0,
                           sd_W = 0,
                           radius = 1,
                           virus_attract = 1,
                           fly_p = 0.05,
                           p_load_alate = 0.5,
                           p_load_plant = 0.5,
                           spat_config = "random",
                           force_N_distr = FALSE,
                           ...) {


    single_integer(n_pseudo, "n_pseudo", .min = 0)
    is_type(spat_config, "spat_config", "character", L = 1L)

    spat_config <- match.arg(spat_config, c("random", "no virus", "diagonal",
                                            "near virus", "far virus", "over virus"))

    n_x <- 3L
    n_y <- 3L

    landscapes <- make_lil_lands(n_x, n_y, n_sims, n_pseudo, spat_config)


    out <- plantscape_shared(n_x = n_x,
                             n_y = n_y,
                             n_sims = n_sims,
                             N0 = N0,
                             sd_N = sd_N,
                             W0 = W0,
                             sd_W = sd_W,
                             Y0 = Y0,
                             force_N_distr = force_N_distr,
                             fly_p = fly_p,
                             radius = radius,
                             virus_attract = virus_attract,
                             pseudo_repel = pseudo_repel,
                             pseudo_surv = pseudo_surv,
                             zeta = zeta,
                             p_load_alate = p_load_alate,
                             p_load_plant = p_load_plant,
                             landscapes = landscapes,
                             ...)

    return(out)

}

