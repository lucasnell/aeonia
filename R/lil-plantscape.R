


# make a little landscape
make_lil_land <- function(n_x, n_y, n_sims, n_pseudo, spat_config) {

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




# make an array of starting aphid abundances (also works for winged aphids)
make_aphids0 <- function(N0, sd_N, n_x, n_y, n_sims, force_N_distr) {

    n_plants <- n_x * n_y

    if (length(N0) == 1) {

        N0 <- array(N0, c(n_x, n_y, n_sims))

        if (sd_N > 0) {

            # Convert from mean and sd of lognormal distribution to parameters
            # to use for lognormal (mean and sd of underlying normal distribution):
            mu_N <- log(N0^2 / sqrt(N0^2 + sd_N^2))
            sigma_N <- sqrt(log(1 + sd_N^2 / N0^2))

            for (i in 1:n_sims) {
                if (force_N_distr) {
                    N0[,,i] <- rnorm(n_plants, 0, 1) |>
                        (\(x) ((x - mean(x)) / sd(x)) * sigma_N + mu_N)() |>
                        exp()
                } else {
                    N0[,,i] <- rlnorm(n_plants, mu_N, sigma_N)
                }
            }

        }
    } else {
        if (!is.array(N0) || ! length(dim(N0)) %in% 2:3)
            stop("If length(N0) > 1, it must be a 2D or 3D array.")
        if (!identical(dim(N0)[1:2], c(n_x, n_y)))
            stop("If length(N0) > 1, it must be an array with 3 rows & 3 cols.")
        if (sd_N > 0) stop("If length(N0) > 1, then sd_N must be 0.")
        if (length(dim(N0)) == 2L) N0 <- array(N0, c(n_x, n_y, n_sims))
        if (dim(N0)[3] != n_sims) stop("dim(N0)[3] must equal n_sims")
    }

    return(N0)


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
#'     [make_disease_ptr()], or [make_insect_ptr()].
#' @inheritParams sim_plantscape
#' @inheritParams make_disease_ptr
#' @inheritParams make_insect_ptr
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
                           force_N_distr = TRUE,
                           ...) {


    stopifnot(is.numeric(n_sims) && length(n_sims) == 1L && n_sims > 0)
    stopifnot(is.numeric(n_pseudo) && length(n_pseudo) == 1L && n_pseudo >= 0)
    stopifnot(is.numeric(N0) && all(N0 >= 0))
    stopifnot(is.numeric(Y0) && all(Y0 >= 0) && length(Y0) %in% c(1, n_sims))
    stopifnot(is.numeric(sd_N) && length(sd_N) == 1L && sd_N >= 0)
    stopifnot(is.numeric(fly_p) && length(fly_p) == 1L)
    stopifnot(is.numeric(epsilon) && length(epsilon) == 1L)
    stopifnot(is.numeric(p_load_alate) && length(p_load_alate) == 1L)
    stopifnot(is.numeric(p_load_plant) && length(p_load_plant) == 1L)
    stopifnot(is.character(spat_config) && length(spat_config) == 1L)
    stopifnot(inherits(force_N_distr, "logical") && length(force_N_distr) == 1L)

    spat_config <- match.arg(spat_config, c("random", "no virus", "diagonal",
                                            "near virus", "far virus", "over virus"))

    n_x <- 3L
    n_y <- 3L
    n_plants <- n_x * n_y

    land <- make_lil_land(n_x, n_y, n_sims, n_pseudo, spat_config)

    N0 <- make_aphids0(N0, sd_N, n_x, n_y, n_sims, force_N_distr)
    W0 <- make_aphids0(W0, sd_W, n_x, n_y, n_sims, force_N_distr)
    if (length(Y0) == 1) Y0 <- rep(Y0, n_sims)

    .args <- list(insect = list(fly_p = fly_p),
                 disease = list(radius = radius,
                                virus_attract = virus_attract,
                                pseudo_repel = pseudo_repel,
                                p_load_alate = p_load_alate,
                                p_load_plant = p_load_plant),
                 plantscape = list(landscapes = land,
                                   N0 = N0,
                                   W0 = W0,
                                   Y0 = Y0,
                                   insect_ptr = NULL,
                                   disease_ptr = NULL))
    other_args <- list(...)
    if (length(other_args) > 0) {

        stopifnot(!is.null(names(other_args)) && all(names(other_args) != ""))

        .formals <- c(insect = make_insect_ptr, disease = make_disease_ptr,
                      plantscape = sim_plantscape) |>
            lapply(\(x) names(formals(x)))
        all_formals <- do.call(c, .formals) |> unname()

        if (!all(names(other_args) %in% all_formals)) {
            print(names(other_args)[!names(other_args) %in% all_formals])
            stop("\nThe names printed above (inside `...`) do not match args ",
                 "in sim_plantscape, make_insect_ptr, or make_disease_ptr")
        }

        for (x in names(.args)) {
            nm_args <- names(other_args)[names(other_args) %in% .formals[[x]]]
            for (n in nm_args) .args[[x]][[n]] <- other_args[[n]]
        }

    }

    .args$plantscape[["insect_ptr"]] <- do.call(make_insect_ptr, .args$insect)
    .args$plantscape[["disease_ptr"]] <- do.call(make_disease_ptr, .args$disease)

    out <- do.call(sim_plantscape, .args$plantscape)

    return(out)

}

