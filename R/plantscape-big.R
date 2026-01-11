




#' Simulate large (100 x 100) plantscape.
#'
#' @param landscape Integer matrix of size 100 x 100 with the types of each
#'     plant in each landscape. It's assumed that rows are x the dimension and
#'     columns are the y dimension.
#'     Values in each cell give the state of the plant:
#'     `0` indicates nothing on plant,
#'     `1` indicates just virus on plant (infectious),
#'     `2` indicates just *Pseudomonas* on plant,
#'     `3` indicates both virus and *Pseudomonas* on plant.
#'     Values < 0 or > 3 are not allowed.
#' @param N0 If a single numeric, it gives the mean for the lognormal
#'     distribution used to generate non-winged aphid densities per patch.
#'     If a numeric matrix, it must be 100 x 100 and gives the abundance of non-winged
#'     aphids per plant that is the same across repetitions.
#'     If a 3D array, it must be 100 x 100 x `n_sims` and gives the abundance of
#'     non-winged aphids per plant across repetitions.
#' @param W0 If a single numeric, it gives the mean for the lognormal
#'     distribution used to generate winged aphid densities per patch.
#'     If a numeric matrix, it must be 100 x 100 and gives the abundance of winged
#'     aphids per plant that is the same across repetitions.
#'     If a 3D array, it must be 100 x 100 x `n_sims` and gives the abundance of
#'     winged aphids per plant across repetitions.
#'     Defaults to `0`.
#' @param radius Max distance that alates will travel between plants.
#'     Must be >= 1.
#'     See "Radius" section in the description for [lil_plantscape()] for details.
#'     Defaults to `pop_info$radius`.
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
big_plantscape <- function(n_sims,
                           landscape,
                           pseudo_repel,
                           pseudo_surv,
                           zeta,
                           N0,
                           sd_N,
                           Y0,
                           W0 = 0,
                           sd_W = 0,
                           radius = pop_info$radius,
                           virus_attract = 1,
                           fly_p = 0.05,
                           p_load_alate = 0.5,
                           p_load_plant = 0.5,
                           ...) {

    n_x <- 100L
    n_y <- 100L

    is_landscape_array(landscape, "landscape", c(n_x, n_y))

    landscape <- array(landscape, dim = c(dim(landscape), 1L))

    out <- plantscape_shared(n_x = n_x,
                             n_y = n_y,
                             n_sims = n_sims,
                             N0 = N0,
                             sd_N = sd_N,
                             W0 = W0,
                             sd_W = sd_W,
                             Y0 = Y0,
                             force_N_distr = FALSE,
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

