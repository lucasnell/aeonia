

# make an array of starting aphid abundances (also works for winged aphids)
make_aphids0 <- function(N0, sd_N, n_x, n_y, n_sims, force_N_distr) {

    n_plants <- n_x * n_y

    if (length(N0) == 1) {

        N0 <- array(N0, c(n_x, n_y, n_sims))

        if (sd_N > 0) {

            # Convert from mean and sd of lognormal distribution to parameters
            # to use for lognormal (mean and sd of underlying normal distribution):
            mu_N <- log(N0[1,1,1]^2 / sqrt(N0[1,1,1]^2 + sd_N^2))
            sigma_N <- sqrt(log(1 + sd_N^2 / N0[1,1,1]^2))

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


plantscape_shared <- function(n_x,
                              n_y,
                              n_sims,
                              N0,
                              sd_N,
                              W0,
                              sd_W,
                              M0,
                              sd_M,
                              Y0,
                              force_N_distr,
                              radius,
                              virus_attract,
                              pseudo_repel,
                              pseudo_surv,
                              zeta,
                              p_load_alate,
                              p_load_plant,
                              landscapes,
                              ...) {

    single_integer(n_x, "n_x", .min = 1, .max = 100e3)
    single_integer(n_y, "n_y", .min = 1, .max = 100e3)
    single_integer(n_sims, "n_sims", .min = 1)
    is_type(N0, "N0", is.numeric, .min = 0)
    single_number(sd_N, "sd_N", .min = 0)
    is_type(W0, "W0", is.numeric, .min = 0)
    single_number(sd_W, "sd_W", .min = 0)
    is_type(M0, "M0", is.numeric, .min = 0)
    single_number(sd_M, "sd_M", .min = 0)
    is_type(Y0, "Y0", is.numeric, L = c(1L, n_sims), .min = 0)
    is_type(force_N_distr, "force_N_distr", "logical", L = 1L)
    single_number(radius, "radius", .min = 1)
    single_number(virus_attract, "virus_attract", .min = 1)
    single_number(pseudo_repel, "pseudo_repel", .min = 1)
    single_number(pseudo_surv, "pseudo_surv", .min = 0, .max = 1)
    single_number(zeta, "zeta", .min = 0, .max = 1)
    single_number(p_load_alate, "p_load_alate", .min = 0, .max = 1)
    single_number(p_load_plant, "p_load_plant", .min = 0, .max = 1)
    is_landscape_array(landscapes, "landscapes", c(n_x, n_y, n_sims))

    N0 <- make_aphids0(N0, sd_N, n_x, n_y, n_sims, force_N_distr)
    W0 <- make_aphids0(W0, sd_W, n_x, n_y, n_sims, force_N_distr)
    M0 <- make_aphids0(M0, sd_M, n_x, n_y, n_sims, force_N_distr)
    if (length(Y0) == 1) Y0 <- rep(Y0, n_sims)



    .args <- list(insects = list(pseudo_surv = pseudo_surv,
                                 zeta = zeta),
                  disease = list(radius = radius,
                                 virus_attract = virus_attract,
                                 pseudo_repel = pseudo_repel,
                                 p_load_alate = p_load_alate,
                                 p_load_plant = p_load_plant),
                  plantscape = list(landscapes = landscapes,
                                    N0 = N0,
                                    W0 = W0,
                                    M0 = M0,
                                    Y0 = Y0,
                                    insect_ptr = NULL,
                                    disease_ptr = NULL))
    other_args <- list(...)
    if (length(other_args) > 0) {

        if (is.null(names(other_args)) || any(names(other_args) == "")) {
            stop("...  must be entirely named items")
        }

        .formals <- c(insects = make_insects_ptr, disease = make_disease_ptr,
                      plantscape = sim_plantscape) |>
            lapply(\(x) names(formals(x)))
        all_formals <- do.call(c, .formals) |> unname()

        if (!all(names(other_args) %in% all_formals)) {
            print(names(other_args)[!names(other_args) %in% all_formals])
            stop("\nThe names printed above (inside `...`) do not match args ",
                 "in sim_plantscape, make_insects_ptr, or make_disease_ptr")
        }

        for (x in names(.args)) {
            nm_args <- names(other_args)[names(other_args) %in% .formals[[x]]]
            for (n in nm_args) .args[[x]][[n]] <- other_args[[n]]
        }

    }

    .args$plantscape[["insect_ptr"]] <- do.call(make_insects_ptr, .args$insects)
    .args$plantscape[["disease_ptr"]] <- do.call(make_disease_ptr, .args$disease)

    out <- do.call(sim_plantscape, .args$plantscape)

    return(out)
}

