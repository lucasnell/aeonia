#!/usr/bin/env -S Rscript --vanilla


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Must have exactly 2 arguments")

n_threads <- suppressWarnings(as.integer(args[[1]]))
if (is.na(n_threads)) stop("Argument #1 must be an integer")

sim_idx <- suppressWarnings(as.integer(args[[2]]))
if (is.na(sim_idx) || sim_idx < 1 || sim_idx > 6) {
    stop("Argument #2 must be an integer between 1 and 6")
}


suppressPackageStartupMessages({
    library(sensobol)
    library(tidyverse)
    library(aeonia)
})


# Set threads for simulations:
options("mc.cores" = n_threads)




#'
#' Parameters that aren't especially interesting or relevant:
#'   - fly_p (set to 0.05 below)
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
                  epsilon = c(0, 4),
                  wasp_disp_m0 = 0.3 * c(0, 2),
                  wasp_disp_m1 = 0.349 * c(0, 2))

par_names <- c(names(struct_pars), names(vary_pars))





one_combo <- function(mu_Y, mu_N, n_pseudo, K,
                      B, alpha, beta, epsilon, wasp_disp_m0, wasp_disp_m1,
                      n_sims, ...) {

    fly_p <- 0.05
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
               B = B, K = K, alpha = alpha, beta = beta, epsilon = epsilon,
               wasp_disp_m0 = wasp_disp_m0, wasp_disp_m1 = wasp_disp_m1)

    return(out)

}



# sensobol arguments
N <- 2^12
R <- 2e3
type <- "norm"
conf <- 0.95



# Takes just a second or two
set.seed(141672556)
sobol_mats <- map(1:nrow(struct_pars), \(i) {
    pars_i <- names(vary_pars)
    no_pseudo <- struct_pars$n_pseudo[[i]] == 0
    no_wasps <- struct_pars$mu_Y[[i]] == -Inf
    if (no_pseudo) pars_i <- pars_i[!pars_i %in% c("beta", "B")]
    if (no_wasps) pars_i <- pars_i[!grepl("^wasp_disp", pars_i)]
    mat <- sobol_matrices(N = N, params = pars_i)
    for (n in pars_i) {
        p1 <- min(vary_pars[[n]])
        p2 <- max(vary_pars[[n]])
        mat[,n] <- p1 + mat[,n] * (p2 - p1)
    }
    if (no_pseudo) {
        mat <- cbind(mat, beta = 0)
        mat <- cbind(mat, B = 0)
    }
    if (no_wasps) {
        mat <- cbind(mat, wasp_disp_m0 = 0)
        mat <- cbind(mat, wasp_disp_m1 = 0)
    }
    return(mat)
})


# Should be TRUE
sobol_mats |>
    map_lgl(\(x) identical(sort(names(vary_pars)), sort(colnames(x)))) |>
    all()





one_struct_sobol_sims <- function(i) {

    mat <- sobol_mats[[i]]
    mu_Y <- struct_pars[["mu_Y"]][[i]]
    mu_N <- struct_pars[["mu_N"]][[i]]
    n_pseudo <- struct_pars[["n_pseudo"]][[i]]
    K <- struct_pars[["K"]][[i]]

    sim_outs <- map(1:nrow(mat), \(j) {
        B <- mat[j,"B"]
        alpha <- mat[j,"alpha"]
        beta <- mat[j,"beta"]
        epsilon <- mat[j,"epsilon"]
        wasp_disp_m0 <- mat[j,"wasp_disp_m0"]
        wasp_disp_m1 <- mat[j,"wasp_disp_m1"]
        sim <- one_combo(mu_Y = mu_Y, mu_N = mu_N, n_pseudo = n_pseudo, K = K,
                         B, alpha, beta, epsilon, wasp_disp_m0, wasp_disp_m1,
                         n_sims = 100)
        out <- sim[1,names(vary_pars)]
        yvars <- c("p_alates", "log_aphids", "aphids", "log_alates", "alates",
                   "log_parasitized", "parasitized", "log_mummies", "mummies",
                   "log_wasps", "wasps", "infect_time", "outbreak_size")
        for (y in yvars) {
            out[[y]] <- mean(sim[[y]])
        }
        return(out)
    }) |>
        list_rbind()

    return(sim_outs)
}



dn <- length(sobol_mats) %/% 6L
start <- (sim_idx - 1L) * dn + 1L
end <- sim_idx * dn


# Takes ~2.6 hrs with 108 / 6 combos with N=2^12 and 20 threads:
t0 <- Sys.time()
sobol_sims <- map(start:end, one_struct_sobol_sims)
write_rds(sobol_sims, sprintf("sobol-sims-%i.rds", sim_idx), compress = "gz")
t1 <- Sys.time()
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(sprintf("Took %.2f hours\n", difftime(t1, t0, units = "hours")))
cat(paste(rep("-", 80), collapse = ""), "\n")
cat(paste(rep("-", 80), collapse = ""), "\n")




