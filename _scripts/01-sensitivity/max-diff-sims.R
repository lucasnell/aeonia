
#'
#' Find maximum effect (increase if positive) of Pseudomonas on outbreak size
#'
#'
#' This was then run on BioHPC in an interactive job started with the following:
#'
#' cd /home2/lan68/max-sims
#' srun -n 1 -c 50 -N 1 --mem=80G --time=1-20:00:00 \
#'   --job-name="max-sims" --pty bash -l
#'
#'
#' Because of issues pasting directly into the R session, I pasted the
#' function `one_combo` below first into the `one-combo.R` file.
#'
#' Then I started R using:
#' R --vanilla
#'
#' Then I ran the following R code interactively:
#'

library(tidyverse)
library(aeonia)
library(nloptr)
library(estimatePMR)


slurm_nt <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (slurm_nt == "") stop("Variable SLURM_CPUS_PER_TASK not found!")
n_threads <- suppressWarnings(as.integer(slurm_nt))
if (is.na(n_threads)) stop("Argument must be an integer")
# Set threads for simulations:
options("mc.cores" = n_threads)

source("one-combo.R")


#
# one_combo <- function(n_pseudo, alate_dens,
#                       Y0, mean_N, sd_N, K, virus_attract, pseudo_repel, pseudo_surv,
#                       zeta, spat_config,
#                       n_sims, ...) {
#
#     stopifnot(alate_dens %in% 0:1)
#
#     fly_p <- 0.05
#     epsilon <- 1
#     p_load_alate <- 0.5
#     p_load_plant <- 0.5
#
#     n_x <- 3L
#     n_y <- 3L
#     n_plants <- n_x * n_y
#
#     # Convert from mean and sd of underlying normal distribution to parameters
#     # to use for lognormal:
#     mu_N <- log(mean_N^2 / sqrt(mean_N^2 + sd_N^2))
#     sigma_N <- sqrt(log(1 + sd_N^2 / mean_N^2))
#
#     if (alate_dens == 0L) {
#         alate_slope <- 0
#         alate_max  <- 0.075
#     } else {
#         alate_slope <- pop_info$alate_slope
#         alate_max <- 1
#     }
#
#
#     land <- array(0L, c(n_x, n_y, n_sims))
#     N0 <- array(0.0, c(n_x, n_y, n_sims))
#     for (i in 1:n_sims) {
#         land[1,1,i] <- 1L
#         if (n_pseudo > 0) {
#             if (spat_config == 0L) {
#                 k <- sample.int(n_plants - 1L, n_pseudo)
#                 x <- k - n_x * (k %/% n_x) + 1L
#                 y <- k %/% n_x + 1L
#             } else if (spat_config == 1L) {
#                 x <- 3:1
#                 y <- 1:3
#             } else if (spat_config == 2L) {
#                 x <- c(1, 2, 2)
#                 y <- c(2, 1, 2)
#             } else if (spat_config == 3L) {
#                 x <- c(2, 3, 3)
#                 y <- c(3, 2, 3)
#             } else if (spat_config == 4L) {
#                 land[1,1,i] <- 3L
#                 x <- c(1, 2)
#                 y <- c(2, 1)
#             }
#             land[cbind(x,y,i)] <- 2L
#         }
#         N0[,,i] <- rlnorm(n_plants, mu_N, sigma_N)
#     }
#
#     insect_args <- list(K = K, pseudo_surv = pseudo_surv, fly_p = fly_p,
#                         zeta = zeta, alate_slope = alate_slope,
#                         alate_max = alate_max)
#     plant_args <- list(landscapes = land,
#                        max_t = 100,
#                        N0 = N0,
#                        W0 = array(0.0, c(n_x, n_y, n_sims)),
#                        Y0 = Y0,
#                        virus_attract = virus_attract,
#                        pseudo_repel = pseudo_repel,
#                        epsilon = epsilon,
#                        infect_time_n = 5,
#                        p_load_alate = p_load_alate,
#                        p_load_plant = p_load_plant,
#                        radius = 1,
#                        infect_stop = FALSE,
#                        summ = "all")
#     other_args <- list(...)
#     if (length(other_args) > 0) {
#         stopifnot(!is.null(names(other_args)) && !any(names(other_args) == ""))
#         stopifnot(all(names(other_args) %in% c(names(formals(sim_plantscape)),
#                                                names(formals(make_insect_ptr)))))
#         not_allowed <- c("landscapes", "N0", "W0", "Y0")
#         if (any(names(other_args) %in% not_allowed)) {
#             not_allowed <- names(other_args)[names(other_args) %in% not_allowed]
#             stop("The following are not allowed in `make_arg_list`: ",
#                  paste(not_allowed, collapse = ", "))
#         }
#
#         nm_insect_args <- names(other_args)[names(other_args) %in%
#                                                 names(formals(make_insect_ptr))]
#         for (n in nm_insect_args) insect_args[[n]] <- other_args[[n]]
#
#         nm_plant_args <- names(other_args)[names(other_args) %in%
#                                                names(formals(sim_plantscape))]
#         for (n in nm_plant_args) {
#             if (n %in% c("N0", "W0", "Y0", "landscapes")) {
#                 plant_args[[n]] <- array(other_args[[n]], c(n_x, n_y, n_sims))
#             } else plant_args[[n]] <- other_args[[n]]
#         }
#     }
#
#     plant_args[["insect_ptr"]] <- do.call(make_insect_ptr, insect_args)
#
#     out <- do.call(sim_plantscape, plant_args)
#
#     # Now add arguments to the output dataframe (not including any in `...`):
#     in_args <- as.list(match.call(expand.dots = FALSE))[-1]  # -1 removes fxn name
#     # remove dots if present:
#     if ("..." %in% names(in_args)) { # << this if statement is important!
#         in_args <- in_args[-which(names(in_args) == "...")]
#     }
#     for (n in names(in_args)) out[[n]] <- in_args[[n]]
#
#     # Make these args come first:
#     out <- out |> select(all_of(names(in_args)), everything())
#
#     return(out)
#
# }




# Force parameters inside bounds:
trans_pars <- function(p) {
    lower_bounds <- c(1, 10, 0, 1250, 1, 1, 0.85, 0)
    upper_bounds <- c(9, 100, 50, 25000, 5, 5, 1, 1)
    x <- exp(-exp(-p)) * (upper_bounds - lower_bounds) + lower_bounds
    return(x)
}
# Convert to version used in objective function:
inv_trans_pars <- function(x) {
    lower_bounds <- c(1, 10, 0, 1250, 1, 1, 0.85, 0)
    upper_bounds <- c(9, 100, 50, 25000, 5, 5, 1, 1)
    p <- -log(-log((x - lower_bounds) / (upper_bounds - lower_bounds)))
    return(p)
}



dob_pred_f <- function(pars) {
    pars <- trans_pars(pars)
    names(pars) <- c("Y0",
                     "mean_N",
                     "sd_N",
                     "K" ,
                     "virus_attract",
                     "pseudo_repel",
                     "pseudo_surv",
                     "zeta")
    args <- c(as.list(pars), list(alate_dens = 1L,
                                  spat_config = 1L,
                                  n_sims = 1e3L,
                                  summ = "all"))
    sims3 <- do.call(one_combo, c(args, list(n_pseudo = 3L))) |>
        getElement("outbreak_size")
    sims0 <- do.call(one_combo, c(args, list(n_pseudo = 0L))) |>
        getElement("outbreak_size")
    return(-1 * (mean(sims3) - mean(sims0)))
}


x0 <- c(Y0 = 2.5,
        mean_N = 25,
        sd_N = 1e-6,
        K = 23e3,
        virus_attract = 1+1e-6,
        pseudo_repel = 1+1e-6,
        pseudo_surv = 0.85+1e-6,
        zeta = 1e-6)


# Takes ~1.4 min with 100 threads
t0 <- Sys.time()
op <- nloptr::neldermead(inv_trans_pars(x0), dob_pred_f)
t1 <- Sys.time()
t1 - t0

# Takes ~33 sec with 100 threads
t0 <- Sys.time()
op2 <- nloptr::bobyqa(inv_trans_pars(x0), dob_pred_f)
t1 <- Sys.time()
t1 - t0

# Takes ~2.2 min with 100 threads
t0 <- Sys.time()
op3 <- optim(inv_trans_pars(x0), dob_pred_f)
t1 <- Sys.time()
t1 - t0

# Takes ~ with 100 threads
t0 <- Sys.time()
wop <- winnowing_optim(dob_pred_f,
                       rep(-10, 8) |> set_names(names(x0)),
                       rep(100, 8) |> set_names(names(x0)),
                       n_bevals = 100L, n_boxes = 100L,
                       fn_args = list(), n_outputs = c(75L, 50L, 25L),
                       optimizers = c(nloptr::bobyqa,
                                      nloptr::bobyqa,
                                      nloptr::neldermead),
                       verbose = TRUE)
t1 <- Sys.time()
t1 - t0


x1 <- c(4.43644021206967, 10.0009697729649, 2.74678259005212e-124,
        24704.5399820386, 1, 1.03461566412611, 0.850005035666408,
        9.98782556332457e-21) |>
    set_names(names(x0))

dob_pred_f(inv_trans_pars(x0))
dob_pred_f(inv_trans_pars(x1))
