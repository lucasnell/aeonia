
#'
#' Sensitivity via Latin hypercube sampling (LHS) with partial rank
#' correlation coefficients (PRCC)
#'

# make sure these come before preamble to avoid masking dplyr::select:
library(lhs)
library(ppcor)

source("_scripts/01-sensitivity/00-preamble.R")


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
                  wasp_disp_m0 = 0.3 * c(0, 2),
                  wasp_disp_m1 = 0.349 * c(0, 2))

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



# if (!file.exists("_scripts/lhs-sims.rds")) {
#     # Takes ~6 min with 108 combos of 1000 reps
#     t0 <- Sys.time()
#     set.seed(891393509)
#     lhs_sims <- lhs_df |>
#         pmap(\(mu_Y, mu_N, n_pseudo, other_pars) {
#             other_pars |>
#                 pmap(one_combo, mu_Y = mu_Y, mu_N = mu_N, n_pseudo = n_pseudo) |>
#                 list_rbind() |>
#                 select(all_of(par_names), everything())
#         }, .progress = .prog_args)
#     write_rds(lhs_sims, "_scripts/lhs-sims.rds", compress = "gz")
#     t1 <- Sys.time()
#     t1 - t0
# } else {
#     lhs_sims <- read_rds("_scripts/lhs-sims.rds")
# }

lhs_sims <- read_rds("~/_globus/lhs-sims.rds") |>
    map(\(x) mutate(x, n_pseudo = as.integer(n_pseudo)))


# Average across each set of 100 simulations
lhs_sims_summs <- lhs_sims |>
    map(\(x) {
        x[1,names(struct_pars)] |>
            mutate(results = x |>
                       group_by(across(all_of(names(vary_pars)))) |>
                       summarize(across(p_alates:outbreak_size, mean), .groups = "drop") |>
                       list())
    }) |>
    list_rbind() |>
    mutate(pcors = map(results,
                       \(.res) {
                           map(c("log_alates", "log_aphids", "infect_time"),
                               \(x) {
                                   vars_x <- c(names(vary_pars), x)
                                   pcor_x = .res |> select(all_of(vars_x)) |>
                                       pcor(method = "kendall")
                                   pcor_x$estimate |>
                                       as.data.frame() |>
                                       as_tibble() |>
                                       select(all_of(x)) |>
                                       slice(-n())
                               }) |>
                               list_cbind() |>
                               mutate(param = names(vary_pars)) |>
                               select(param, everything())
                       }))




lhs_sims_summs |>
    filter(mu_Y == max(mu_Y), n_pseudo == min(n_pseudo), mu_N == min(mu_N), K == median(K)) |>
    getElement("pcors") |>
    getElement(1) |>
    pivot_longer(-param, names_to = "output", values_to = "pcor") |>
    ggplot(aes(pcor, param)) +
    geom_vline(xintercept = 0, color = "gray70") +
    geom_segment(aes(xend = 0, yend = param)) +
    geom_point() +
    facet_wrap(~output, scales = "free_x")


pretty_params <- function(x) {
    case_when(x == "alpha" ~ "&alpha;",
              x == "beta" ~ "&beta;",
              x == "wasp_disp_m0" ~ "&delta;<sub>0</sub>",
              x == "wasp_disp_m1" ~ "&delta;<sub>z</sub>",
              .default = x)
}

lhs_sims_summs |>
    filter(mu_Y == -Inf, mu_N == median(mu_N)) |>
    select(-results, -mu_Y, -mu_N) |>
    unnest(pcors) |>
    rename(`n<sub>P</sub>` = n_pseudo) |>
    mutate(param = pretty_params(param) |>
               factor(levels = pretty_params(names(vary_pars))),
           K = factor(K)) |>
    ggplot(aes(infect_time, param, fill = K)) +
    geom_vline(xintercept = 0, color = "gray60", linewidth = 1) +
    # geom_segment(aes(xend = 0, yend = param), linewidth = 1,
    #              position = position_dodge(width = 0.2)) +
    # geom_point(size = 3, position = position_dodge(width = 0.2)) +
    geom_col(position = position_dodge(width = 0.5)) +
    facet_wrap(~ `n<sub>P</sub>`, labeller = \(x) label_both(x, sep = " = ")) +
    labs(x = "Partial Kendall's rank correlation (&tau;)") +
    scale_fill_viridis_d(begin = 0.1, end = 0.9) +
    theme(strip.text.x = element_markdown(angle = 0),
          axis.title.y = element_blank(),
          axis.text.y = element_markdown(),
          axis.title.x = element_markdown(),
          panel.grid.major.y = element_line(linewidth = 0.5, color = "gray80"))



lhs_sims_summs |>
    filter(mu_Y == -Inf, mu_N == median(mu_N)) |>
    select(-pcors, -mu_Y, -mu_N) |>
    unnest(results) |>
    select(n_pseudo, K, all_of(names(vary_pars)), infect_time) |>
    pivot_longer(all_of(names(vary_pars)), names_to = "param",
                 values_to = "x") |>
    mutate(n_pseudo = paste("n<sub>P</sub> = ", n_pseudo),
           param = pretty_params(param) |>
               factor(levels = pretty_params(names(vary_pars))),
           K = factor(K, levels = sort(unique(K)),
                      labels = sort(unique(K)) |> prettyNum(big.mark = ","))) |>
    ggplot(aes(x, infect_time, color = K)) +
    geom_point(shape = 1) +
    facet_grid(n_pseudo ~ param, scales = "free", switch = "x", axes = "all") +
    labs(y = "Days to 5 plants infected") +
    guides(color = guide_legend(override.aes = list(shape = 19, size = 3))) +
    scale_color_viridis_d(begin = 0.1, end = 0.9) +
    theme(strip.text.x = element_markdown(),
          strip.text.y = element_markdown(angle = 0, hjust = 0),
          axis.title.y = element_markdown(),
          axis.text.y = element_markdown(),
          axis.title.x = element_blank(),
          strip.placement = 'outside')





