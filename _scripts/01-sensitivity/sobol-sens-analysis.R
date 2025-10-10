
#'
#' Analysis of sobol sensitivity output from cluster
#'

library(sensobol)

source("_scripts/01-sensitivity/00-preamble.R")


#' Directory containing output from `sobol-sens.R`:
sobol_dir <- "~/_globus"


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


#' Should be TRUE
sobol_mats |>
    map_lgl(\(x) identical(sort(names(vary_pars)), sort(colnames(x)))) |>
    all()


#' Output from `sobol-sens.R`:
sobol_sims <- list.files(sobol_dir, "^sobol.*.rds", full.names = TRUE) |>
    map(read_rds) |>
    do.call(what = c)


j <- 35

sim_outs <- sobol_sims[[j]]
mat <- sobol_mats[[j]]


scatter <- function(sim_outs, mat) {

    y <- sim_outs[["infect_time"]]
    params <- names(vary_pars)
    for (p in params) {
        if (sd(sim_outs[[p]]) == 0) params <- params[params != p]
    }

    sim_outs |>
        rename(y = infect_time) |>
        select(all_of(params), y) |>
        pivot_longer(all_of(params), names_to = "param", values_to = "x") |>
        mutate(param = param |>
                   pretty_params() |>
                   factor(levels = pretty_params(params))) |>
        ggplot(aes(x, y)) +
        geom_point(size = 0.7, alpha = 0.1) +
        stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE) +
        facet_wrap(~ param, scales = "free_x") +
        labs(x = "Value", y = "Days to 5 plants infected") +
        theme(strip.text = element_markdown())
}


# rm(dt, out, da, i, cols)

multi_scatter <- function(mat, sim_outs, N) {

    Y <- sim_outs[["infect_time"]]
    params <- names(vary_pars)
    for (p in params) {
        if (sd(sim_outs[[p]]) == 0) params <- params[params != p]
    }
    dt <- as_tibble(mat)
    out <- t(combn(params, 2))
    output <- map(1:nrow(out), \(i) {
        cols <- out[i, ]
        dt[1:N, cols] |>
            set_names(c("xvar", "yvar")) |>
            mutate(x = cols[1], y = cols[2], output = Y[1:N])
    }) |>
        list_rbind()

    y_range <- output$output |> range() |> floor() |> (\(x) x + c(0, 1))()

    plot_list <- output |>
        mutate(across(x:y, \(x) factor(pretty_params(x),
                                       levels = pretty_params(params)))) |>
        split(~ x + y, drop = TRUE) |>
        map(\(z) {
            z |>
                ggplot(aes(xvar, yvar, color = output)) +
                geom_point(size = 0.5) +
                scale_colour_gradientn(colours = grDevices::terrain.colors(10),
                                       limits = y_range, name = "y") +
                scale_x_continuous(z$x[[1]], breaks = scales::pretty_breaks(n = 3)) +
                scale_y_continuous(z$y[[1]], breaks = scales::pretty_breaks(n = 3)) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.title.x = element_markdown(),
                      axis.title.y = element_markdown(),
                      legend.background = element_rect(fill = "transparent", color = NA),
                      legend.key = element_rect(fill = "transparent", color = NA))
        })

    do.call(wrap_plots, c(plot_list, list(guides = "collect"))) &
        theme(legend.position = "top")
}






# plot_uncertainty(Y = y, N = N) + labs(y = "Counts", x = "$y$")


plot_multiscatter(data = mat, N = N, Y = y, params = params)

ind <- sobol_indices(Y = y, N = N, params = params, boot = TRUE, R = R,
                     type = type, conf = conf)

ind
