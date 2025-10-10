
#'
#' Small-scale sensitivity via Sobol indices
#'
#'


library(sensobol)

source("_scripts/01-sensitivity/00-preamble.R")


#' Directory containing output from `sobol-sens-paired.R`:
sobol_dir <- "~/_globus"



struct_pars <- crossing(mu_Y = log(c(0.1, 1)),
                        mu_N = log(c(2, 4)),
                        n_pseudo = 2L)

vary_pars <- list(B = c(0, 0.15), # >= 0.1655339 results in carrying capacity of ~0
                  K = 12500 * c(0.1, 2),
                  alpha = c(0, 5),
                  beta = c(-5, 0),
                  wasp_disp_m0 = 0.3 * c(0, 2),
                  wasp_disp_m1 = 0.349 * c(0, 2))

N <- 2^8

# Takes just a second or two
set.seed(641272456)
sobol_mats <- map(1:nrow(struct_pars), \(i) {
    mat <- sobol_matrices(N = N, params = names(vary_pars))
    for (n in names(vary_pars)) {
        p1 <- min(vary_pars[[n]])
        p2 <- max(vary_pars[[n]])
        mat[,n] <- p1 + mat[,n] * (p2 - p1)
    }
    return(mat)
})


#' Should be TRUE
sobol_mats |>
    map_lgl(\(x) identical(sort(names(vary_pars)), sort(colnames(x)))) |>
    all()


#' Output from `sobol-sens.R`:
sobol_sims <- paste0(sobol_dir, "/sobol-sims-paired.rds") |>
    read_rds()



scatter <- function(sim_outs, .struct_pars, yvar = "infect_time") {

    if ("K" %in% colnames(sim_outs)) sim_outs[["K"]] <- sim_outs[["K"]] / 1000

    y <- sim_outs[[yvar]]
    params <- names(vary_pars)

    .title <- sprintf("%s = %.1f", pretty_params(names(.struct_pars)[1:2]),
                      as.numeric(.struct_pars[1:2])) |>
        paste(collapse = "; ")

    sim_outs |>
        rename(y = infect_time) |>
        select(all_of(params), y) |>
        pivot_longer(all_of(params), names_to = "param", values_to = "x") |>
        mutate(param = param |>
                   pretty_params() |>
                   factor(levels = pretty_params(params))) |>
        ggplot(aes(x, y)) +
        geom_point(size = 0.7, alpha = 0.1) +
        geom_hline(yintercept = 0, linetype = "22", color = "firebrick",
                   linewidth = 1) +
        stat_smooth(method = "gam", formula = y ~ s(x), se = FALSE,
                    linewidth = 1, color = "dodgerblue") +
        facet_wrap(~ param, scales = "free_x") +
        labs(title = .title,
             x = "Value",
             y = paste0("Effect of *Pseudomonas* on days to 5 plants infected<br>",
                        "(< 0 means *Pseudomonas* speeds virus spread)")) +
        theme(strip.text = element_markdown(),
              axis.title.y = element_markdown(),
              plot.title = element_markdown())
}



sobol_plots <- crossing(j = 1:nrow(struct_pars),
                        p = c("infect_time", "infect_time_Inf")) |>
    pmap(\(j, p) scatter(sobol_sims[[j]], struct_pars[j,], p) )

# for (j in 1:length(sobol_plots)) {
#     save_plot(sprintf("_plots/sobol-paired/sobol-paired-%i.pdf", j), sobol_plots[[j]],
#               width = 6, height = 6)
#     save_plot(sprintf("_plots/sobol-paired/sobol-paired-%i.svg", j), sobol_plots[[j]],
#               width = 6, height = 6)
#     save_plot(sprintf("_plots/sobol-paired/sobol-paired-simple-%i.png", j),
#               sobol_plots[[j]] + theme(axis.title.y = element_blank(),
#                                        plot.title = element_blank()),
#               width = 6, height = 6, dpi = 150)
#     cat(j, "\n")
# }; rm(j)


library(survival)

x <- rlnorm(1e6, 1, 1)

# hist(x)

xt <- quantile(x, 0.8) |> unname()
z <- x <= xt
xx <- x
xx[xx > xt] <- xt

xs <- Surv(xx, event = z, type = "right")
xsr <- survreg(xs ~ 1, dist = "lognormal")
predict(xsr, newdata = data.frame(zzz = 1)) |> unname(); mean(x)


plot(survfit(xs~1),xlim=c(0,xt))


crossing(x = 1:100, h = c(0.5, 1:5)) |>
    group_by(h) |>
    mutate(z = x^h / sum(x^h),
           z2 = 1- z) |>
    ungroup() |>
    mutate(h = factor(h)) |>
    ggplot(aes(x, z2, color = h)) +
    geom_point() +
    scale_color_viridis_d(option = "plasma", end = 0.95, begin = 0.1)



curve((1 + pop_info$a * 10 / (pop_info$k * (pop_info$h * x + 1)))^(- 1 ), 0, 1000)





multi_scatter <- function(sim_outs, mat, .struct_pars, .N = N, .params = NULL) {

    if ("K" %in% colnames(mat)) mat[,"K"] <- mat[,"K"] / 1000

    Y <- sim_outs[["infect_time"]]
    if (is.null(.params)) {
        params <- names(vary_pars)
        for (p in params) {
            if (sd(sim_outs[[p]]) == 0) params <- params[params != p]
        }
    } else params <- .params
    dt <- as_tibble(mat)
    out <- t(combn(params, 2))
    output <- map(1:nrow(out), \(i) {
        cols <- out[i, ]
        dt[1:.N, cols] |>
            set_names(c("xvar", "yvar")) |>
            mutate(x = cols[1], y = cols[2], output = Y[1:.N])
    }) |>
        list_rbind()

    y_range <- output$output |> range() |> floor() |> (\(x) x + c(0, 1))()

    .title <- sprintf("%s = %.0f", pretty_params(names(.struct_pars)),
                      as.numeric(.struct_pars)) |>
        paste(collapse = "; ")

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

    do.call(wrap_plots, c(plot_list, list(guides = "collect"))) +
        plot_annotation(title = .title,
                        theme = theme(plot.title = element_markdown())) &
        theme(legend.position = "top")
}



# j = 5
# multi_scatter(sobol_sims[[j]], sobol_mats[[j]], struct_pars[j,],
#               .params = c("K", "alpha", "beta", "wasp_disp_m0", "wasp_disp_m1"))






one_combo <- function(mu_Y, mu_N, n_pseudo,
                      B, K, alpha, beta, wasp_disp_m0, wasp_disp_m1,
                      n_sims, ...) {

    fly_p <- 0.05
    sigma_Y <- 1
    sigma_N <- 1
    epsilon <- 1

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
                       infect_stop = TRUE,
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


# Rcpp::cppFunction("IntegerVector foo() {
#     arma::vec out(2);
#     out(0) = arma::datum::nan;
#     out(1) = 1.0;
#     IntegerVector out2 = wrap(out);
#     return out2;
# }", depends = "RcppArmadillo", includes = "#include <RcppArmadillo.h>")
#
# foo()


j = 7
sobol_sims[[j]] |>
    arrange(infect_time) |>
    select(infect_time, everything())

set.seed(45678)
ss <- sobol_sims[[j]] |>
    arrange(infect_time) |>
    slice_head(n = 100) |>
    slice_sample(n = 1)

.mu_Y = 2
.mu_N = 4
.B = ss[["B"]][[1]]  # 0.144
.K = ss[["K"]][[1]]  # 24e3
.alpha = ss[["alpha"]][[1]]  # 1.16
.beta = ss[["beta"]][[1]]  # -0.52
.wasp_disp_m0 = ss[["wasp_disp_m0"]][[1]]  # 0.146
.wasp_disp_m1 = ss[["wasp_disp_m1"]][[1]]  # 0.0395


sims0 <- one_combo(mu_Y = .mu_Y, mu_N = .mu_N, n_pseudo = 0, B = .B,
                   K = .K, alpha = .alpha, beta = .beta,
                   wasp_disp_m0 = .wasp_disp_m0, wasp_disp_m1 = .wasp_disp_m1,
                   max_t = 1e3,
                   n_sims = 1e3, summ = "all")
sims <- one_combo(mu_Y = .mu_Y, mu_N = .mu_N, n_pseudo = 2, B = .B,
                  K = .K, alpha = .alpha, beta = .beta,
                  wasp_disp_m0 = .wasp_disp_m0, wasp_disp_m1 = .wasp_disp_m1,
                  max_t = 1e3,
                  n_sims = 1e3, summ = "all")

# If first is greater than second, then Pseudomonas is bad for plants:
sims0$infect_time |> mean(na.rm = TRUE)
sims$infect_time |> mean(na.rm = TRUE)

sum(is.na(sims0$infect_time))
sum(is.na(sims$infect_time))





sims_ts <- one_combo(mu_Y = .mu_Y, mu_N = .mu_N, n_pseudo = 2, B = .B,
                     K = .K, alpha = .alpha, beta = .beta,
                     wasp_disp_m0 = .wasp_disp_m0, wasp_disp_m1 = .wasp_disp_m1,
                     n_sims = 6, summ = "time")
sims_ts |>
    mutate(rep = factor(rep)) |>
    select(rep:wasps) |>
    pivot_longer(virus:wasps) |>
    mutate(id = interaction(rep, name, drop = TRUE)) |>
    ggplot(aes(time, value, group = id)) +
    geom_line(aes(color = name)) +
    facet_wrap(~ rep, nrow = 2, scales = "free") +
    scale_color_viridis_d(begin = 0.2, end = 0.95)



dd <- sims_ts |>
    filter(rep == 2) |>
    select(time:wasps) |>
    mutate(aphids = aphids + alates + parasitized) |>
    select(-alates, -parasitized, -mummies) |>
    pivot_longer(virus:wasps) |>
    mutate(plant = interaction(x, y, drop = TRUE),
           is_virus = factor(name == "virus"))

trans <- \(x) x * max(dd$value) / max(x)

# LEFT OFF HERE ----

dd |>
    group_by(name) |>
    mutate(value = ifelse(name == "aphids", value, )) |>
    ungroup() |>
    ggplot(aes(time, value)) +
    geom_line(aes(color = name, linetype = is_virus), linewidth = 1) +
    facet_wrap(~ plant, nrow = 3) +
    scale_color_viridis_d(begin = 0.2, end = 0.95, option = "plasma") +
    scale_linetype_manual(values = c(`TRUE` = "22", `FALSE` = "solid"))


