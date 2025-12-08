
#'
#' Basic time series simulations
#'

source("_scripts/00-preamble.R")
source("_scripts/01-basic-sims/00-time-series-plotter.R")









n_sims <- 12L

landscapes <- array(0L, c(3L, 3L, n_sims))
for (i in 1:(dim(landscapes)[3])) {
    landscapes[1,1,i] <- 1L
    # if (n_pseudo > 0) {
    #     k <- sample.int(n_plants - 1L, n_pseudo)
    #     x <- k - n_x * (k %/% n_x) + 1L
    #     y <- k %/% n_x + 1L
    #     land[cbind(x,y,i)] <- 2L
    # }
    landscapes[1,3,i] <- landscapes[2,2,i] <- landscapes[3,1,i] <- 2L
}


zeta <- 1
xi <- 1
psi <- 0.9
nu <- 1
rho <- 1
mean_N <- 5
sd_N  <- 5
Y0 <- 1
# p_title <- sprintf(paste(pretty_params(c("zeta", "virus_attract", "pseudo_surv",
#                                          "pseudo_repel")),
#                          c("= %.2f&emsp;&emsp; | &emsp;", "= %.1f<br>",
#                            "= %.2f&emsp; | &emsp;", "= %.1f"), collapse = ""),
#                    zeta, nu, psi, rho)
p_title <- sprintf(paste(pretty_params(c("zeta", "mean_N", "sd_N", "Y0")),
                         c("= %.2f&emsp; | &emsp;", "= %.1f<br>",
                           "= %.2f&emsp; | &emsp;", "= %.1f"), collapse = ""),
                   zeta, mean_N, sd_N, Y0)



sims <- sim_plantscape(landscapes = landscapes, max_t = 100L,
                       insect_ptr = make_insect_ptr(pseudo_surv = psi, fly_p = 0.05, zeta = zeta),
                       N0 = array(rlnorm(prod(dim(landscapes)),
                                         log(mean_N^2 / sqrt(mean_N^2 + sd_N^2)),
                                         sqrt(log(1 + sd_N^2 / mean_N^2))),
                                  dim(landscapes)),
                       Y0 = Y0,
                       W0 = array(0, dim(landscapes)),
                       virus_attract = nu,
                       pseudo_repel = rho,
                       epsilon = 1,
                       p_load_alate = 0.5,
                       p_load_plant = 0.5,
                       radius = 1,
                       infect_stop = FALSE)


sim_plotter(sims, .title = p_title) +
    facet_wrap( ~ rep, nrow = 3) +
    theme(strip.text = element_blank(),
          strip.background = element_blank(),
          legend.position = "none")


# sims |> sim_plotter(TRUE, .title = p_title)


#
#
#
#
# a <- pop_info$a
# h <- pop_info$h
# k <- pop_info$k
#
# crossing(x = (0:100) * 10, zh = seq(0, 1, length.out = 101),
#          zeta = c(0, 1)) |>
#     mutate(gamma_i = 100 * ((1-zeta) * 1/9 + zeta * zh),
#            A = (1 + a / (k * (h * x + 1)) * gamma_i)^(-k),
#            zeta = factor(zeta, levels = c(0, 1),
#                          labels = c("uniform", "&prop; density"))) |>
#     ggplot(aes(x, zh)) +
#     geom_raster(aes(fill = A)) +
#     geom_contour(aes(z = A), color = "gray70") +
#     scale_fill_viridis_c() +
#     facet_wrap(~ zeta) +
#     theme(strip.text = element_markdown(family = "STIX Two Math", size = 12))
