
library(tidyverse)
library(aeonia)
library(progressr)
library(ggtext)


.n_threads <- max(1L, parallel::detectCores()-2)


#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists("LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}


.B <- 0.1
# Carrying capacity with no B or natural enemies:
CC <- with(list(L = rbind(c(pop_info$surv_j, pop_info$fecund),
                          c(pop_info$recruit, pop_info$surv_a)),
                K = pop_info$K),
           {K * (max(abs(eigen(L)[["values"]])) - 1)})

P = 25
a = 5e-3
h = 0.1
k = 1
z = 100

curve((1 + a * P / (k * (x * z + 1)))^(-k), 0, 1)


# insect_sims <-
test_insect_pops(max_sim_t = 1000,
                                B = .B,
                                a = 5e-3,
                                h = 0.02,
                                k = 1,
                                s = 0.85,
                                # alate_0 = -Inf,
                                # alate_1 = 0,
                                A0 = 100,
                                W0 = 0,
                                P0 = 1) |>
    mutate(aphids = aphids + alates) |>
    select(-alates) |>
    pivot_longer(aphids:preds, names_to = "species", values_to = "N") |>
    mutate(species = factor(species, levels = c("aphids", "preds"))) |>
    # filter(species != "alates") |>
    ggplot(aes(time, N, color = species)) +
    geom_hline(yintercept = CC, linetype = "22", color = "gray70", linewidth = 1) +
    geom_line(linewidth = 1) +
    scale_color_viridis_d(begin = 0.2, end = 0.9) +
    scale_y_continuous("Abundance", limits = c(0, CC)) +
    theme_minimal()




# Number of replicate simulations:
# n_sims <- 100L
n_sims <- 12L





#' Vary:
#' - Pseudomonas abundance
#' - effect of Pseudomonas on abundance (`B` via `give_B()`)
#' - carrying capacity (`K`)


one_sim_combo <- function(pseudo, B, K, alpha, beta, epsilon, ...) {
    insect_ptr <- make_insect_ptr(K = K,
                                  B = B,
                                  a = 5e-3,
                                  h = 0.008,
                                  k = 0.1,
                                  s = 0.1,
                                  fly_p = 0.05)
    # land <- matrix(c(rep(0L, 4), 1L, rep(0L, 4)), 3, 3)
    # if (pseudo >= 2) land[1,1] <- land[3,3] <- 2L
    # if (pseudo >= 4) land[1,3] <- land[3,1] <- 2L
    others <- rep(0L, 10e3-1)
    if (pseudo > 0) others[sample.int(10e3-1, pseudo*1e3)] <- 2L
    land <- matrix(c(1L, others), 100, 100)
    land <- array(land, c(dim(land), n_sims))
    m0 <- matrix(0, nrow(land), ncol(land)) # for abundances
    out <- sim_plantscape(land,
                          max_sim_t = 80,
                          insect_ptr,
                          A0 = m0+100,
                          W0 = m0,
                          P0 = m0,
                          alpha = alpha,
                          beta = beta,
                          epsilon = epsilon,
                          delta_a = 0.5,
                          delta_p = 0.5,
                          radius = 1.5,
                          n_threads = .n_threads,
                          out_by_plant = FALSE,
                          ...) |>
        mutate(pseudo = pseudo, B = B, K = K, alpha = alpha,
               beta = beta, epsilon = epsilon)
    return(out)
}



# # Takes ~24 min on a 100x100 landscape (1.5 min for 3x3 landscape)
# ps_sims <- crossing(pseudo = c(0, 2, 4),
#                     B = c(0.1, 0.05, 0.01, 0),
#                     # K = 12500 * (-1:1 * 0.25 + 1),
#                     K = 12500,
#                     alpha = c(0, 0.25, 1, 2),
#                     beta = -1 * c(0, 0.25, 1, 2),
#                     epsilon = c(0.25, 1, 2)) |>
#     (function(x){
#         x |>
#             pmap(one_sim_combo, .progress = TRUE)
#     })()

# write_rds(ps_sims, "_building/ps_sims.rds", compress = "xz", compression = 9L)
ps_sims <- read_rds("_building/ps_sims.rds")



# How long does it take to reach fully infected?
ps_tti <- ps_sims |>
    map(\(x) {
        x |>
            group_by(pseudo, B, K, alpha, beta, epsilon, rep) |>
            summarize(time = (\(t,v) {
                if (any(v == 9)) return(t[v == 9][1])
                return(Inf)
            })(time, virus),
            .groups = "drop")
    }) |>
    list_rbind() |>
    mutate(across(pseudo:rep, factor))

ps_tti |>
    filter(epsilon == 1, K == 12500) |>
    ggplot(aes(pseudo, log10(time))) +
    geom_hline(aes(yintercept = min(log10(.data[["time"]]))),
               linewidth = 0.75, color = "black") +
    # geom_hline(aes(yintercept = max(log10(.data[["time"]]))),
    #            linewidth = 1, color = "gray70") +
    geom_jitter(aes(color = B), alpha = 0.25, width = 0.2, height = 0) +
    stat_summary(fun.data = "mean_cl_boot") +
    # ggplot(aes(time)) +
    # geom_freqpoly(aes(color = pseudo), bins = 10) +
    scale_color_viridis_d("*Pseudomonas*-<br>induced<br>mortality",
                          option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ylab("log<sub>10</sub>(Time to fully infected)") +
    xlab("Number of *Pseudomonas* patches") +
    # facet_wrap(~ alpha + beta, labeller = \(x) label_both(x, FALSE, "=")) +
    facet_grid(alpha ~ beta, labeller = \(x) label_both(x, FALSE, "=")) +
    theme_minimal() +
    theme(strip.text = element_text(size = 8),
          axis.title = element_markdown(),
          legend.title = element_markdown())





ps_alates <- ps_sims |>
    map(\(x) {
        x |>
            group_by(pseudo, B, K, alpha, beta, epsilon, rep) |>
            summarize(log_alates = mean(log10(alates+1)),
                      alates = mean(alates),
                      .groups = "drop")
    }) |>
    list_rbind() |>
    mutate(across(pseudo:rep, factor))



ps_alates |>
    filter(epsilon == 1, K == 12500) |>
    ggplot(aes(pseudo, log_alates)) +
    geom_hline(aes(yintercept = min(.data[["log_alates"]])),
               linewidth = 0.75, color = "black") +
    # geom_hline(aes(yintercept = max(.data[["log_alates"]])),
    #            linewidth = 1, color = "gray70") +
    geom_jitter(aes(color = B), alpha = 0.25, width = 0.2, height = 0) +
    # stat_summary(aes(color = B), fun.data = "mean_cl_boot") +
    # ggplot(aes(time)) +
    # geom_freqpoly(aes(color = pseudo), bins = 10) +
    scale_color_viridis_d("*Pseudomonas*-<br>induced<br>mortality",
                          option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ylab("Mean log<sub>10</sub>(total alates)") +
    xlab("Number of *Pseudomonas* patches") +
    # facet_wrap(~ alpha + beta, labeller = \(x) label_both(x, FALSE, "=")) +
    facet_grid(alpha ~ beta, labeller = \(x) label_both(x, FALSE, "=")) +
    theme_minimal() +
    theme(strip.text = element_text(size = 8),
          axis.title = element_markdown(),
          legend.title = element_markdown())






# =============================================================================*
# =============================================================================*
# OLD CODE ========================
# =============================================================================*
# =============================================================================*


ps_sims |>
    mutate(rep = factor(rep), pseudo = factor(pseudo)) |>
    group_by(pseudo, rep, time) |>
    summarize(virus = sum(virus), .groups = "drop") |>
    mutate(id = interaction(rep, pseudo, drop = TRUE)) |>
    ggplot(aes(time, virus)) +
    geom_hline(yintercept = c(0, length(land0[,,1])),
               linetype = "22", color = "gray70") +
    geom_line(aes(group = id, color = pseudo), alpha = 0.25) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    # facet_wrap(~ pseudo, ncol = 1) +
    theme_minimal()

ps_sims |>
    mutate(rep = factor(rep), pseudo = factor(pseudo)) |>
    group_by(pseudo, rep, time) |>
    summarize(virus = sum(virus), .groups = "drop") |>
    group_by(pseudo, rep) |>
    summarize(time = (\(t,v) {
        if (any(v == 9)) return(t[v == 9][1])
        return(Inf)
    })(time, virus),
              .groups = "drop") |>
    ggplot(aes(pseudo, time)) +
    geom_jitter(aes(color = pseudo), alpha = 0.25, width = 0.2, height = 0) +
    stat_summary(fun.data = "mean_cl_boot") +
    # ggplot(aes(time)) +
    # geom_freqpoly(aes(color = pseudo), bins = 10) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ylab("Time to fully infected") +
    theme_minimal()

ps_sims |>
    filter(rep == 1) |>
    select(-virus) |>
    pivot_longer(aphids:preds, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("aphids", "alates", "preds")),
           id = interaction(pseudo, rep, type, x, y, drop = TRUE)) |>
    ggplot(aes(time, density, color = type)) +
    geom_line(aes(group = id), alpha = 0.1) +
    scale_color_viridis_d(begin = 0.2, end = 0.9) +
    theme_minimal() +
    facet_wrap(~ type, scales = "free_y")

ps_sims |>
    filter(rep == 1) |>
    pivot_longer(virus:preds, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("virus", "aphids", "alates", "preds")),
           plant = interaction(x, y, drop = TRUE),
           id = interaction(pseudo, rep, type, plant, drop = TRUE)) |>
    filter(type != "preds", type != "virus") |>
    ggplot(aes(time, density, color = type)) +
    geom_line(aes(group = id)) +
    scale_color_viridis_d(begin = 0.2, end = 0.9) +
    theme_minimal() +
    facet_grid(pseudo ~ plant)





Rcpp::sourceCpp(code =
"#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// Retrieve dataset named `dataset` of type `T` from this package
template <typename T>
T retrieve_dataset(const std::string& dataset) {

    Environment pkg(\"package:aeonia\");
    T data = pkg[dataset];
    return data;
}

//[[Rcpp::export]]
double foo(double surv_j = NA_REAL) {
    List pop_info = retrieve_dataset<List>(\"pop_info\");
    if (!arma::is_finite(surv_j)) surv_j = pop_info[\"surv_j\"];
    return surv_j;
}
"
)
