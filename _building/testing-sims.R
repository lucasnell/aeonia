library(tidyverse)
library(aeonia)
library(future.apply)
library(progressr)
library(ggtext)

# Set threads for simulations:
options("mc.cores" = max(1L, parallel::detectCores()-2))

plan(multisession, workers = options()[["mc.cores"]])
handlers(global = TRUE)
handlers("progress")


#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists("LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}


# Carrying capacity with no B or natural enemies:
CC <- with(list(L = rbind(c(pop_info$surv_j, pop_info$fecund),
                          c(pop_info$recruit, pop_info$surv_a)),
                K = pop_info$K),
           {K * (max(abs(eigen(L)[["values"]])) - 1)})


# # insect_sims <-
# test_insect_pops(max_t = 1000,
#                                 B = 0.1,
#                                 a = 5e-3,
#                                 h = 0.02,
#                                 k = 1,
#                                 s = 0.85,
#                                 # alate_0 = -Inf,
#                                 # alate_1 = 0,
#                                 A0 = 100,
#                                 W0 = 0,
#                                 P0 = 1) |>
#     mutate(aphids = aphids + alates) |>
#     select(-alates) |>
#     pivot_longer(aphids:preds, names_to = "species", values_to = "N") |>
#     mutate(species = factor(species, levels = c("aphids", "preds"))) |>
#     # filter(species != "alates") |>
#     ggplot(aes(time, N, color = species)) +
#     geom_hline(yintercept = CC, linetype = "22", color = "gray70", linewidth = 1) +
#     geom_line(linewidth = 1) +
#     scale_color_viridis_d(begin = 0.2, end = 0.9) +
#     scale_y_continuous("Abundance", limits = c(0, CC)) +
#     theme_minimal()






#' Vary:
#' - Pseudomonas abundance
#' - effect of Pseudomonas on abundance (`B` via `give_B()`)
#' - carrying capacity (`K`)


one_sim_combo <- function(pseudo, B, K, alpha, beta, epsilon,
                          max_t, n_x, n_y, radius, n_sims,
                          .full_inf_time = TRUE) {

    insect_ptr <- make_insect_ptr(K = K,
                                  B = B,
                                  a = 5e-3,
                                  h = 0.008,
                                  k = 0.1,
                                  s = 0.1,
                                  fly_p = 0.05)
    stopifnot(round(n_x) == n_x && round(n_y) == n_y && round(n_sims) == n_sims)
    n_x <- as.integer(n_x)
    n_y <- as.integer(n_y)
    n_plants <- n_x * n_y
    stopifnot(round(pseudo) == pseudo && pseudo >= 0 && pseudo <= (n_plants-1L))

    land <- array(0L, c(n_x, n_y, n_sims))
    for (i in 1:n_sims) {
        land[1,1,i] <- 1L
        if (pseudo > 0) {
            k <- sample.int(n_plants - 1L, pseudo)
            x <- k - n_x * (k %/% n_x) + 1L
            y <- k %/% n_x + 1L
            land[cbind(x,y,i)] <- 2L
        }
    }
    # Define function for infection summary:
    if (.full_inf_time) {
        # Time to full infection:
        inf_summ <- function(t, v, p) {
            v <- v[!p] + v[p]
            t <- t[!p]
            if (any(v == n_plants)) return(t[v == n_plants][1])
            return(Inf)
        }
        # Same but for only non-Pseudomonas plants:
        inf_summ2 <- function(t, v, p) {
            t <- t[!p]
            v <- v[!p]
            np <- n_plants - pseudo
            if (any(v == np)) return(t[v == np][1])
            return(Inf)
        }
    } else {
        # Number of plants infected (i.e., outbreak size):
        inf_summ <- function(t, v, p) {
            v <- v[!p] + v[p]
            return(max(v))
        }
        # Same but for only non-Pseudomonas plants:
        inf_summ2 <- function(t, v, p) {
            v <- v[!p]
            return(max(v))
        }
    }

    one_rep <- function(x) {
        xs <- x |> group_by(time) |> summarize(across(aphids:alates, sum))
        o <- tibble(rep = x$rep[[1]])
        o[["p_alates"]] = mean(xs$alates / (xs$alates + xs$aphids))
        o[["aphids"]] = mean(xs$alates + xs$aphids)
        o[["alates"]] = mean(xs$alates)
        o[["infect_time"]] = inf_summ(x$time, x$virus, x$pseudo)
        o[["infect_time2"]] = inf_summ2(x$time, x$virus, x$pseudo)
        return(o)
    }

    out <- sim_plantscape(land,
                          max_t = max_t,
                          insect_ptr,
                          A0 = matrix(10, n_x, n_y),
                          W0 = matrix(0, n_x, n_y),
                          P0 = matrix(0, n_x, n_y),
                          alpha = alpha,
                          beta = beta,
                          epsilon = epsilon,
                          delta_a = 0.5,
                          delta_p = 0.5,
                          radius = radius,
                          summ = "pseudo") |>
        split(~ rep) |>
        future_lapply(one_rep, future.seed = TRUE, future.packages = "dplyr") |>
        list_rbind() |>
        mutate(pseudo = pseudo, B = B, K = K, alpha = alpha,
               beta = beta, epsilon = epsilon)
    # if (.full_inf_time && any(is.infinite(out$infect_time))) {
    #     warning(sprintf(paste("Did not reach fully infected with",
    #                           "pseudo=%.2f, B=%.2f, K=%.2f, alpha=%.2f,",
    #                           "beta=%.2f, epsilon=%.2f"),
    #                     pseudo, B, K, alpha, beta, epsilon))
    # }
    if (!.full_inf_time) {
        out <- out |>
            rename(outbreak_size = infect_time,
                   outbreak_size2 = infect_time2)
    }
    return(out)
}






if (!file.exists("_building/ps_sims.rds")) {
    # Takes ~20 sec for 3x3 landscape
    set.seed(1114260777)
    ps_sims <- crossing(pseudo = c(0L, 3L),
                        B = c(0.1, 0.05, 0.01, 0),
                        K = 12500 * (-1:1 * 0.25 + 1),
                        alpha = c(0, 0.25, 1, 2),
                        beta = -1 * c(0, 0.25, 1, 2),
                        epsilon = c(0.25, 1, 2)) |>
        # Below, `max_t` is set to a value way higher than I think I'll need
        # so that all reps reach fully infected.
        mutate(n_x = 3L,
               n_y = ifelse(pseudo == 0L, 2L, 3L),
               radius = 1.5, max_t = 10e3, n_sims = 100) |>
        pmap(one_sim_combo, .progress = TRUE) |>
        list_rbind() |>
        select(pseudo, B, K, alpha, beta, epsilon, rep, everything()) |>
        mutate(across(pseudo:rep, factor))
    write_rds(ps_sims, "_building/ps_sims.rds", compress = "xz")
} else {
    ps_sims <- read_rds("_building/ps_sims.rds")
}





pretty_facet_factors <- function(.df, .pars, .greek) {

    if (length(.greek) == 1) .greek <- rep(.greek, length(.pars))
    stopifnot(length(.pars) == length(.greek))
    stopifnot(is.logical(.greek))

    .greek <- as.list(.greek)
    names(.greek) <- .pars

    for (p in .pars) {
        if (is.factor(.df[[p]])) .df[[p]] <- as.numeric(paste(.df[[p]]))
        .lvls <- sort(unique(.df[[p]]))
        .name <- ifelse(.greek[[p]], paste0("&", p, ";"), p)
        .fmt <- ifelse(all(.lvls == round(.lvls)), "%i", "%.2f")
        .labs <- sprintf(paste0(.name, " = ", .fmt), .lvls) |>
            str_replace_all("-", "&minus;")
        .df[[p]] <- factor(.df[[p]], levels = sort(unique(.df[[p]])),
                           labels = .labs)
    }

    return(.df)

}








ps_sims |>
    filter(epsilon == 1,
           K == levels(K)[2]) |>
           # alpha == levels(alpha)[1]) |>
    pretty_facet_factors(c("alpha", "beta", "epsilon", "K"),
                         .greek = c(rep(TRUE, 3), FALSE)) |>
    ggplot(aes(B, (infect_time2), color = pseudo)) +
    geom_hline(aes(yintercept = min((.data[["infect_time2"]]))),
               linewidth = 0.75, color = "black") +
    geom_violin(position = position_dodge(0.5), fill = NA) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0.5)) +
    scale_color_viridis_d("Proportion of<br>*Pseudomonas*<br>patches",
                          option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ylab("Days to fully infected") +
    xlab("*Pseudomonas*-induced mortality") +
    facet_grid(alpha ~ beta, scales = "free_y") +
    # facet_grid(K ~ beta) +
    theme_minimal() +
    theme(strip.text = element_markdown(size = 8),
          axis.title = element_markdown(),
          legend.title = element_markdown())



ps_sims |>
    filter(epsilon == 1,
           K == median(as.numeric(paste(K)))) |>
           # alpha == levels(alpha)[2]) |>
    pretty_facet_factors(c("alpha", "beta", "epsilon", "K"),
                         .greek = c(rep(TRUE, 3), FALSE)) |>
    ggplot(aes(B, p_alates, color = pseudo)) +
    geom_hline(aes(yintercept = min(.data[["p_alates"]])),
               linewidth = 0.75, color = "black") +
    # geom_hline(aes(yintercept = max(.data[["log_alates"]])),
    #            linewidth = 1, color = "gray70") +
    geom_violin(position = position_dodge(0.5), fill = NA) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0.5)) +
    scale_color_viridis_d("Proportion of<br>*Pseudomonas*<br>patches",
                          option = "plasma", end = 0.8) +
    ylab("Mean proportion alates") +
    xlab("*Pseudomonas*-induced mortality") +
    facet_grid(alpha ~ beta) +
    # facet_grid(K ~ beta) +
    theme_minimal() +
    theme(strip.text = element_markdown(size = 8),
          axis.title = element_markdown(),
          legend.title = element_markdown())



# ---------------------------*
# How often does Pseudomonas increase time to fully infected? ----
# ---------------------------*


ps_inc_sims <- ps_sims |>
    split(~ B + K + alpha + beta + epsilon) |>
    lapply(\(d) {
        .it0 <- d$infect_time[d$pseudo == 0]
        .it20 <- d$infect_time2[d$pseudo == 0]
        .pa0 <- d$p_alates[d$pseudo == 0]
        d <- d[d$pseudo != 0,]
        d[["p_better"]] <- sapply(d$infect_time, \(it) mean(.it0 < it))
        d[["p_better2"]] <- sapply(d$infect_time2, \(it2) mean(.it20 < it2))
        d[["p_more"]] <- sapply(d$p_alates, \(pa) mean(.pa0 < pa))
        return(d)
    }) |>
    list_rbind() |>
    mutate(pseudo = fct_drop(pseudo))


ps_inc_sims |>
    # filter(pseudo == levels(pseudo)[2]) |>
    filter(beta %in% levels(beta)[c(1,4)],
           alpha %in% levels(alpha)[c(1,4)]) |>
    filter(# epsilon == 0.25,
           K == levels(K)[2]) |>
    # alpha == levels(alpha)[1]) |>
    # pretty_facet_factors(c("alpha", "beta", "epsilon", "K"),
                         # .greek = c(rep(TRUE, 3), FALSE)) |>
    pretty_facet_factors(c("alpha", "beta", "K"),
                         .greek = c(rep(TRUE, 2), FALSE)) |>
    ggplot(aes(B, p_better2 * 100, color = epsilon)) +
    geom_hline(aes(yintercept = 50),
               linewidth = 0.75, color = "gray70", linetype = "22") +
    geom_violin(position = position_dodge(0.6), fill = NA) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0.6), size = 3) +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",
                 position = position_dodge(0.6), width = 0.3) +
    scale_color_viridis_d("&epsilon;",
                          option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ylab("Percent longer to infected than without *Pseudomonas*") +
    xlab("*Pseudomonas*-induced mortality") +
    facet_grid(alpha ~ beta) +
    # facet_grid(K ~ beta) +
    theme_minimal() +
    theme(strip.text = element_markdown(size = 8),
          axis.title = element_markdown(),
          legend.title = element_markdown())



ps_inc_sims |>
    filter(beta %in% levels(beta)[c(1,4)],
           alpha %in% levels(alpha)[c(1,4)]) |>
    filter(pseudo == levels(pseudo)[1],
           K == levels(K)[2]) |>
    # alpha == levels(alpha)[1]) |>
    # pretty_facet_factors(c("alpha", "beta", "epsilon", "K"),
                         # .greek = c(rep(TRUE, 3), FALSE)) |>
    pretty_facet_factors(c("alpha", "beta"), .greek = TRUE) |>
    (\(.data) {
        .title <<- sprintf("%i *Pseudomonas* patches<br>K = %s",
                           as.integer(paste(.data[["pseudo"]][[1]])),
                           .data[["K"]][[1]])
        return(.data)
    })() |>
    ggplot(aes(B, p_more * 100, color = epsilon)) +
    geom_hline(aes(yintercept = 50),
               linewidth = 0.75, color = "gray70", linetype = "22") +
    geom_violin(position = position_dodge(0.6), fill = NA) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0.6), size = 3) +
    stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",
                 position = position_dodge(0.6), width = 0.3) +
    scale_color_viridis_d("&epsilon;",
                          option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ggtitle(.title) +
    ylab("Percent more alates than without *Pseudomonas*") +
    xlab("*Pseudomonas*-induced mortality") +
    facet_grid(alpha ~ beta) +
    # facet_grid(K ~ beta) +
    theme_minimal() +
    theme(strip.text = element_markdown(size = 8),
          axis.title = element_markdown(),
          legend.title = element_markdown(),
          plot.title = element_markdown())





# =============================================================================*
# =============================================================================*
# LARGER SIMULATIONS ========================
# =============================================================================*
# =============================================================================*

if (!file.exists("_building/big_ps_sims.rds")) {
    # Takes ~8 hrs for 134x134 landscape (and 12 instead of 100 sims)
    t0 <- Sys.time()
    set.seed(1952926471)
    big_ps_sims <- crossing(pseudo = c(0.0, 0.2, 0.5),
                            B = c(0.1, 0.05, 0.01, 0),
                            K = 12500 * (-1:1 * 0.25 + 1),
                            alpha = c(0, 0.25, 1, 2),
                            beta = -1 * c(0, 0.25, 1, 2),
                            epsilon = c(0.25, 1, 2)) |>
        #' Below...
        #' `max_t` is set to a reasonable value for pea growing season.
        #' `n_x` and `n_y` approximate a square hectare with 0.75 m plant spacing
        mutate(n_x = 134L, n_y = 134L,
               radius = formals(sim_plantscape)[["radius"]],
               max_t = 100, n_sims = 12,
               # We are not outputting the time to full infection (all plants
               # infected). Instead we're returning outbreak size.
               .full_inf_time = FALSE) |>
        pmap(one_sim_combo, .progress = TRUE) |>
        list_rbind() |>
        select(pseudo, B, K, alpha, beta, epsilon, rep, everything()) |>
        mutate(pseudo = round(pseudo, digits = 1)) |>
        mutate(across(pseudo:rep, factor))
    # write_rds(big_ps_sims, "_building/big_ps_sims.rds", compress = "xz")
    t1 <- Sys.time()
    print(t1 - t0); # rm(t0, t1)

} else {
    big_ps_sims <- read_rds("_building/big_ps_sims.rds")
}


big_ps_sims |>
    filter(outbreak_size == max(outbreak_size))


big_ps_sims |>
    filter(epsilon == 2,
           # K == median(as.numeric(paste(K)))) |>
           alpha == levels(alpha)[1]) |>
    pretty_facet_factors(c("alpha", "beta", "epsilon", "K"),
                         .greek = c(rep(TRUE, 3), FALSE)) |>
    ggplot(aes(B, (outbreak_size), color = pseudo)) +
    # geom_hline(aes(yintercept = min((.data[["outbreak_size"]]))),
    #            linewidth = 0.75, color = "black") +
    geom_violin(position = position_dodge(0.5), fill = NA) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0.5)) +
    scale_color_viridis_d("Proportion of<br>*Pseudomonas*<br>patches",
                          option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ylab("Outbreak size") +
    xlab("*Pseudomonas*-induced mortality") +
    # facet_grid(alpha ~ beta, scales = "free_y") +
    facet_grid(K ~ beta, scales = "free_y") +
    theme_minimal() +
    theme(strip.text = element_markdown(size = 8),
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

class A {

    arma::vec X;

public:

    A(const arma::vec& X_) : X(X_) {}

    double& winged_adults() {
        return X.back();
    }

    arma::vec all() {
        return X;
    }

};


//[[Rcpp::export]]
arma::vec foo(const arma::vec& y, const double& z) {

    A x(y);
    x.winged_adults() += z;
    return x.all();
}

")
