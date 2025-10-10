

library(tidyverse)
library(aeonia)
library(future.apply)
library(progressr)
library(patchwork)
library(ggtext)

# Set threads for simulations:
options("mc.cores" = max(1L, parallel::detectCores()-2))

plan(multisession, workers = options()[["mc.cores"]])
handlers(global = TRUE)
handlers("progress")

# For purrr progress bar:
.prog_args <- list(clear = FALSE,
                   format = paste("{cli::pb_bar}",
                                  "{cli::pb_percent}",
                                  "[{cli::pb_elapsed}] |",
                                  "ETA: {cli::pb_eta}"))

# Transparent theme:
.trans_theme <- theme(panel.background = element_rect(fill="transparent"),
                      plot.background = element_rect(fill="transparent", color=NA),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      legend.background = element_rect(fill="transparent"),
                      legend.box.background = element_rect(fill="transparent"))

.no_axes <- theme(axis.title = element_blank(),
                  axis.text = element_blank(),
                  legend.position = "none",
                  strip.text = element_blank(),
                  plot.title = element_blank())

#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists(".LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}


# Carrying capacity with no alates or natural enemies:
CC <- function(K, B, m){
    L  <- rbind(c(pop_info$surv_j, pop_info$fecund),
                c(pop_info$recruit, pop_info$surv_a))
    K * (max(abs(eigen(L)[["values"]])) * (1 - B) * (1 - m) - 1)
}


spp_pal <- viridisLite::plasma(101)[c(10, 50, 80)] |>
    set_names(c("aphids", "alates", "enemies"))
pseudo_pal <- c(`3` = "#1E90FF", `0` = "gray60")


# # >> Insect pops ----
# p <- test_insect_pops(max_t = 100,
#                  B = 0,
#                  demog_error = FALSE,
#                  N0 = 10,
#                  W0 = 0,
#                  Y0 = 0) |>
#     select(-enemies) |>
#     pivot_longer(aphids:alates, names_to = "morph", values_to = "N") |>
#     mutate(morph = factor(morph, levels = c("aphids", "alates"))) |>
#     # filter(morph == "aphids") |>
#     ggplot(aes(time, N, color = morph)) +
#     # geom_hline(yintercept = CC(12500, 0, 0.1), color = "gray80", linewidth = 1) +
#     geom_line(linewidth = 1) +
#     # scale_color_viridis_d(begin = 0.2, end = 0.9, drop = FALSE) +
#     scale_color_manual(values = spp_pal, guide = "none") +
#     ylab("Abundance") +
#     coord_cartesian(ylim = c(0, 1600)) +
#     .no_axes +
#     .trans_theme
#
# ggsave("~/Desktop/p.svg", p, width = 3, height = 2, bg = "transparent")




make_arg_list <- function(n_pseudo, B, K, alpha, beta, epsilon, with_P,
                          max_t, n_x, n_y, radius, n_sims, ...) {


    stopifnot(round(n_x) == n_x && round(n_y) == n_y && round(n_sims) == n_sims)
    n_x <- as.integer(n_x)
    n_y <- as.integer(n_y)
    n_plants <- n_x * n_y
    stopifnot(round(n_pseudo) == n_pseudo && n_pseudo >= 0 && n_pseudo <= (n_plants-1L))

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

        # if (with_P) {
        #     idx <- sample.int(length(starts$N0), n_plants, replace = TRUE)
        #     N0[,,i] <- starts$N0[idx]
        #     Y0[,,i] <- starts$Y0[idx]
        # } else N0[,,i] <- sample(starts$N0, n_plants, replace = TRUE)\

        # N0[,,i] <- 10
        N0[,,i] <- exp(runif(n_plants, -3, 5))
        # N0[,,i] <- runif(n_plants, 0, 100)
        if (with_P) {
            # Y0[,,i] <- runif(n_plants, 0, 3)
            Y0[,,i] <- N0[,,i] * exp(runif(length(N0[,,i]), -10, -1))
        }

    }

    insect_args <- list(K = K, B = B, h = 5, fly_p = 0.05, wasp_disp_m0 = 0.3)
    plant_args <- list(landscapes = land,
                       max_t = max_t,
                       N0 = N0,
                       W0 = array(0.0, c(n_x, n_y, n_sims)),
                       Y0 = Y0,
                       alpha = alpha,
                       beta = beta,
                       epsilon = epsilon,
                       delta_a = 0.5,
                       delta_p = 0.5,
                       radius = radius,
                       infect_stop = FALSE,
                       summ = "all")
    other_args <- list(...)
    if (length(other_args) > 0) {
        stopifnot(!is.null(names(other_args)) && !any(names(other_args) == ""))
        stopifnot(all(names(other_args) %in% c(names(formals(sim_plantscape)),
                                               names(formals(make_insect_ptr)))))
        # not_allowed <- c("landscapes")
        # if (any(names(other_args) %in% not_allowed)) {
        #     not_allowed <- names(other_args)[names(other_args) %in% not_allowed]
        #     stop("The following are not allowed in `make_arg_list`: ",
        #          paste(not_allowed, collapse = ", "))
        # }

        nm_insect_args <- names(other_args)[names(other_args) %in%
                                                names(formals(make_insect_ptr))]
        for (n in nm_insect_args) insect_args[[n]] <- other_args[[n]]

        nm_plant_args <- names(other_args)[names(other_args) %in%
                                               names(formals(sim_plantscape))]
        for (n in nm_plant_args) {
            if (n %in% c("N0", "W0", "Y0", "landscape")) {
                plant_args[[n]] <- array(other_args[[n]], c(n_x, n_y, n_sims))
            } else plant_args[[n]] <- other_args[[n]]
        }
    }


    plant_args[["insect_ptr"]] <- do.call(make_insect_ptr, insect_args)

    return(plant_args)

}





one_sim_combo <- function(n_pseudo, B, K, alpha, beta, epsilon, with_P,
                          max_t, n_x, n_y, radius, n_sims, ...,
                          summarize = TRUE, doomsday = FALSE, p_N0 = 1,
                          time_inf_np = NULL) {

    # n_pseudo = 2; B = 0.1; K = 12500; alpha = 1; beta = -1; epsilon = 1;
    # with_P = TRUE; max_t = 100; n_x = 3; n_y = 3; radius = 1; n_sims = 100;
    # summarize = TRUE; doomsday = FALSE; p_N0 = 1; time_inf_np = NULL

    n_plants <- n_x * n_y
    if (p_N0 < 1 && doomsday) stop("Ambiguous: p_N0 < 1 && doomsday")
    if (doomsday) {
        .N0 <- c(10, rep(0, n_plants-1L))
    } else {
        # .N0 <- rep(10, n_plants * n_sims) # runif(n_plants * n_sims, 0, 1000)
        .N0 <- rep(10, n_plants * n_sims)  # exp(runif(n_plants * n_sims, -3, 5))
        if (p_N0 < 1) {
            for (i in 1:n_sims) {
                idx0 <- (i - 1L) * n_plants
                idx <- sample.int(n_plants, round(p_N0 * n_plants))
                .N0[idx0 + idx] <- 0
            }
        }
    }
    .Y0 <- 0
    if (with_P) .Y0 <- .N0 * exp(runif(n_plants * n_sims, -10, -1))

    arg_list <- make_arg_list(n_pseudo, B, K, alpha, beta, epsilon, with_P,
                              max_t, n_x, n_y, radius, n_sims,
                              N0 = .N0, Y0 = .Y0, ...)

    if (!summarize) arg_list[["summ"]] <- "none"


    # Time to `p` plants infected:
    time_inf_fun <- function(p, t, v) {
        out <- map_dbl(p, \(.p) {
            if (any(v >= .p)) return(t[v >= .p][1])
            return(Inf)
        })
        if (length(out) > 1) out <- list(tibble(np = p, infect_time = out))
        return(out)
    }

    if (is.null(time_inf_np)) time_inf_np <- n_plants

    out <- do.call(sim_plantscape, arg_list) |>
        mutate(n_pseudo = n_pseudo, B = B, K = K, alpha = alpha,
               beta = beta, epsilon = epsilon, with_P = with_P)
    return(out)
}







# >> Lower alates with Predators ----

#' This somehow results in more alates with Pseudomonas:
#' list(B = 0.1, K = 12500 * 1, alpha = 0,
#' beta = -2, epsilon = 1, with_P = TRUE,
#' max_t = 100, n_x = 3, n_y = 3, radius = 1.5, n_sims = 1e3,
#' time_inf_np = np, demog_error = FALSE, disaster_p = 0)


{
    np <- 5L
    args <- list(max_t = 100, n_x = 3, n_y = 3, radius = 1, n_sims = 1e3,
                 epsilon = 1,
                 alpha = 0,
                 beta = -2,
                 B = 0.1,
                 K = 12500 * 1,
                 time_inf_np = np,
                 with_P = TRUE,
                 # doomsday = TRUE,
                 # p_N0 = 0.5,
                 # wasp_disp_m0 = 0,
                 wasp_disp_m1 = 0.349 * 0.1,
                 demog_error = FALSE,
                 disaster_p = 0)
    d <- bind_rows(do.call(one_sim_combo, c(list(n_pseudo = 3L), args)),
                   do.call(one_sim_combo, c(list(n_pseudo = 0L), args))) |>
        # mutate(infect_time = ifelse(is.infinite(infect_time), args$max_t * 1.5,
        #                             infect_time)) |>
        mutate(n_pseudo = factor(n_pseudo)) |>
        select(rep:outbreak_size, n_pseudo)
    # if (any(is.infinite(d$infect_time))) warning("Some times are infinite.")
    p1 <- d |>
        ggplot(aes(n_pseudo, log_alates, color = n_pseudo)) +
        geom_violin(aes(fill = n_pseudo), alpha = 0.25) +
        stat_summary(fun = mean, geom = "point") +
        scale_color_manual(values = pseudo_pal,
                           guide = "none", aesthetics = c("color", "fill")) +
        # scale_y_continuous(breaks = log10(2^(c(1,3,5)) + 1),
        #                    labels = 2^(c(1,3,5))) +
        labs(x = "Number of *Pseudomonas* patches",
             y = "Mean log<sub>10</sub>(total alates)") +
        theme(axis.title.x = element_markdown(),
              axis.title.y = element_markdown(),
              axis.ticks.x = element_blank()) +
        .trans_theme
    # p2 <- d |>
    #     ggplot(aes(n_pseudo, infect_time, color = n_pseudo)) +
    #     geom_violin(aes(fill = n_pseudo), alpha = 0.25) +
    #     stat_summary(fun = mean, geom = "point") +
    #     scale_color_manual(values = pseudo_pal,
    #                        guide = "none", aesthetics = c("color", "fill")) +
    #     scale_y_continuous("Time for 5 plants infected") +
    #     xlab("Number of *Pseudomonas* patches") +
    #     theme(axis.title = element_markdown(),
    #           axis.ticks.x = element_blank()) +
    #     .trans_theme
    p2 <- d |>
        ggplot(aes(n_pseudo, outbreak_size, color = n_pseudo)) +
        geom_violin(aes(fill = n_pseudo), alpha = 0.25) +
        stat_summary(fun = mean, geom = "point") +
        scale_color_manual(values = pseudo_pal,
                           guide = "none", aesthetics = c("color", "fill")) +
        # scale_y_continuous(breaks = c(1,5,9)) +
        labs(x = "Number of *Pseudomonas* patches",
             y = "Outbreak size") +
        theme(axis.title.x = element_markdown(),
              axis.title.y = element_markdown(),
              axis.ticks.x = element_blank()) +
        .trans_theme
    p1 + p2 + plot_layout(nrow = 1)
    # rm(p1, p2, d, args, np)

}

ggsave("~/Desktop/lower_violin_log_alates_p.svg", p1 + .no_axes,
       width = 3, height = 2, bg = "transparent")
ggsave("~/Desktop/lower_violin_outbreak_size_p.svg", p2 + .no_axes,
       width = 3, height = 2, bg = "transparent")




# Time series for above:

dts <- bind_rows(do.call(one_sim_combo, c(list(n_pseudo = 3L, summarize = FALSE,
                                               out_pseudo = TRUE),
                                          args |> modify_at("n_sims", \(x) 100))),
                 do.call(one_sim_combo, c(list(n_pseudo = 0L, summarize = FALSE,
                                               out_pseudo = TRUE),
                                          args |> modify_at("n_sims", \(x) 100)))) |>
    select(rep:n_pseudo) |>
    mutate(pseudo = factor(pseudo, labels = c("noPseudo", "Pseudo"))) |>
    mutate(n_pseudo = factor(n_pseudo))




dts |>
    mutate(plant = interaction(x, y, pseudo),
           id = interaction(rep, n_pseudo, plant)) |>
    select(-x, -y, -virus) |>
    pivot_longer(aphids:enemies, names_to = "type", values_to = "density") |>
    ggplot(aes(time, density, color = n_pseudo)) +
    geom_line(aes(group = id), alpha = 0.1) +
    # stat_summary(geom = "line", fun = mean, linewidth = 1) +
    facet_grid(type ~ n_pseudo, scales = "free_y") +
    scale_color_manual(values = pseudo_pal,
                       guide = "none", aesthetics = c("color", "fill"))

p <- levels(dts$n_pseudo) |>
    set_names() |>
    map(\(.n_pseudo){
        dts |>
            filter(n_pseudo == .n_pseudo) |>
            mutate(aphids = aphids + alates) |>
            select(-alates) |>
            pivot_longer(aphids:enemies, names_to = "type",
                         values_to = "density") |>
            mutate(plant = interaction(x, y, pseudo),
                   id = interaction(rep, n_pseudo, x, y, type)) |>
            select(-x, -y, -virus) |>
            ggplot(aes(time, density, color = type)) +
            geom_line(aes(group = id), alpha = 0.1) +
            scale_color_manual(values = spp_pal) +
            guides(color = guide_legend(override.aes = list(alpha = 1,
                                                            linewidth = 1.5))) +
            scale_y_continuous(NULL, limits = c(0, 2000), breaks = 0:2*1000) +
            .trans_theme
    })

do.call(wrap_plots, p) +
    plot_layout(nrow = 1, guides = "collect")

# ggsave("~/Desktop/lower_ts_0.svg", p[[1]] + .no_axes, width = 2, height = 1.5)
# ggsave("~/Desktop/lower_ts_3.svg", p[[2]] + .no_axes, width = 2, height = 1.5)




# dts |>
#     filter(rep == sample.int(100, 1)) |>
#     filter(virus == 1) |>
#     distinct(rep, n_pseudo, x, y)


ns_spat_p <- levels(dts$n_pseudo) |>
    set_names() |>
    map(\(.n_pseudo){
        dd <- dts |>
            filter(rep == 27, n_pseudo == .n_pseudo) |>
            mutate(plant = interaction(x, y, pseudo, lex.order = TRUE)) |>
            mutate(aphids = aphids + alates) |>
            select(-alates, -x, -y)
        ddv <- dd |>
            group_by(plant) |>
            filter(virus == 1) |>
            filter(time == min(time)) |>
            ungroup()
        y_max <- dts |>
            filter(rep == dd$rep[[1]]) |>
            mutate(aphids = aphids + alates) |>
            getElement("aphids") |>
            max()
        yax_max <- (y_max %/% 500) * 500
        dd |>
            select(-virus) |>
            pivot_longer(aphids:enemies, names_to = "type",
                         values_to = "density") |>
            ggplot(aes(time, density, color = type)) +
            geom_line(linewidth = 1) +
            geom_vline(data = ddv, aes(xintercept = time), color = "#EC008C",
                       linetype = "22") +
            geom_point(data = ddv, aes(y = y_max), shape = 8, color = "#EC008C") +
            ggtitle(sprintf("%s *Pseudomonas* plants", .n_pseudo)) +
            # stat_summary(geom = "line", fun = mean, linewidth = 1) +
            facet_wrap( ~ plant, nrow = 3) +
            scale_y_continuous(limits = c(0, y_max),
                               breaks = c(0, yax_max/2, yax_max)) +
            scale_color_manual(values = spp_pal) +
            theme(plot.title = element_markdown()) +
            .trans_theme
    })





do.call(wrap_plots, ns_spat_p) +
    plot_layout(nrow = 1, guides = "collect")

ggsave("~/Desktop/lower_ts_spat_0.svg", ns_spat_p[[1]] + .no_axes,
       width = 3, height = 3, bg = "transparent")
ggsave("~/Desktop/lower_ts_spat_3.svg", ns_spat_p[[2]] + .no_axes,
       width = 3, height = 3, bg = "transparent")









levels(dts$n_pseudo) |>
    set_names() |>
    map(\(.n_pseudo){
        dts |>
            filter(rep == 1, n_pseudo == .n_pseudo) |>
            mutate(aphids = aphids + alates) |>
            rename(parasitoid = enemies) |>
            group_by(time) |>
            summarize(parasitoid = sum(parasitoid),
                      aphids = sum(aphids)) |>
            pivot_longer(aphids:parasitoid, names_to = "type",
                         values_to = "density") |>
            ggplot(aes(time, density, color = type)) +
            geom_line(linewidth = 1) +
            ggtitle(sprintf("%s *Pseudomonas* plants", .n_pseudo)) +
            # coord_cartesian(ylim = c(0, 2e3)) +
            scale_color_viridis_d(end = 0.8, option = "plasma") +
            theme(plot.title = element_markdown())
    }) |>
    do.call(what = wrap_plots) +
    plot_layout(nrow = 1, guides = "collect")



dts |>
    group_by(n_pseudo, rep) |>
    summarize(enemies = mean(log(enemies)), .groups = "drop") |>
    ggplot(aes(n_pseudo, enemies, fill = n_pseudo, color = n_pseudo)) +
    geom_violin(alpha = 0.25) +
    stat_summary(fun = mean, geom = "point") +
    scale_color_manual(values = pseudo_pal,
                       guide = "none", aesthetics = c("color", "fill"))





# =============================================================================*
# Small simulations ----
# =============================================================================*


# Takes ~1 min
ps_sims <- crossing(n_pseudo = c(0L, 3L),
                    B = c(0.05, 0.01, 0),
                    K = 12500 * c(0.25, 0.5, 1),
                    alpha = c(0, 1, 2),
                    beta = -1 * c(0, 1, 2),
                    epsilon = 1,
                    with_P = FALSE) |>
    mutate(n_x = 3L, n_y = 3L, radius = 1, n_sims = 1000, max_t = 100,
           time_inf_np = list(2:6)) |>
    pmap(one_sim_combo, .progress = .prog_args) |>
    list_rbind() |>
    select(-epsilon, -with_P) |>
    select(n_pseudo, B, K, alpha, beta, rep, everything()) |>
    mutate(across(n_pseudo:rep, factor))


# Parameter names that differ in ps_sims:
par_names <- ps_sims |>
    select(n_pseudo:rep) |>
    select(-rep) |>
    map_int(\(x, i) length(unique(x))) |>
    discard(\(x) x <= 1) |>
    names()





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





ps_sims_plotter <- function(yvar, .np = 4L, .K_lvl = 3L) {

    # yvar = "infect_time"; .np = 4L

    stopifnot(length(yvar) == 1L && is.character(yvar))
    stopifnot(length(.np) == 1L && is.numeric(.np) && round(.np) == .np)

    yvar <- match.arg(yvar, c("infect_time", "log_alates", "alates",
                              "p_alates", "log_aphids", "aphids"))

    dd <- ps_sims |>
        filter(K == levels(K)[.K_lvl]) |>
        pretty_facet_factors(c("alpha", "B", "K"),
                             .greek = c(TRUE, FALSE, FALSE)) |>
        mutate(beta = fct_rev(beta))

    if (yvar == "infect_time") {

        .trans <- log10
        .ylab <- sprintf("log<sub>10</sub>(Days to %i plants infected)", .np)
        dd <- dd |>
            unnest(infect_time) |>
            filter(np == .np)

    } else {

        .trans <- identity
        .ylab <- list(log_alates = "Mean log<sub>10</sub>(Total alates + 1)",
                       alates = "Mean total alates",
                       p_alates = "Mean proportion alates",
                       log_aphids = "Mean log<sub>10</sub>(Total aphids + 1)",
                       aphids = "Mean total aphids")[[yvar]]

    }

    dd |>
        ggplot(aes(beta, .trans(.data[[yvar]]), color = n_pseudo)) +
        # geom_hline(aes(yintercept = min(.trans(.data[[yvar]]))),
        #            linewidth = 0.75, color = "black") +
        geom_violin(position = position_dodge(0.5), fill = NA) +
        stat_summary(fun = mean, geom = "point", position = position_dodge(0.5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",
                     position = position_dodge(0.5), width = 0.3) +
        scale_color_manual(values = pseudo_pal,
                           guide = "none", aesthetics = c("color", "fill")) +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        ylab(.ylab) +
        xlab("Effect of *Pseudomonas* on alates alighting (&beta;)") +
        facet_wrap(B ~ alpha, nrow = 3, scales = "free_y") +
        # facet_grid(K ~ beta) +
        theme(strip.text = element_markdown(size = 8),
              strip.text.y = element_markdown(angle = 0),
              axis.title = element_markdown(),
              legend.title = element_markdown(),
              panel.spacing = unit(1, "lines"))

}



ps_sims_plotter("infect_time") +
    ps_sims_plotter("log_alates") +
    plot_layout(ncol = 1, guides = "collect")

ps_sims_plotter("infect_time", .K_lvl = 2L) +
    ps_sims_plotter("log_alates", .K_lvl = 2L) +
    plot_layout(ncol = 1, guides = "collect")

# ps_sims_plotter("log_alates", .K_lvl = 1L) +
ps_sims_plotter("log_alates", .K_lvl = 2L) +
    ps_sims_plotter("log_alates", .K_lvl = 3L) +
    plot_layout(ncol = 1, guides = "collect")

# ps_sims_plotter("log_aphids")



# >> First results talk plot ----

# yvar <- "log_alates"
yvar <- "infect_time"
# yvar <- "outbreak_size"

p <- ps_sims |>
    unnest(infect_time) |>
    filter(np == 5) |>
    filter(K == levels(K)[3], alpha == 1, beta == -2, B == 0.05) |>
    ggplot(aes(n_pseudo, .data[[yvar]], color = n_pseudo)) +
    geom_violin(aes(fill = n_pseudo), alpha = 0.25) +
    stat_summary(fun = mean, geom = "point") +
    # scale_color_viridis_d(option = "viridis", end = 0.8,
    scale_color_manual(values = pseudo_pal,
                       guide = "none", aesthetics = c("color", "fill")) +
    scale_y_continuous(list(infect_time = "Days to 5 plants infected",
                            alates = "Mean total alates",
                            log_alates = "Mean log<sub>10</sub>(total alates)",
                            outbreak_size = "Outbreak size")[[yvar]],
                       breaks = if (str_starts(yvar, "log")) {
                           log10(c(800, 900, 1000)) } else waiver(),
                       labels = if (str_starts(yvar, "log")) {
                           c(800, 900, 1000) } else waiver()) +
    xlab("Number of *Pseudomonas* patches") +
    # facet_grid(n_pseudo ~ ., scales = "free_y") +
    # facet_grid(K ~ beta) +
    theme(axis.ticks.x = element_blank(),
          axis.title = element_markdown()) +
    .trans_theme


p

ggsave(sprintf("~/Desktop/violin_%s_p.svg", yvar), p + .no_axes,
       width = 3, height = 2, bg = "transparent")



ps_sims |>
    filter(n_pseudo %in% levels(n_pseudo)[c(1,3)],
           epsilon == levels(epsilon)[2],
           K == levels(K)[3]) |>
           # alpha == levels(alpha)[1]) |>
    pretty_facet_factors(c("alpha", "B", "K"),
                         .greek = c(rep(TRUE, 2), FALSE)) |>
    unnest(infect_time) |>
    filter(np == .np) |>
    ggplot(aes(beta, log10(infect_time), color = n_pseudo)) +
    geom_hline(aes(yintercept = min(log10(.data[["infect_time"]]))),
               linewidth = 0.75, color = "black") +
    geom_violin(position = position_dodge(0.5), fill = NA) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0.5)) +
    scale_color_viridis_d("Number of<br>*Pseudomonas*<br>patches",
                          option = "viridis", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ylab(sprintf("log<sub>10</sub>(Days to %i plants infected)", .np)) +
    xlab("*Pseudomonas*-induced mortality") +
    facet_grid(alpha ~ B) +
    # facet_grid(K ~ beta) +
    theme(strip.text = element_markdown(size = 8),
          axis.title = element_markdown(),
          legend.title = element_markdown())




ps_sims |>
    filter(n_pseudo %in% levels(n_pseudo)[c(1,3)],
           with_P == FALSE,
           K == levels(K)[3]) |>
           # alpha == levels(alpha)[2]) |>
    pretty_facet_factors(c("alpha", "beta", "K"),
                         .greek = c(rep(TRUE, 2), FALSE)) |>
    ggplot(aes(B, log_alates, color = n_pseudo)) +
    # geom_hline(aes(yintercept = 0), #min(.data[["alates"]])),
    #            linewidth = 0.75, color = "black") +
    geom_violin(position = position_dodge(0.5), fill = NA) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0.5)) +
    scale_color_viridis_d("Number of<br>*Pseudomonas*<br>patches",
                          option = "viridis", end = 0.8) +
    ylab("Mean log<sub>10</sub>(Total alates)") +
    xlab("Effect of *Pseudomonas* on alates alighting (&beta;)") +
    facet_grid(alpha ~ beta) +
    # facet_grid(K ~ beta) +
    theme(strip.text = element_markdown(size = 8),
          axis.title = element_markdown(),
          legend.title = element_markdown())



















# ---------------------------*
# How often does Pseudomonas increase time to fully infected? ----
# ---------------------------*


ps_inc_sims <- ps_sims |>
    split(as.formula(paste("~", paste(par_names[-1], collapse = "+")))) |>
    future_lapply(\(d) {
        .a0 <- d$alates[d$n_pseudo == 0]
        .pa0 <- d$p_alates[d$n_pseudo == 0]
        .it20 <- sapply(d$infect_time[d$n_pseudo == 0], \(x) x[[2]][[1]])
        .it30 <- sapply(d$infect_time[d$n_pseudo == 0], \(x) x[[2]][[2]])
        .it40 <- sapply(d$infect_time[d$n_pseudo == 0], \(x) x[[2]][[3]])

        d1 <- d[d$n_pseudo != 0,]

        d1[["p_more"]] <- sapply(d1$alates, \(a) mean(a > .a0))
        d1[["p_more_p"]] <- sapply(d1$p_alates, \(pa) mean(pa > .pa0))

        it2 <- sapply(d1$infect_time, \(x) x[[2]][[1]])
        it3 <- sapply(d1$infect_time, \(x) x[[2]][[2]])
        it4 <- sapply(d1$infect_time, \(x) x[[2]][[3]])

        d1[["p_longer2"]] <- sapply(it2, \(it) mean(it > .it20))
        d1[["p_longer3"]] <- sapply(it3, \(it) mean(it > .it30))
        d1[["p_longer4"]] <- sapply(it4, \(it) mean(it > .it40))
        return(d1)
    }) |>
    list_rbind() |>
    mutate(n_pseudo = fct_drop(n_pseudo))

pc_inc_plotter <- function(yvar) {
    if (str_starts(yvar, "p_longer")) {
        .ylab <- sprintf("Percent longer to infect %s plants",
                         str_remove(yvar, "p_longer"))
    } else {
        .ylab <- list(p_more = "Percent more alates",
                      p_more_p = "Percent more alate proportion")[[yvar]]
    }
    ps_inc_sims |>
        filter(# beta %in% levels(beta)[c(1,4)],
            alpha %in% levels(alpha)[c(1,3)]) |>
        # filter(K %in% levels(K)[c(1,3)]) |>
        # alpha == levels(alpha)[1]) |>
        pretty_facet_factors(c("alpha", "K"), .greek = c(TRUE, FALSE)) |>
        mutate(beta = fct_rev(beta)) |>
        ggplot(aes(beta, .data[[yvar]] * 100, color = B)) +
        geom_hline(aes(yintercept = 50),
                   linewidth = 0.75, color = "gray70", linetype = "22") +
        geom_violin(position = position_dodge(0.6), fill = NA) +
        stat_summary(fun = mean, geom = "point", position = position_dodge(0.6), size = 3) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",
                     position = position_dodge(0.6), width = 0.3) +
        scale_color_viridis_d("*Pseudomonas*-<br>induced<br>mortality",
                              option = "plasma", end = 0.8) +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        ylab(.ylab) +
        xlab("Effect of *Pseudomonas* on alates alighting (&beta;)") +
        # facet_grid(. ~ alpha) +
        facet_grid(K ~ alpha) +
        theme_minimal() +
        theme(strip.text = element_markdown(size = 10),
              strip.text.y = element_markdown(size = 10, angle = 0),
              axis.title = element_markdown(),
              legend.title = element_markdown(),
              legend.text = element_markdown())
}

pc_inc_plotter("p_longer4") +
    pc_inc_plotter("p_more") +
    plot_layout(ncol = 1, guides = "collect") &
    theme(panel.spacing = unit(2, "lines"))



ps_inc_sims |>
    filter(beta %in% levels(beta)[c(1,4)],
           alpha %in% levels(alpha)[c(1,4)]) |>
    filter(n_pseudo == levels(n_pseudo)[1],
           K == levels(K)[2]) |>
    # alpha == levels(alpha)[1]) |>
    # pretty_facet_factors(c("alpha", "beta", "epsilon", "K"),
                         # .greek = c(rep(TRUE, 3), FALSE)) |>
    pretty_facet_factors(c("alpha", "beta"), .greek = TRUE) |>
    (\(.data) {
        .title <<- sprintf("%i *Pseudomonas* patches<br>K = %s",
                           as.integer(paste(.data[["n_pseudo"]][[1]])),
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


