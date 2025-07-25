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


#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists("LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}


# Carrying capacity with no alates or natural enemies:
CC <- function(K, B){
    L  <- rbind(c(pop_info$surv_j, pop_info$fecund),
                c(pop_info$recruit, pop_info$surv_a))
    K * (max(abs(eigen(L)[["values"]])) * (1 - B) - 1)
}




insect_sims <- map(1:10, \(i){
    test_insect_pops(max_t = 100,
                     B = 0,
                     h = 5,
                     alate_0 = -Inf,
                     alate_1 = 0,
                     demog_error = TRUE,
                     disaster_p = 0.02,
                     disaster_s = 0.1,
                     A0 = 10,
                     W0 = 0,
                     P0 = 0) |>
        mutate(aphids = aphids + alates) |>
        select(-alates) |>
        pivot_longer(aphids:enemies, names_to = "species", values_to = "N") |>
        mutate(species = factor(species, levels = c("aphids", "enemies")),
               rep = i)
}) |>
    list_rbind() |>
    mutate(rep = factor(rep))
insect_sims2 <- test_insect_pops(max_t = 100,
                                B = 0,
                                h = 5,
                                alate_0 = -Inf,
                                alate_1 = 0,
                                demog_error = FALSE,
                                A0 = 10,
                                W0 = 0,
                                P0 = 0) |>
    mutate(aphids = aphids + alates) |>
    select(-alates) |>
    pivot_longer(aphids:enemies, names_to = "species", values_to = "N") |>
    mutate(species = factor(species, levels = c("aphids", "enemies")))

# insect_time_pts <- c(375, 434, 500, 577)
# insect_sims |>
#     filter(time %in% insect_time_pts)

insect_sims |>
    # filter(species != "alates") |>
    ggplot(aes(time, N, color = species)) +
    geom_hline(yintercept = CC(12500, 0), color = "gray80", linewidth = 1) +
    geom_line(aes(group = interaction(rep, species)), linewidth = 1, alpha = 0.25) +
    geom_line(data = insect_sims2 |> filter(species == "aphids"), linewidth = 1,
              color = "red") +
    # geom_vline(xintercept = insect_time_pts, linetype = "22",
    #            color = "gray70", linewidth = 1) +
    scale_color_viridis_d(begin = 0.2, end = 0.9) +
    scale_y_continuous("Abundance") +
    theme_minimal()

# starts <- list(A0 = insect_sims$aphids[insect_sims$time %in% insect_time_pts],
#                P0 = insect_sims$enemies[insect_sims$time %in% insect_time_pts])



make_arg_list <- function(pseudo, B, K, alpha, beta, epsilon, with_P,
                          max_t, n_x, n_y, radius, n_sims, ...) {


    stopifnot(round(n_x) == n_x && round(n_y) == n_y && round(n_sims) == n_sims)
    n_x <- as.integer(n_x)
    n_y <- as.integer(n_y)
    n_plants <- n_x * n_y
    stopifnot(round(pseudo) == pseudo && pseudo >= 0 && pseudo <= (n_plants-1L))

    land <- array(0L, c(n_x, n_y, n_sims))
    A0 <- array(0.0, c(n_x, n_y, n_sims))
    P0 <- array(0.0, c(n_x, n_y, n_sims))
    for (i in 1:n_sims) {
        land[1,1,i] <- 1L
        if (pseudo > 0) {
            k <- sample.int(n_plants - 1L, pseudo)
            x <- k - n_x * (k %/% n_x) + 1L
            y <- k %/% n_x + 1L
            land[cbind(x,y,i)] <- 2L
        }

        # if (with_P) {
        #     idx <- sample.int(length(starts$A0), n_plants, replace = TRUE)
        #     A0[,,i] <- starts$A0[idx]
        #     P0[,,i] <- starts$P0[idx]
        # } else A0[,,i] <- sample(starts$A0, n_plants, replace = TRUE)\

        # A0[,,i] <- 10
        A0[,,i] <- exp(runif(n_plants, -3, 5))
        # A0[,,i] <- runif(n_plants, 0, 100)
        if (with_P) {
            # P0[,,i] <- runif(n_plants, 0, 3)
            P0[,,i] <- A0[,,i] * exp(runif(length(A0[,,i]), -10, -1))
        }

    }

    insect_args <- list(K = K, B = B, h = 5, fly_p = 0.05, wasp_d_p = 0.1)
    plant_args <- list(landscapes = land,
                       max_t = max_t,
                       A0 = A0,
                       W0 = array(0.0, c(n_x, n_y, n_sims)),
                       P0 = P0,
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
            if (n %in% c("A0", "W0", "P0", "landscape")) {
                plant_args[[n]] <- array(other_args[[n]], c(n_x, n_y, n_sims))
            } else plant_args[[n]] <- other_args[[n]]
        }
    }


    plant_args[["insect_ptr"]] <- do.call(make_insect_ptr, insect_args)

    return(plant_args)

}





one_sim_combo <- function(pseudo, B, K, alpha, beta, epsilon, with_P,
                          max_t, n_x, n_y, radius, n_sims, ...,
                          time_inf_np = NULL) {

    n_plants <- n_x * n_y
    # .A0 <- rep(10, n_plants * n_sims) # runif(n_plants * n_sims, 0, 1000)
    .A0 <- 10  # exp(runif(n_plants * n_sims, -3, 5))
    .P0 <- 0
    if (with_P) .P0 <- .A0 * exp(runif(length(.A0), -10, -1))

    arg_list <- make_arg_list(pseudo, B, K, alpha, beta, epsilon, with_P,
                              max_t, n_x, n_y, radius, n_sims,
                              A0 = .A0, P0 = .P0, ...)


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
        mutate(aphids = alates + aphids) |>
        group_by(rep) |>
        summarize(p_alates = mean(alates[aphids > 0] / aphids[aphids > 0]),
                  aphids = mean(aphids),
                  log_aphids = mean(log10(aphids)),
                  alates = mean(alates),
                  log_alates = mean(log10(alates)),
                  infect_time = time_inf_fun(time_inf_np, time, virus),
                  outbreak_size = max(virus)) |>
        mutate(pseudo = pseudo, B = B, K = K, alpha = alpha,
               beta = beta, epsilon = epsilon, with_P = with_P)
    return(out)
}







# LEFT OFF ----
#'
#' Does Pseudomonas do ANYTHING to total alates??
#'

{
    np <- 4L
    args <- list(B = 0, K = 12500 * 1, alpha = 0,
                 beta = -10, epsilon = 1, with_P = TRUE,
                 max_t = 50, n_x = 3, n_y = 3, radius = 1, n_sims = 10e3,
                 time_inf_np = np, demog_error = TRUE, disaster_p = 0.02)
    d <- bind_rows(do.call(one_sim_combo, c(list(pseudo = 3L), args)),
                   do.call(one_sim_combo, c(list(pseudo = 0L), args))) |>
        # mutate(infect_time = ifelse(is.infinite(infect_time), 200, infect_time)) |>
        mutate(pseudo = factor(pseudo))
    if (any(is.infinite(d$infect_time))) warning("Some times are infinite.")
    p1 <- d |>
        ggplot(aes(pseudo, log10(infect_time), color = pseudo)) +
        # geom_jitter(shape = 1, alpha = 0.1) +
        geom_violin() +
        stat_summary(fun = mean, geom = "point") +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
        scale_color_viridis_d(guide = "none") +
        ylab(sprintf("log<sub>10</sub>(Days to %i plants infected)", np)) +
        xlab("Number of *Pseudomonas* patches") +
        theme(axis.title = element_markdown())
    p2 <- d |>
        ggplot(aes(pseudo, log_alates, color = pseudo)) +
        # geom_jitter(shape = 1, alpha = 0.1) +
        geom_violin() +
        stat_summary(fun = mean, geom = "point") +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
        scale_color_viridis_d(guide = "none") +
        ylab("Mean log<sub>10</sub>(Total alates)") +
        xlab("Number of *Pseudomonas* patches") +
        theme(axis.title = element_markdown())
    p1 + p2 + plot_layout(nrow = 1)
    # rm(p1, p2, d, args, np)
}






# =============================================================================*
# Small simulations ----
# =============================================================================*


if (!file.exists("_building/ps_sims.rds")) {
    # Takes ~15 sec for 3x3 landscape
    set.seed(1114260777)
    ps_sims <- crossing(pseudo = c(0L, 1L, 3L),
                        B = c(0.05, 0.01, 0),
                        K = 12500 * c(0.25, 0.5, 1),
                        alpha = c(0, 1, 2),
                        beta = -1 * c(0, 1, 2),
                        epsilon = 1,
                        with_P = FALSE) |>
        mutate(n_x = 3L, n_y = 3L, radius = 1, n_sims = 1000, max_t = 50,
               time_inf_np = list(2:6)) |>
        pmap(one_sim_combo, .progress = .prog_args) |>
        list_rbind() |>
        select(-epsilon, -with_P) |>
        select(pseudo, B, K, alpha, beta, rep, everything()) |>
        mutate(across(pseudo:rep, factor))
    write_rds(ps_sims, "_building/ps_sims.rds", compress = "xz")
} else {
    ps_sims <- read_rds("_building/ps_sims.rds")
}




# Parameter names that differ in ps_sims:
par_names <- ps_sims |>
    select(pseudo:rep) |>
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
        filter(pseudo %in% levels(pseudo)[c(1,3)],
               K == levels(K)[.K_lvl]) |>
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
        .ylab <- list(log_alates = "Mean log<sub>10</sub>(Total alates)",
                       alates = "Mean total alates",
                       p_alates = "Mean proportion alates",
                       log_aphids = "Mean log<sub>10</sub>(Total aphids)",
                       aphids = "Mean total aphids")[[yvar]]

    }

    dd |>
        ggplot(aes(beta, .trans(.data[[yvar]]), color = pseudo)) +
        # geom_hline(aes(yintercept = min(.trans(.data[[yvar]]))),
        #            linewidth = 0.75, color = "black") +
        # geom_violin(position = position_dodge(0.5), fill = NA) +
        stat_summary(fun = mean, geom = "point", position = position_dodge(0.5)) +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar",
                     position = position_dodge(0.5), width = 0.3) +
        scale_color_viridis_d("Number of<br>*Pseudomonas*<br>patches",
                              option = "viridis", end = 0.8) +
        guides(color = guide_legend(override.aes = list(alpha = 1))) +
        ylab(.ylab) +
        xlab("Effect of *Pseudomonas* on alates alighting (&beta;)") +
        facet_grid(B ~ alpha, scales = "free_y") +
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
ps_sims_plotter("alates", .K_lvl = 2L) +
    ps_sims_plotter("alates", .K_lvl = 3L) +
    plot_layout(ncol = 1, guides = "collect")

# ps_sims_plotter("log_aphids")









ps_sims |>
    filter(pseudo %in% levels(pseudo)[c(1,3)],
           epsilon == levels(epsilon)[2],
           K == levels(K)[3]) |>
           # alpha == levels(alpha)[1]) |>
    pretty_facet_factors(c("alpha", "B", "K"),
                         .greek = c(rep(TRUE, 2), FALSE)) |>
    unnest(infect_time) |>
    filter(np == .np) |>
    ggplot(aes(beta, log10(infect_time), color = pseudo)) +
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
    filter(pseudo %in% levels(pseudo)[c(1,3)],
           with_P == FALSE,
           K == levels(K)[3]) |>
           # alpha == levels(alpha)[2]) |>
    pretty_facet_factors(c("alpha", "beta", "K"),
                         .greek = c(rep(TRUE, 2), FALSE)) |>
    ggplot(aes(B, log_alates, color = pseudo)) +
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
    filter(pseudo %in% levels(pseudo)[c(1,3)]) |>
    split(as.formula(paste("~", paste(par_names[-1], collapse = "+")))) |>
    future_lapply(\(d) {
        .a0 <- d$alates[d$pseudo == 0]
        .pa0 <- d$p_alates[d$pseudo == 0]
        .it20 <- sapply(d$infect_time[d$pseudo == 0], \(x) x[[2]][[1]])
        .it30 <- sapply(d$infect_time[d$pseudo == 0], \(x) x[[2]][[2]])
        .it40 <- sapply(d$infect_time[d$pseudo == 0], \(x) x[[2]][[3]])

        d1 <- d[d$pseudo != 0,]

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
    mutate(pseudo = fct_drop(pseudo))

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
        filter(epsilon == 1) |>
        # filter(K %in% levels(K)[c(1,3)]) |>
        # alpha == levels(alpha)[1]) |>
        pretty_facet_factors(c("alpha", "epsilon", "K"),
                             .greek = c(rep(TRUE, 2), FALSE)) |>
        mutate(beta = fct_rev(beta)) |>
        # mutate(with_P = factor(with_P, levels = c("without parasitoids", "with parasitoids"),
        #                        labels = c("absent", "present"))) |>
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
        pmap(one_sim_combo, .progress = .prog_args) |>
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






