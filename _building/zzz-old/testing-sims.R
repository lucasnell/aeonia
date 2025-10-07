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


spp_pal <- viridisLite::plasma(101)[c(10, 50, 90)] |>
    set_names(c("aphids", "alates", "enemies"))


p <- test_insect_pops(max_t = 100,
                 B = 0,
                 demog_error = FALSE,
                 A0 = 10,
                 W0 = 0,
                 P0 = 0) |>
    select(-enemies) |>
    pivot_longer(aphids:alates, names_to = "morph", values_to = "N") |>
    mutate(morph = factor(morph, levels = c("aphids", "alates"))) |>
    # filter(morph == "aphids") |>
    ggplot(aes(time, N, color = morph)) +
    # geom_hline(yintercept = CC(12500, 0, 0.1), color = "gray80", linewidth = 1) +
    geom_line(linewidth = 1) +
    # scale_color_viridis_d(begin = 0.2, end = 0.9, drop = FALSE) +
    scale_color_manual(values = spp_pal, guide = "none") +
    ylab("Abundance") +
    coord_cartesian(ylim = c(0, 1600)) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          panel.background = element_rect(fill="transparent"), #transparent panel bg
        plot.background = element_rect(fill="transparent", color=NA), #transparent plot bg
        panel.grid.major = element_blank(), #remove major gridlines
        panel.grid.minor = element_blank(), #remove minor gridlines
        legend.background = element_rect(fill="transparent"), #transparent legend bg
        legend.box.background = element_rect(fill="transparent") #transparent legend panel
    )

ggsave("~/Desktop/p.svg", p, width = 3, height = 2, bg = "transparent")





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
    geom_hline(yintercept = CC(12500, 0, 0.1), color = "gray80", linewidth = 1) +
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



make_arg_list <- function(n_pseudo, B, K, alpha, beta, epsilon, with_P,
                          max_t, n_x, n_y, radius, n_sims, ...) {


    stopifnot(round(n_x) == n_x && round(n_y) == n_y && round(n_sims) == n_sims)
    n_x <- as.integer(n_x)
    n_y <- as.integer(n_y)
    n_plants <- n_x * n_y
    stopifnot(round(n_pseudo) == n_pseudo && n_pseudo >= 0 && n_pseudo <= (n_plants-1L))

    land <- array(0L, c(n_x, n_y, n_sims))
    A0 <- array(0.0, c(n_x, n_y, n_sims))
    P0 <- array(0.0, c(n_x, n_y, n_sims))
    for (i in 1:n_sims) {
        land[1,1,i] <- 1L
        if (n_pseudo > 0) {
            k <- sample.int(n_plants - 1L, n_pseudo)
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

    insect_args <- list(K = K, B = B, h = 5, fly_p = 0.05, wasp_disp_m0 = 0.3)
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





one_sim_combo <- function(n_pseudo, B, K, alpha, beta, epsilon, with_P,
                          max_t, n_x, n_y, radius, n_sims, ...,
                          summarize = TRUE, doomsday = FALSE, p_A0 = 1,
                          time_inf_np = NULL) {

    n_plants <- n_x * n_y
    if (p_A0 < 1 && doomsday) stop("Ambiguous: p_A0 < 1 && doomsday")
    if (doomsday) {
        .A0 <- c(10, rep(0, n_plants-1L))
    } else {
        # .A0 <- rep(10, n_plants * n_sims) # runif(n_plants * n_sims, 0, 1000)
        .A0 <- rep(10, n_plants * n_sims)  # exp(runif(n_plants * n_sims, -3, 5))
        if (p_A0 < 1) {
            for (i in 1:n_sims) {
                idx0 <- (i - 1L) * n_plants
                idx <- sample.int(n_plants, round(p_A0 * n_plants))
                .A0[idx0 + idx] <- 0
            }
        }
    }
    .P0 <- 0
    if (with_P) .P0 <- .A0 * exp(runif(n_plants * n_sims, -10, -1))

    arg_list <- make_arg_list(n_pseudo, B, K, alpha, beta, epsilon, with_P,
                              max_t, n_x, n_y, radius, n_sims,
                              A0 = .A0, P0 = .P0, ...)

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

    out <- do.call(sim_plantscape, arg_list)


    if (summarize) {
        out <- out |>
            mutate(aphids = alates + aphids) |>
            group_by(rep) |>
            summarize(p_alates = mean(alates[aphids > 0] / aphids[aphids > 0]),
                      log_aphids = mean(log10(aphids)),
                      aphids = mean(aphids),
                      log_alates = mean(log10(alates)),
                      alates = mean(alates),
                      infect_time = time_inf_fun(time_inf_np, time, virus),
                      outbreak_size = max(virus))
    }
    out <- out |>
        mutate(n_pseudo = n_pseudo, B = B, K = K, alpha = alpha,
               beta = beta, epsilon = epsilon, with_P = with_P)
    return(out)
}







# LEFT OFF ----
#'
#' Does Pseudomonas do ANYTHING to total alates??
#'


#' This somehow results in more alates with Pseudomonas:
#' list(B = 0.1, K = 12500 * 1, alpha = 0,
#' beta = -2, epsilon = 1, with_P = TRUE,
#' max_t = 100, n_x = 3, n_y = 3, radius = 1.5, n_sims = 1e3,
#' time_inf_np = np, demog_error = FALSE, disaster_p = 0)


{
    np <- 2L
    args <- list(max_t = 100, n_x = 3, n_y = 3, radius = 1, n_sims = 1e3,
                 epsilon = 1,
                 alpha = 0,
                 beta = -2,
                 B = 0.1,
                 K = 12500 * 1,
                 time_inf_np = np,
                 with_P = TRUE,
                 # doomsday = TRUE,
                 # p_A0 = 0.5,
                 wasp_disp_m0 = 0.1,
                 wasp_disp_m1 = 0.349 * 0,
                 demog_error = FALSE,
                 disaster_p = 0)
    d <- bind_rows(do.call(one_sim_combo, c(list(n_pseudo = 3L), args)),
                   do.call(one_sim_combo, c(list(n_pseudo = 0L), args))) |>
        mutate(infect_time = ifelse(is.infinite(infect_time), args$max_t * 1.5,
                                    infect_time)) |>
        mutate(n_pseudo = factor(n_pseudo))
    if (any(is.infinite(d$infect_time))) warning("Some times are infinite.")
    p1 <- d |>
        ggplot(aes(n_pseudo, log10(infect_time), color = n_pseudo)) +
        # geom_jitter(shape = 1, alpha = 0.1) +
        geom_violin() +
        stat_summary(fun = mean, geom = "point") +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
        scale_color_viridis_d(end = 0.8, guide = "none") +
        ylab(sprintf("log<sub>10</sub>(Days to %i plants infected)", np)) +
        xlab("Number of *Pseudomonas* patches") +
        theme(axis.title = element_markdown())
    p2 <- d |>
        ggplot(aes(n_pseudo, log10(alates), color = n_pseudo)) +
        # geom_jitter(shape = 1, alpha = 0.1) +
        geom_violin() +
        stat_summary(fun = mean, geom = "point") +
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", width = 0.1) +
        scale_color_viridis_d(end = 0.8, guide = "none") +
        ylab("Mean log<sub>10</sub>(Total alates)") +
        xlab("Number of *Pseudomonas* patches") +
        theme(axis.title = element_markdown())
    p1 + p2 + plot_layout(nrow = 1)
    # rm(p1, p2, d, args, np)
}


# Time series for above:

dts <- bind_rows(do.call(one_sim_combo, c(list(n_pseudo = 3L, summarize = FALSE,
                                               out_pseudo = TRUE),
                                          args |> modify_at("n_sims", \(x) 100))),
                 do.call(one_sim_combo, c(list(n_pseudo = 0L, summarize = FALSE,
                                               out_pseudo = TRUE),
                                          args |> modify_at("n_sims", \(x) 100)))) |>
    select(rep:n_pseudo) |>
    mutate(pseudo = factor(pseudo, labels = c("noPseudo", "Pseudo")))

dts |>
    mutate(n_pseudo = factor(n_pseudo),
           plant = interaction(x, y, pseudo),
           id = interaction(rep, n_pseudo, plant)) |>
    select(-x, -y, -virus) |>
    pivot_longer(aphids:enemies, names_to = "type", values_to = "density") |>
    ggplot(aes(time, density, color = n_pseudo)) +
    geom_line(aes(group = id), alpha = 0.1) +
    # stat_summary(geom = "line", fun = mean, linewidth = 1) +
    facet_grid(type ~ n_pseudo, scales = "free_y") +
    scale_color_viridis_d(end = 0.8, guide = "none")


dts |>
    filter(rep == 1, n_pseudo == 3) |>
    mutate(plant = interaction(x, y, pseudo)) |>
    # filter(plant == "3.1.Pseudo") |>
    select(-x, -y, -virus) |>
    pivot_longer(aphids:enemies, names_to = "type", values_to = "density") |>
    ggplot(aes(time, density, color = type)) +
    geom_line(linewidth = 1) +
    # stat_summary(geom = "line", fun = mean, linewidth = 1) +
    facet_wrap( ~ plant, nrow = 3) +
    scale_color_viridis_d(end = 0.8, option = "plasma")





# =============================================================================*
# Small simulations ----
# =============================================================================*


if (!file.exists("_building/ps_sims.rds")) {
    # Takes ~15 sec for 3x3 landscape
    set.seed(1114260777)
    ps_sims <- crossing(n_pseudo = c(0L, 1L, 3L),
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
        select(n_pseudo, B, K, alpha, beta, rep, everything()) |>
        mutate(across(n_pseudo:rep, factor))
    write_rds(ps_sims, "_building/ps_sims.rds", compress = "xz")
} else {
    ps_sims <- read_rds("_building/ps_sims.rds")
}




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
        filter(n_pseudo %in% levels(n_pseudo)[c(1,3)],
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
        ggplot(aes(beta, .trans(.data[[yvar]]), color = n_pseudo)) +
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
    filter(n_pseudo %in% levels(n_pseudo)[c(1,3)]) |>
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





# =============================================================================*
# =============================================================================*
# LARGER SIMULATIONS ========================
# =============================================================================*
# =============================================================================*

if (!file.exists("_building/big_ps_sims.rds")) {
    # Takes ~8 hrs for 134x134 landscape (and 12 instead of 100 sims)
    t0 <- Sys.time()
    set.seed(1952926471)
    big_ps_sims <- crossing(n_pseudo = c(0.0, 0.2, 0.5),
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
        select(n_pseudo, B, K, alpha, beta, epsilon, rep, everything()) |>
        mutate(n_pseudo = round(n_pseudo, digits = 1)) |>
        mutate(across(n_pseudo:rep, factor))
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
    ggplot(aes(B, (outbreak_size), color = n_pseudo)) +
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
    mutate(rep = factor(rep), n_pseudo = factor(n_pseudo)) |>
    group_by(n_pseudo, rep, time) |>
    summarize(virus = sum(virus), .groups = "drop") |>
    mutate(id = interaction(rep, n_pseudo, drop = TRUE)) |>
    ggplot(aes(time, virus)) +
    geom_hline(yintercept = c(0, length(land0[,,1])),
               linetype = "22", color = "gray70") +
    geom_line(aes(group = id, color = n_pseudo), alpha = 0.25) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    # facet_wrap(~ n_pseudo, ncol = 1) +
    theme_minimal()

ps_sims |>
    mutate(rep = factor(rep), n_pseudo = factor(n_pseudo)) |>
    group_by(n_pseudo, rep, time) |>
    summarize(virus = sum(virus), .groups = "drop") |>
    group_by(n_pseudo, rep) |>
    summarize(time = (\(t,v) {
        if (any(v == 9)) return(t[v == 9][1])
        return(Inf)
    })(time, virus),
              .groups = "drop") |>
    ggplot(aes(n_pseudo, time)) +
    geom_jitter(aes(color = n_pseudo), alpha = 0.25, width = 0.2, height = 0) +
    stat_summary(fun.data = "mean_cl_boot") +
    # ggplot(aes(time)) +
    # geom_freqpoly(aes(color = n_pseudo), bins = 10) +
    scale_color_viridis_d(option = "plasma", end = 0.8) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    ylab("Time to fully infected") +
    theme_minimal()

ps_sims |>
    filter(rep == 1) |>
    select(-virus) |>
    pivot_longer(aphids:preds, names_to = "type", values_to = "density") |>
    mutate(type = factor(type, levels = c("aphids", "alates", "preds")),
           id = interaction(n_pseudo, rep, type, x, y, drop = TRUE)) |>
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
           id = interaction(n_pseudo, rep, type, plant, drop = TRUE)) |>
    filter(type != "preds", type != "virus") |>
    ggplot(aes(time, density, color = type)) +
    geom_line(aes(group = id)) +
    scale_color_viridis_d(begin = 0.2, end = 0.9) +
    theme_minimal() +
    facet_grid(n_pseudo ~ plant)






