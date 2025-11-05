
#'
#' Small-scale sensitivity via Sobol indices
#'
#'


source("_scripts/01-sensitivity/00-preamble.R")
source("_scripts/01-sensitivity/sobol-preamble.R")


#' Directory containing output from `sobol-sens-paired.R`:
sobol_dir <- "~/_globus"



if (!file.exists("_scripts/interm-data/sobol-sims-summs.rds")) {

    #' Output from `sobol-sens-paired.R`:
    sobol_sims <- paste0(sobol_dir, "/sobol-sims-paired.rds") |>
        read_rds()
    # Summarize each set of simulations:
    # Takes ~1 min  (multithreading doesn't help)
    sobol_summs <- imap(sobol_sims, \(sim_set, i) {
        out_df <- sim_set |>
            mutate(combo = factor(i, levels = 1:length(sobol_sims)),
                   sims = map(sims, \(s) s[1, names(vary_pars)])) |>
            unnest(sims) |>
            select(combo, everything())
        for (yv in c(yvars, "p_outbreak", "outbreak_size2", "sd_outbreak_size")) {
            out_df[[yv]] <- 0.0
        }
        for (j in 1:nrow(sim_set)) {
            for (yv in yvars) {
                y <- sim_set$sims[[j]][[yv]]
                if (yv != "infect_time" && any(is.na(y))) stop(y, " has NA values")
                out_df[[yv]][[j]] <- mean(y, na.rm = TRUE)
            }
            y <- sim_set$sims[[j]][["outbreak_size"]]
            # now do prob. outbreak happened:
            out_df[["p_outbreak"]][[j]] <- mean(y > 1)
            # and outbreak size when there was one:
            out_df[["outbreak_size2"]][[j]] <- mean(y[y > 1])
            # Lastly, SD(outbreak size):
            out_df[["sd_outbreak_size"]][[j]] <- sd(y)
        }
        return(out_df)

    }, .progress = .prog_args) |>
        list_rbind()

    write_rds(sobol_summs, "_scripts/interm-data/sobol-sims-summs.rds",
              compress = "gz")
    rm(sobol_sims); gc()

} else {

    sobol_summs <- read_rds("_scripts/interm-data/sobol-sims-summs.rds")

}




# Same thing but looking at differences between with and without Pseudo:
diff_sobol_summs <- sobol_summs |>
    group_by(combo, alate_dens, across(all_of(names(vary_pars)))) |>
    summarize(across(all_of(c(yvars, "p_outbreak", "outbreak_size2",
                              "sd_outbreak_size")),
                     \(x) x[n_pseudo > 0] - x[n_pseudo == 0]),
              .groups = "drop") |>
    mutate(across(starts_with("outbreak_size"), \(x) round(x, 2)))









scatter <- function(sim_outs,
                    .filter_vars = NULL,
                    .filter_conds = NULL,
                    yvar = "outbreak_size",
                    .facet_nrow = NULL,
                    .title = NULL) {

    # sim_outs = sobol_summs; yvar = "outbreak_size"; .facet_nrow = NULL; .title = NULL
    # .filter_vars = "alate_dens"; .filter_conds = 1
    # rm(sim_outs, .filter_vars, .filter_conds, yvar, .facet_nrow, .title, .vary_pars, .df, has_n_pseudo, np_name, .ylab, .y_range, .y_breaks, plot_list)

    # if (yvar == "infect_time") sim_outs <- sim_outs |>
    #     mutate(infect_time = ifelse(is.na(infect_time), 101, infect_time))

    if (!is.null(.filter_vars)) {
        stopifnot(length(.filter_vars) == length(.filter_conds))
        for (i in 1:length(.filter_vars)) {
            sim_outs <- sim_outs |>
                filter(.data[[.filter_vars[i]]] == .filter_conds[i])
        }; rm(i)
    }

    if ("K" %in% colnames(sim_outs)) sim_outs[["K"]] <- sim_outs[["K"]] / 1000

    .vary_pars <- c("n_pseudo", names(vary_pars)) |>
        keep(\(x) x %in% colnames(sim_outs))

    .df <- sim_outs |>
        rename(y = !!yvar) |>
        select(all_of(.vary_pars), y)
    has_n_pseudo <- "n_pseudo" %in% colnames(.df)
    if (has_n_pseudo) {
        .df$n_pseudo <- factor(.df$n_pseudo)
        np_name <- "n_pseudo"
        .ylab <- yvar_desc[[yvar]] |> first_cap()
    } else {
        np_name <- NULL
        .ylab <- paste("Effect of *Pseudomonas* on", yvar_desc[[yvar]])
    }

    if (yvar == "outbreak_size" && has_n_pseudo) {
        .y_range <- c(1, 9)
        .y_breaks <- c(1, 5, 9)
    } else {
        .y_range <- range(.df$y)
        .y_breaks <- waiver()
    }

    plot_list <- .vary_pars |>
        discard(\(x) x == "n_pseudo") |>
        map(\(p){
            # rm(p, z, .pos, plt)
            z <- .df |>
                rename(x = !!p) |>
                select(all_of(c(np_name, "y", "x")))
            # This parameter is a categorical integer, so must be treated differently:
            if (p == "spat_config") {
                z <- z |>
                    mutate(x = x |> as.integer() |> factor())
                .pos <- position_jitter(width = 0.4, height = 0, seed = 871522125)
            } else .pos <- "identity"
            if (has_n_pseudo) {
                plt <- z |>
                    ggplot(aes(x, y, color = n_pseudo, fill = n_pseudo)) +
                    scale_color_manual(pretty_params("n_pseudo", TRUE),
                                       values = c("goldenrod", "dodgerblue")) +
                    scale_fill_manual(pretty_params("n_pseudo", TRUE),
                                      values = c("goldenrod", "dodgerblue")) +
                    guides(color = guide_legend(override.aes = list(alpha = 1)))
            } else {
                plt <- z |>
                    ggplot(aes(x, y))
            }
            plt <- plt +
                geom_point(size = 0.7, alpha = 0.05, position = .pos) +
                scale_y_continuous(limits = .y_range, breaks = .y_breaks) +
                labs(x = pretty_params(p) |> first_cap(), y = .ylab) +
                theme(legend.title = element_markdown(),
                      axis.title.y = element_markdown(),
                      axis.title.x = element_markdown())
            if (p == "spat_config") {
                plt <- plt +
                    stat_summary(fun = "mean", size = 4, geom = "point",
                                 color = "dodgerblue") +
                    theme(legend.position = "none")
            } else {
                plt <- plt +
                    stat_smooth(method = "gam",
                                color = "dodgerblue",
                                formula = y ~ s(x, bs = "cs"),
                                se = TRUE, linewidth = 1)
            }
            return(plt)
        })

        do.call(wrap_plots, c(plot_list, list(nrow = .facet_nrow, guides = "collect", axes = "collect"))) +
            plot_annotation(title = .title, theme = theme(plot.title = element_markdown()))

}







scatter2 <- function(sim_outs,
                     .filter_vars = NULL,
                     .filter_conds = NULL,
                     xvars = c("log_aphids", "log_alates"),
                     yvar = "outbreak_size",
                     .title = NULL) {

    if (!is.null(.filter_vars)) {
        stopifnot(length(.filter_vars) == length(.filter_conds))
        for (i in 1:length(.filter_vars)) {
            sim_outs <- sim_outs |>
                filter(.data[[.filter_vars[i]]] == .filter_conds[i])
        }; rm(i)
    }

    if ("K" %in% colnames(sim_outs)) sim_outs[["K"]] <- sim_outs[["K"]] / 1000

    plot_list <- map(xvars, \(xvar) {
        # if (yvar == "outbreak_size") sim_outs[[yvar]] <- asin(sqrt(sim_outs[[yvar]] / 9))
        # if (yvar == "outbreak_size") sim_outs[[yvar]] <- logit(sim_outs[[yvar]] / 9)
        sim_outs |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            ggplot(aes(.data[[xvar]], .data[[yvar]], color = n_pseudo)) +
            geom_point(size = 0.7, alpha = 0.1) +
            # geom_hline(yintercept = 0, linetype = "22", color = "black", linewidth = 1) +
            stat_smooth(aes(fill = n_pseudo), method = "gam",
                        formula = y ~ s(x, bs = "cs"),
                        se = TRUE, linewidth = 1) +
            scale_color_manual(values = c("goldenrod", "dodgerblue")) +
            scale_fill_manual(values = c("goldenrod", "dodgerblue")) +
            labs(x = yvar_desc[[xvar]] |> str_to_sentence(),
                 y = yvar_desc[[yvar]] |> str_to_sentence()) +
            theme(axis.title.y = element_markdown(),
                  axis.title.x = element_markdown(),
                  plot.title = element_markdown())
    })

    do.call(wrap_plots, c(plot_list, list(nrow = 1, guides = "collect", axes = "collect"))) +
        plot_annotation(title = .title, theme = theme(plot.title = element_markdown()))
}





# # diff_sens_scat_p <-
# scatter(diff_sobol_summs[[2L]], "outbreak_size")
# # scatter(diff_sobol_summs[[1L]], "outbreak_size")
#
# scatter(diff_sobol_summs[[2L]], "p_outbreak")
# # scatter(diff_sobol_summs[[1L]], "p_outbreak")


sum(diff_sobol_summs[["outbreak_size"]] > 0)
sum(diff_sobol_summs[["outbreak_size"]] == 0)
sum(diff_sobol_summs[["outbreak_size"]] < 0)




scatter(sobol_summs, .filter_vars = "alate_dens", .filter_conds = 1, "outbreak_size")
scatter(sobol_summs, .filter_vars = "alate_dens", .filter_conds = 1, "sd_outbreak_size")

sobol_summs |>
    mutate(n_pseudo = factor(n_pseudo)) |>
    ggplot(aes(outbreak_size, sd_outbreak_size, color = n_pseudo, fill = n_pseudo)) +
    geom_point(alpha = 0.1) +
    # stat_smooth(method = "lm", formula = y ~ poly(x,2,raw = TRUE)) +
    stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
                se = TRUE, linewidth = 1) +
    scale_color_manual(pretty_params("n_pseudo", TRUE),
                       values = c("goldenrod", "dodgerblue"),
                       aesthetics = c("color", "fill")) +
    labs(x = yvar_desc[["outbreak_size"]] |> first_cap(),
         y = yvar_desc[["sd_outbreak_size"]] |> first_cap()) +
    facet_wrap(~ alate_dens) +
    theme(axis.title = element_markdown(),
          legend.title = element_markdown())




# QUESTION #1 ----
#' For parameter combos that seemed to result in a negative
#' effect of Pseudomonas, does alate ~ density affect outcomes?

alate_dens_p <- diff_sobol_summs |>
    group_by(combo) |>
    summarize(without = outbreak_size[alate_dens == 0],
              with = outbreak_size[alate_dens == 1]) |>
    ggplot(aes(with, without)) +
    geom_point(alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "22") +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = 0, color = "gray70") +
    coord_equal(xlim = range(diff_sobol_summs$outbreak_size),
                ylim = range(diff_sobol_summs$outbreak_size)) +
    labs(x = "Effect of *Pseudomonas* on outbreak size - alates ~ density",
         y = "Effect of *Pseudomonas* on outbreak size - constant alates") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())
# alate_dens_p
# save_plot("_plots/alate-dens.pdf", alate_dens_p, width = 6, height = 6)

# Answer: Yes, having the alate ~ density effect makes typically makes
# *Pseudomonas* less beneficial to plants






# LEFT OFF #2 ----
#' What parameter values are associated with Pseudomonas being bad for plants?


# Paired plot functions ----

# Paired time series for with and without Pseudomonas
# Note: only takes output summarized by plant
paired_timeseries <- function(sim_df_p,
                              sim_df_np,
                              .rep = NULL,
                              .alate_max = NULL,
                              .aphid_max = NULL,
                              .wasp_max = NULL,
                              .title = waiver(),
                              .tag = waiver()) {

    # sim_df_p = max_diff_sims2; sim_df_np = max_diff_sims02; .rep = NULL; .title = waiver()
    # rm(sim_df_p, sim_df_np, .rep, .title, all_sims, itrans_w, trans_w, itrans_a, trans_a, aphid_breaks, aphid_labels, aphid_ylab, plot_list)

    stopifnot(identical(colnames(sim_df_p), colnames(sim_df_np)))
    stopifnot("n_pseudo" %in% colnames(sim_df_p))
    stopifnot("x" %in% colnames(sim_df_p) && "y" %in% colnames(sim_df_p))
    if (is.null(.rep)) .rep <- 1L
    stopifnot(length(.rep) == 1 && is.numeric(.rep))
    stopifnot(all(.rep %in% sim_df_p$rep) && all(.rep %in% sim_df_np$rep))

    stopifnot(!is.null(sim_df_p[["zeta"]]) && !is.null(sim_df_np[["zeta"]]))

    all_sims <- list(sim_df_p, sim_df_np) |>
        map(\(x) {
            zeta <- x$zeta[[1]]
            x |>
                filter(rep %in% .rep) |>
                mutate(x = ifelse(x == 0, NA, x),
                       y = ifelse(y == 0, NA, y),
                       plant = interaction(x, y),
                       n_pseudo = factor(n_pseudo)) |>
                mutate(aphids = aphids + parasitized) |>
                select(n_pseudo, plant, time, aphids, alates, wasps) |>
                split(~ time, drop = TRUE) |>
                map(\(df_t) {
                    Y <- df_t$wasps[!is.na(df_t$wasps)]
                    stopifnot(length(Y) == 1)
                    z_i <- df_t$aphids[!is.na(df_t$plant)] + df_t$alates[!is.na(df_t$plant)]
                    hat_z_i <- z_i / sum(z_i)
                    q <- length(z_i)
                    gamma_i <- Y * {(1 - zeta) / q + zeta * hat_z_i}
                    df_t_out <- df_t[is.na(df_t$wasps),]
                    df_t_out$wasps <- gamma_i
                    return(df_t_out)
                }) |>
                list_rbind() |>
                pivot_longer(aphids:wasps, names_to = "species",
                             values_to = "density")
        })

    dens_maxes <- c("aphids", "alates", "wasps") |>
        set_names() |>
        map(\(spp) {
            map(all_sims, \(x) x$density[x$species == spp]) |>
                do.call(what = c) |>
                max(na.rm = TRUE)
        })
    if (!is.null(.alate_max)) {
        if (.alate_max < dens_maxes$alates) {
            stop("\nIncrease .alate_max to at least ", dens_maxes$alates)
        }
        dens_maxes$alates <- .alate_max
    }
    if (!is.null(.aphid_max)) {
        if (.aphid_max < dens_maxes$aphids) {
            stop("\nIncrease .aphid_max to at least ", dens_maxes$aphids)
        }
        dens_maxes$aphids <- .aphid_max
    }
    if (!is.null(.wasp_max)) {
        if (.wasp_max < dens_maxes$wasps) {
            stop("\nIncrease .wasp_max to at least ", dens_maxes$wasps)
        }
        dens_maxes$wasps <- .wasp_max
    }

    # convert from wasps --> aphids:
    trans_w <- \(x) x * dens_maxes$aphids / dens_maxes$wasps
    # convert from aphids --> wasps:
    itrans_w <- \(x) x * dens_maxes$wasps / dens_maxes$aphids
    # convert from alates --> aphids:
    trans_a <- \(x) x * dens_maxes$aphids / dens_maxes$alates
    # convert from aphids --> alates:
    itrans_a <- \(x) x * dens_maxes$alates / dens_maxes$aphids

    aphid_breaks <- scales::breaks_extended(n = 4)(c(0, dens_maxes$aphids))
    aphid_labels <- sprintf(paste0("%s (<span style=\"color: ", spp_pal[["alates"]],
                                   ";\">%.1f</span>)"), aphid_breaks,
                            itrans_a(aphid_breaks))
    aphid_ylab <- paste0("Aphid density (<span style=\"color: ", spp_pal[["alates"]],
                        ";\">alate density</span>)")


    plot_list <- map(1:length(all_sims), \(i) {
        sims <- all_sims[[i]]
        tag__ <- .tag
        if (i > 1) tag__ <- waiver()
        sims |>
            mutate(density = case_when(species == "aphids" ~ density,
                                       species == "alates" ~ trans_a(density),
                                       species == "wasps" ~ trans_w(density),
                                       .default = NA)) |>
            ggplot(aes(time, density)) +
            geom_line(aes(color = species), linewidth = 1) +
            facet_wrap( ~ plant, nrow = 3) +
            scale_y_continuous(aphid_ylab,
                               breaks = aphid_breaks,
                               labels = aphid_labels,
                               sec.axis = sec_axis(itrans_w, "Wasp density")) +
            scale_color_manual(values = spp_pal, guide = "none") +
            labs(title = sprintf("%s *Pseudomonas* patches", sims$n_pseudo[[1]]),
                 x = "Time (days)",
                 tag = tag__) +
            coord_cartesian(ylim = c(0, dens_maxes$aphids)) +
            theme(plot.title = element_markdown(hjust = 0.5),
                  legend.title = element_markdown(),
                  plot.tag = element_markdown(size = 16),
                  plot.tag.location = "plot",
                  axis.title.y.left = element_markdown(color = spp_pal[["aphids"]],
                                                       face = "bold"),
                  axis.title.y.right = element_markdown(color = spp_pal[["wasps"]],
                                                        face = "bold"),
                  axis.text.y.left = element_markdown(color = spp_pal[["aphids"]]),
                  axis.text.y.right = element_markdown(color = spp_pal[["wasps"]]))
    })


    do.call(wrap_plots, plot_list) +
        plot_layout(design = "A#B", widths = c(1, 0.05, 1), axes = "collect") +
        plot_annotation(title = .title,
                        theme = theme(plot.title = element_markdown()))
}

# Paired histogram of P(alates) for with and without Pseudomonas
# Note: only takes output summarized by plant
paired_p_alate_hist <- function(sim_df_p,
                                sim_df_np,
                                .rep = NULL,
                                alate_infl = NULL,
                                alate_slope = NULL,
                                .title = waiver(),
                                .tag = waiver()) {

    # sim_df_p = max_diff_sims; sim_df_np = max_diff_sims0; .rep = NULL; alate_infl = NULL; alate_slope = NULL; .title = waiver(); .tag = NULL
    # rm(sim_df_p, sim_df_np, .rep, alate_infl, alate_slope, .title, .tag, all_sims)

    stopifnot(identical(colnames(sim_df_p), colnames(sim_df_np)))
    stopifnot("n_pseudo" %in% colnames(sim_df_p))
    stopifnot("x" %in% colnames(sim_df_p) && "y" %in% colnames(sim_df_p))
    if (is.null(.rep)) .rep <- unique(c(sim_df_p$rep, sim_df_np$rep))
    stopifnot(length(.rep) == 1 && is.numeric(.rep))
    stopifnot(all(.rep %in% sim_df_p$rep) && all(.rep %in% sim_df_np$rep))

    if (is.null(alate_infl)) alate_infl <- pop_info$alate_infl
    if (is.null(alate_slope)) alate_slope <- pop_info$alate_slope

    all_sims <- list(sim_df_p, sim_df_np) |>
        map(\(x) {
            x |>
                filter(rep %in% .rep) |>
                filter(x > 0) |>
                mutate(plant = interaction(x, y),
                       n_pseudo = factor(n_pseudo)) |>
                mutate(z = aphids + parasitized + alates) |>
                select(n_pseudo, rep, plant, time, z) |>
                mutate(p_alate = 1 / (1 + 10^((alate_infl - z) * alate_slope)))
        }) |>
        list_rbind() |>
        mutate(n_pseudo = factor(paste(n_pseudo),
                                 levels = paste(sort(unique(n_pseudo))),
                                 labels = sprintf("%s *Pseudomonas* patches",
                                                  paste(sort(unique(n_pseudo))))))


    all_sims |>
        # filter(p_alate > 1e-2) |>
        mutate(id = interaction(rep, plant, n_pseudo, drop = TRUE)) |>
        # ggplot(aes(time, p_alate)) +
        # geom_line(aes(group = id), linewidth = 1) +
        ggplot(aes((p_alate), after_stat(density))) +
        geom_histogram(aes(fill = n_pseudo), bins = 50) +
        scale_fill_manual(values = c("#1E90FF", "gray60"), guide = "none") +
        facet_wrap( ~ n_pseudo)+
        scale_y_sqrt() +
        labs(x = "Proportion alate offspring", y = "Density") +
        # geom_freqpoly(aes(color = n_pseudo), bins = 20, linewidth = 1) +
        # scale_color_manual(pretty_params("n_pseudo", TRUE),
        #                    values = c("#1E90FF", "gray60")) +
        theme(plot.title = element_markdown(hjust = 0.5),
              legend.title = element_markdown(),
              plot.tag = element_markdown(size = 16),
              plot.tag.location = "plot",
              strip.text = element_markdown(size = 14),
              axis.title.y = element_markdown(),
              axis.text.x = element_markdown())
}



# Paired plots of wasp attack survivals for with and without Pseudomonas
# Note: only takes output summarized by plant
paired_attack_plots <- function(sim_df_p,
                                sim_df_np,
                                .rep = NULL,
                                .title = waiver(),
                                .tag = waiver()) {

    # sim_df_p = max_diff_sims2; sim_df_np = max_diff_sims02; .rep = NULL; .title = waiver(); .tag = waiver()
    # rm(sim_df_p, sim_df_np, .rep, .title, .tag, all_sims, plot_list)

    stopifnot(identical(colnames(sim_df_p), colnames(sim_df_np)))
    stopifnot("n_pseudo" %in% colnames(sim_df_p))
    stopifnot("x" %in% colnames(sim_df_p) && "y" %in% colnames(sim_df_p))
    if (is.null(.rep)) .rep <- 1L
    stopifnot(length(.rep) == 1 && is.numeric(.rep))
    stopifnot(all(.rep %in% sim_df_p$rep) && all(.rep %in% sim_df_np$rep))

    stopifnot(!is.null(sim_df_p[["zeta"]]) && !is.null(sim_df_np[["zeta"]]))

    all_sims <- list(sim_df_p, sim_df_np) |>
        map(\(x) {
            zeta <- x$zeta[[1]]
            x |>
                filter(rep %in% .rep) |>
                mutate(x = ifelse(x == 0, NA, x),
                       y = ifelse(y == 0, NA, y),
                       plant = interaction(x, y),
                       n_pseudo = factor(n_pseudo)) |>
                mutate(aphids = aphids + parasitized) |>
                select(n_pseudo, plant, time, aphids, alates, wasps) |>
                split(~ time, drop = TRUE) |>
                map(\(df_t) {
                    Y <- df_t$wasps[!is.na(df_t$wasps)]
                    stopifnot(length(Y) == 1)
                    z_i <- df_t$aphids[!is.na(df_t$plant)] + df_t$alates[!is.na(df_t$plant)]
                    hat_z_i <- z_i / sum(z_i)
                    q <- length(z_i)
                    gamma_i <- Y * {(1 - zeta) / q + zeta * hat_z_i}
                    df_t_out <- df_t[is.na(df_t$wasps),]
                    df_t_out$wasps <- gamma_i
                    return(df_t_out)
                }) |>
                list_rbind() |>
                pivot_longer(aphids:wasps, names_to = "species",
                             values_to = "density")
        })

    dens_maxes <- c("aphids", "alates", "wasps") |>
        set_names() |>
        map(\(spp) {
            map(all_sims, \(x) x$density[x$species == spp]) |>
                do.call(what = c) |>
                max(na.rm = TRUE)
        })
    if (!is.null(.alate_max)) {
        if (.alate_max < dens_maxes$alates) {
            stop("\nIncrease .alate_max to at least ", dens_maxes$alates)
        }
        dens_maxes$alates <- .alate_max
    }
    if (!is.null(.aphid_max)) {
        if (.aphid_max < dens_maxes$aphids) {
            stop("\nIncrease .aphid_max to at least ", dens_maxes$aphids)
        }
        dens_maxes$aphids <- .aphid_max
    }
    if (!is.null(.wasp_max)) {
        if (.wasp_max < dens_maxes$wasps) {
            stop("\nIncrease .wasp_max to at least ", dens_maxes$wasps)
        }
        dens_maxes$wasps <- .wasp_max
    }

    # convert from wasps --> aphids:
    trans_w <- \(x) x * dens_maxes$aphids / dens_maxes$wasps
    # convert from aphids --> wasps:
    itrans_w <- \(x) x * dens_maxes$wasps / dens_maxes$aphids
    # convert from alates --> aphids:
    trans_a <- \(x) x * dens_maxes$aphids / dens_maxes$alates
    # convert from aphids --> alates:
    itrans_a <- \(x) x * dens_maxes$alates / dens_maxes$aphids

    aphid_breaks <- scales::breaks_extended(n = 4)(c(0, dens_maxes$aphids))
    aphid_labels <- sprintf(paste0("%s (<span style=\"color: ", spp_pal[["alates"]],
                                   ";\">%.1f</span>)"), aphid_breaks,
                            itrans_a(aphid_breaks))
    aphid_ylab <- paste0("Aphid density (<span style=\"color: ", spp_pal[["alates"]],
                         ";\">alate density</span>)")


    plot_list <- map(1:length(all_sims), \(i) {
        sims <- all_sims[[i]]
        tag__ <- .tag
        if (i > 1) tag__ <- waiver()
        sims |>
            mutate(density = case_when(species == "aphids" ~ density,
                                       species == "alates" ~ trans_a(density),
                                       species == "wasps" ~ trans_w(density),
                                       .default = NA)) |>
            ggplot(aes(time, density)) +
            geom_line(aes(color = species), linewidth = 1) +
            facet_wrap( ~ plant, nrow = 3) +
            scale_y_continuous(aphid_ylab,
                               breaks = aphid_breaks,
                               labels = aphid_labels,
                               sec.axis = sec_axis(itrans_w, "Wasp density")) +
            scale_color_manual(values = spp_pal, guide = "none") +
            labs(title = sprintf("%s *Pseudomonas* patches", sims$n_pseudo[[1]]),
                 x = "Time (days)",
                 tag = tag__) +
            coord_cartesian(ylim = c(0, dens_maxes$aphids)) +
            theme(plot.title = element_markdown(hjust = 0.5),
                  legend.title = element_markdown(),
                  plot.tag = element_markdown(size = 16),
                  plot.tag.location = "plot",
                  axis.title.y.left = element_markdown(color = spp_pal[["aphids"]],
                                                       face = "bold"),
                  axis.title.y.right = element_markdown(color = spp_pal[["wasps"]],
                                                        face = "bold"),
                  axis.text.y.left = element_markdown(color = spp_pal[["aphids"]]),
                  axis.text.y.right = element_markdown(color = spp_pal[["wasps"]]))
    })


    do.call(wrap_plots, plot_list) +
        plot_layout(design = "A#B", widths = c(1, 0.05, 1), axes = "collect") +
        plot_annotation(title = .title,
                        theme = theme(plot.title = element_markdown()))
}


# For this function, it's assumed that you want to summarize by rep if it's
# not summarized by anything (and `.all = FALSE`) and that you want to
# summarize across all reps if it's already summarized by rep (or if `.all = TRUE`)
calc_metrics <- function(sims, .all = TRUE) {
    if ("x" %in% colnames(sims)) {
        out <- sims |>
            filter(!is.na(wasps)) |>
            group_by(rep) |>
            summarize(outbreak_size = max(virus),
                      log_aphids = mean(log10(1 + aphids + parasitized)),
                      log_alates = mean(log10(1 + alates)),
                      log_wasps = mean(log10(1 + wasps)))
        if (!.all) return(out)
    } else {
        out <- sims
    }

    out <- out |>
        summarize(across(all_of(c("outbreak_size", "log_aphids",
                                  "log_alates", "log_wasps")),
                         mean))

    return(out)
}


paired_target_sims <- function(.combo, .n_sims = 4L, .summ = "none",
                               ...) {
    # .combo = 7112L
    # rm(.combo, arg_list, other_args, out)
    arg_list <- sobol_summs |>
        filter(combo == .combo) |>
        filter(alate_dens == 1, n_pseudo > 0) |>
        select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
        as.list()
    other_args <- list(...)
    stopifnot(length(other_args) >= 1L)
    stopifnot(!is.null(names(other_args)) && !any(names(other_args) == ""))
    stopifnot(all(names(other_args) %in% c(names(formals(one_combo)),
                                           names(formals(sim_plantscape)),
                                           names(formals(make_insect_ptr)))))
    out <- crossing(others = c(FALSE, TRUE),
             n_pseudo = c(0L, arg_list$n_pseudo)) |>
        pmap(\(others, n_pseudo) {
            args <- c(arg_list, list(n_sims = .n_sims, summ = .summ))
            if (others) {
                for (n in names(other_args)) args[[n]] <- other_args[[n]]
            }
            if (n_pseudo == 0L) args$n_pseudo <- 0L
            .sims <- do.call(one_combo, args)
            .out <- select(.sims, n_pseudo, rep:last_col())
            for (n in names(other_args)) {
                if (n %in% colnames(.sims)) {
                    .out[[n]] <- .sims[[n]]
                } else if (n %in% names(formals(sim_plantscape))) {
                    .out[[n]] <- formals(sim_plantscape)[[n]]
                } else if (n %in% names(formals(make_insect_ptr))) {
                    .out[[n]] <- formals(make_insect_ptr)[[n]]
                } else stop("\nCannot find ", n)
            }
            return(.out)
        }) |>
        list_rbind() |>
        select(all_of(names(other_args)), n_pseudo, everything())

    return(out)

}



one_combo_par_manip_simmer <- function(.combo) {

    # .combo = 5624
    # rm(.combo, paras, sim_list)

    pars <- list(Y0 = 1:20,
                 mean_N = 1:20 * 10,
                 sd_N = 0:20 * 5,
                 K = 1:25 * 1000,
                 virus_attract = seq(1, 5, 0.25),
                 pseudo_repel = seq(1, 5, 0.25),
                 pseudo_surv = seq(0.85, 1, 0.01),
                 zeta = seq(0, 1, 0.05),
                 spat_config = 0:4) |>
        map(\(x) round(x, 2))

    sim_list <- names(pars) |>
        map(\(x_name) {
            # x_name = "K"
            # rm(x_name, df_og, df_x, df_max_diff)
            df_og <- sobol_summs |>
                filter(combo == .combo) |>
                filter(alate_dens == 1) |>
                select(n_pseudo, all_of(x_name), outbreak_size) |>
                mutate(n_pseudo = factor(n_pseudo))
            df_x <- pars[[x_name]] |>
                map(\(x) {
                    args <- sobol_summs |>
                        filter(combo == .combo) |>
                        filter(alate_dens == 1, n_pseudo > 0) |>
                        select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
                        as.list()
                    args[[x_name]] <- x
                    args <- c(args, list(n_sims = 1000, summ = "all"))
                    np <- args[["n_pseudo"]]
                    .sims <- do.call(one_combo, args) |>
                        getElement("outbreak_size") |>
                        mean()
                    args[["n_pseudo"]] <- 0L
                    .sims0 <- do.call(one_combo, args) |>
                        getElement("outbreak_size") |>
                        mean()
                    return(tibble(!! x_name := x,
                                  n_pseudo = factor(c(0L, np)),
                                  outbreak_size = c(.sims0, .sims)))
                }) |>
                list_rbind() |>
                mutate(combo = .combo) |>
                select(combo, everything())
            return(df_x)
        })

    return(sim_list)
}



one_combo_par_manip_plotter <- function(sim_list) {

    .combo <- sim_list[[1]][["combo"]][[1]]

    plot_list <- sim_list |>
        map(\(df_x) {
            # x_name = "K"
            # rm(x_name, df_og, df_x, df_max_diff)
            x_name <- colnames(df_x)[2]
            df_og <- sobol_summs |>
                filter(combo == .combo) |>
                filter(alate_dens == 1) |>
                select(n_pseudo, all_of(x_name), outbreak_size) |>
                mutate(n_pseudo = factor(n_pseudo))
            df_max_diff <- df_x |>
                group_by(across(all_of(x_name))) |>
                summarize(dos = outbreak_size[as.integer(n_pseudo) == 2L] -
                              outbreak_size[as.integer(n_pseudo) == 1L],
                          .groups = "drop") |>
                filter(dos == max(dos)) |>
                slice(1) |>
                getElement(x_name)
            df_x |>
                ggplot(aes(.data[[x_name]], outbreak_size, color = n_pseudo)) +
                geom_hline(yintercept = c(1, 9), color = "gray70") +
                geom_vline(xintercept = df_max_diff, color = "gray70",
                           linetype = "22") +
                geom_point() +
                geom_line() +
                geom_point(data = df_og, size = 3, shape = 8) +
                labs(x = pretty_params(x_name) |> first_cap(),
                     y = yvar_desc[["outbreak_size"]] |> first_cap()) +
                scale_y_continuous(breaks = c(1, 5, 9)) +
                coord_cartesian(ylim = c(1, 9)) +
                scale_color_manual(pretty_params("n_pseudo", TRUE),
                                   values = c("goldenrod", "dodgerblue")) +
                guides(color = guide_legend(override.aes = list(shape = 19))) +
                theme(axis.title.x = element_markdown(),
                      axis.title.y = element_markdown(),
                      legend.title = element_markdown())
        })

    plot_list |>
        do.call(what = wrap_plots) +
        plot_layout(guides = "collect", axes = "collect") +
        plot_annotation(title = sprintf("combo %i", .combo))
}





# *LEFT OFF 5-Nov ----
diff_sobol_summs |>
    filter(alate_dens == 1) |>
    arrange(desc(outbreak_size)) |>
    select(combo, all_of(names(vary_pars)), outbreak_size)


# Highest outbreak size difference:
combo_sims1 <- one_combo_par_manip_simmer(7112)
# higher virus_attract and pseudo_repel:
combo_sims2 <- one_combo_par_manip_simmer(5624)
# higher zeta, spat_config = 3 (vs 1 or 2):
combo_sims3 <- one_combo_par_manip_simmer(10940)


combo_p1 <- one_combo_par_manip_plotter(combo_sims1)
combo_p2 <- one_combo_par_manip_plotter(combo_sims2)
combo_p3 <- one_combo_par_manip_plotter(combo_sims3)

combo_p1
combo_p2
combo_p3


# save_plot("_plots/combo-plot1.png", combo_p1, width = 8, height = 5, dpi = 150)
# save_plot("_plots/combo-plot2.png", combo_p2, width = 8, height = 5, dpi = 150)
# save_plot("_plots/combo-plot3.png", combo_p3, width = 8, height = 5, dpi = 150)



max_diff_args <- list(alate_dens = 1L,
                      Y0 = 2.5,
                      mean_N = 25,
                      sd_N = 0,
                      K = 23e3,
                      virus_attract = 1,
                      pseudo_repel = 1,
                      pseudo_surv = 0.85,
                      zeta = 0.0,
                      spat_config = 1L,
                      n_sims = 1e3L,
                      n_pseudo = 3L,
                      summ = "none")
max_diff_sims <- do.call(one_combo, max_diff_args)
max_diff_sims0 <- do.call(one_combo, list_assign(max_diff_args, n_pseudo = 0L))

paired_timeseries(max_diff_sims, max_diff_sims0, .title = "Original")


max_diff_sims2 <- do.call(one_combo, list_assign(max_diff_args, pseudo_surv = 0.95))
max_diff_sims02 <- do.call(one_combo, list_assign(max_diff_args, pseudo_surv = 0.95, n_pseudo = 0L))

bind_rows(calc_metrics(max_diff_sims), calc_metrics(max_diff_sims0),
          calc_metrics(max_diff_sims2), calc_metrics(max_diff_sims02)) |>
    mutate(n_pseudo = c(3, 0, 3, 0),
           type = c(rep("orig", 2), rep("modified", 2))) |>
    select(type, n_pseudo, everything())

paired_timeseries(max_diff_sims, max_diff_sims0, .tag = "Original",
                  .aphid_max = 1389, .alate_max = 94, .wasp_max = 39) /
    paired_timeseries(max_diff_sims2, max_diff_sims02,
                      .tag = serify("", "&psi;", " = 0.95"),
                      .aphid_max = 1389, .alate_max = 94, .wasp_max = 39) &
    theme(panel.grid.major.y = element_line(color = "gray80"))











max_diff_sims_ob <- max_diff_sims |>
    select(rep:last_col()) |>
    filter(!is.na(wasps)) |>
    group_by(rep) |>
    summarize(outbreak_size = max(virus)) |>
    getElement("outbreak_size") |>
    mean()
max_diff_sims0_ob <- max_diff_sims0 |>
    select(rep:last_col()) |>
    filter(!is.na(wasps)) |>
    group_by(rep) |>
    summarize(outbreak_size = max(virus)) |>
    getElement("outbreak_size") |>
    mean()

mean(max_diff_sims_ob); mean(max_diff_sims0_ob)
print_diff_mean_w_boot_ci(max_diff_sims_ob, max_diff_sims0_ob)





# largest so far: 5.29   (5.15 - 5.42)
#' alate_dens = 1L,
#' Y0 = 2.5,
#' mean_N = 25,
#' sd_N = 0,
#' K = 23e3,
#' virus_attract = 1,
#' pseudo_repel = 1,
#' pseudo_surv = 0.85,
#' zeta = 0.0,
#' spat_config = 1L


print_diff_mean_w_boot_ci <- function(a, b, B = 2000L, alpha = 0.05) {
    n <- length(a)
    m <- length(b)
    stopifnot(n == m)
    boots <- sapply(1:B, \(i) {
        ai <- sample(a, replace = TRUE)
        bi <- sample(b, replace = TRUE)
        return(mean(ai) - mean(bi))
    })
    boot_ci <- quantile(boots, c(alpha, 1-alpha))
    cat(sprintf("%.2f\t(%s)\n", mean(a) - mean(b),
                paste(format(boot_ci, digits = 3), collapse = " - ")))
    invisible(NULL)
}

# 5.071




paired_target_sims(7112L, .n_sims = 1000, .summ = "all",
                   pseudo_repel = 2,
                   zeta = 0.2) |>
    split(~ zeta +
              pseudo_repel +
              n_pseudo, drop = TRUE) |>
    map(\(x) calc_metrics(x) |>
            mutate(zeta = x$zeta[[1]],
                   pseudo_repel = x$pseudo_repel[[1]],
                   n_pseudo = x$n_pseudo[[1]])) |>
    list_rbind() |>
    arrange(zeta,
            pseudo_repel,
            n_pseudo) |>
    select(zeta,
           pseudo_repel,
           n_pseudo, everything())


# With just changing zeta = 0.2

# # A tibble: 4 × 6
#     zeta n_pseudo outbreak_size log_aphids log_alates log_wasps
#    <dbl>    <int>         <dbl>      <dbl>      <dbl>     <dbl>
# 1 0.0322        0          4.24       5.04       1.81      3.01
# 2 0.0322        3          7.27       5.22       2.06      3.11
# 3 0.2           0          4.25       5.03       1.80      3.01
# 4 0.2           3          5.92       5.09       1.90      3.00


sims <- paired_target_sims(7112L, zeta = 0.5)




# set.seed(1032439295)
sim_df <- sobol_summs |>
    filter(combo == i) |>
    filter(n_pseudo > 0, alate_dens == 1) |>
    select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
    as.list() |>
    c(list(n_sims = 4, summ = "none")) |>
    do.call(what = one_combo) |>
    select(n_pseudo, rep:wasps)
sim0_df <- sobol_summs |>
    filter(combo == i) |>
    filter(n_pseudo == 0, alate_dens == 1) |>
    select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
    as.list() |>
    c(list(n_sims = 4, summ = "none")) |>
    do.call(what = one_combo) |>
    select(n_pseudo, rep:wasps)

sim_df |> calc_metrics()
sim0_df |> calc_metrics()

paired_timeseries(sim_df, sim0_df, .title = sprintf("combo %i", i))
paired_timeseries(sim_df, sim0_df, TRUE, 4, .title = sprintf("combo %i", i))



diff_sobol_summs |>
    filter(alate_dens == 1, spat_config == 1) |>
    ggplot(aes(zeta, Y0)) +
    geom_point(aes(color = outbreak_size)) +
    scale_color_viridis_c()

diff_sobol_summs |>
    filter(alate_dens == 1, spat_config == 1) |>
    ggplot(aes(zeta, outbreak_size)) +
    geom_point(alpha = 0.1) +
    # stat_smooth(method = "gam",
    #             formula = y ~ s(x, bs = "cs"),
    #             se = TRUE, linewidth = 1) +
    labs(x = pretty_params("zeta"), y = yvar_desc[["outbreak_size"]]) +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown())




# LEFT OFF #3 ----
# Why does the effect of Pseudomonas on P(outbreak) not coincide closely with
# its effect on outbreak size?


# names(vary_pars)
# [1] "Y0"            "mean_N"        "sd_N"          "K"             "virus_attract"
# [6] "pseudo_repel"  "pseudo_surv"   "zeta"          "spat_config"

dd <- diff_sobol_summs[[2L]] |>
    filter(spat_config == 0) |>
    mutate(spat_config = factor(spat_config),
           resid = residuals(lm(outbreak_size ~ p_outbreak)))

dd |>
    ggplot(aes(log(outbreak_size), p_outbreak)) +
    geom_point() +
    # facet_wrap(~ spat_config) +
    scale_color_viridis_c()

mod <- lm(resid ~ K * Y0, dd)
mod |> summary()


dd |>
    # group_by(spat_config) |>
    # mutate(resid = residuals(lm(outbreak_size ~ p_outbreak))) |>
    # ungroup() |>
    mutate(resid = residuals(mod)) |>
    ggplot(aes(resid, Y0)) +
    geom_point()




# sens_diff_scat_p <- scatter(diff_sobol_summs, .title = "")
# save_plot("_plots/sens-diff-scatter.pdf", sens_diff_scat_p,
#           width = 8, height = 5)
# save_plot("_plots/sens-diff-scatter.png", sens_diff_scat_p, dpi = 150,
#           width = 8, height = 5)


# sens_scat_p <- scatter(sobol_summs[[2L]], "outbreak_size", .title = "")
# save_plot("_plots/sens-scatter.pdf", sens_scat_p, width = 8, height = 5)
#
# sens_scat0_p <- scatter(sobol_summs[[1L]], "outbreak_size", .title = "")
# save_plot("_plots/sens-scatter-b0.pdf", sens_scat0_p, width = 8, height = 5)
#
# # scatter2(sobol_summs[[1L]])
# sens_scat2_p <- scatter2(sobol_summs[[2L]], .title = "")
#
# save_plot("_plots/sens-scatter2.pdf", sens_scat2_p, width = 8, height = 5)



spat_config_p <- sobol_summs[[2]] |>
    mutate(spat_config = factor(spat_config), n_pseudo = factor(n_pseudo)) |>
    ggplot(aes(spat_config, outbreak_size, color = n_pseudo)) +
    geom_violin() +
    geom_hline(data = sobol_summs[[2]] |>
                   mutate(spat_config = factor(spat_config), n_pseudo = factor(n_pseudo)) |>
                   filter(spat_config == 1) |>
                   group_by(spat_config, n_pseudo) |>
                   summarize(outbreak_size = mean(outbreak_size), .groups = "drop"),
               aes(yintercept = outbreak_size, color = n_pseudo),
               linetype = "22") +
    stat_summary(fun = "mean", size = 4, geom = "point") +
    labs(x = pretty_params("spat_config") |> str_to_sentence(),
         y = yvar_desc[["outbreak_size"]] |> str_to_sentence()) +
    scale_color_manual(pretty_params("n_pseudo", TRUE),
                       values = c("goldenrod", "dodgerblue")) +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown(),
          legend.title = element_markdown(),
          plot.title = element_markdown())

# save_plot("_plots/spat-config.pdf", spat_config_p, width = 6, height = 4)

diff_spat_config_p <- diff_sobol_summs[[2]] |>
    mutate(spat_config = factor(spat_config)) |>
    ggplot(aes(spat_config, outbreak_size)) +
    geom_violin() +
    geom_hline(data = diff_sobol_summs[[2]] |>
                   mutate(spat_config = factor(spat_config)) |>
                   filter(spat_config == 1) |>
                   summarize(outbreak_size = mean(outbreak_size), .groups = "drop"),
               aes(yintercept = outbreak_size),
               linetype = "22") +
    stat_summary(fun = "mean", size = 4, geom = "point") +
    labs(x = pretty_params("spat_config") |> str_to_sentence(),
         y = paste("Effect of *Pseudomonas* on", yvar_desc[["outbreak_size"]])) +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown(),
          legend.title = element_markdown(),
          plot.title = element_markdown())

# save_plot("_plots/spat-config-diff.pdf", diff_spat_config_p, width = 6, height = 4)




multi_scatter <- function(sim_outs, np = NULL, .yvar = "outbreak_size", .N = NULL, ...) {

    # sim_outs = diff_sobol_summs[[1]]; np = NULL; .yvar = "outbreak_size"; .N = N
    # rm(sim_outs, np, .yvar, .N, has_n_pseudo, params, Y, dt, out, output, y_range, .title, plot_list)

    if ("K" %in% colnames(sim_outs)) sim_outs[["K"]] <- sim_outs[["K"]] / 1000

    has_n_pseudo <- "n_pseudo" %in% colnames(sim_outs)
    if (!is.null(np) && has_n_pseudo) {
        sim_outs <- sim_outs |> filter(n_pseudo == np)
    }

    if (is.null(.N)) .N <- nrow(sim_outs)

    stopifnot(nrow(sim_outs) >= .N)

    params <- names(vary_pars)

    Y <- sim_outs[[.yvar]]
    dt <- sim_outs[,params]
    out <- t(combn(params, 2))
    output <- map(1:nrow(out), \(i) {
        cols <- out[i, ]
        if (.N < nrow(dt)) {
            idx <- sample.int(nrow(dt), .N)
            ddt <- dt[idx, cols]
            Yi <- Y[idx]
        } else {
            ddt <- dt[, cols]
            Yi <- Y
        }
        ddt |>
            set_names(c("xvar", "yvar")) |>
            mutate(x = cols[1], y = cols[2], output = Yi)
    }) |>
        list_rbind()

    if (.yvar == "outbreak_size" && has_n_pseudo) {
        y_range <- c(1L, 9L)
    } else y_range <- output$output |> range() |> floor() |> (\(x) x + c(0, 1))()

    .title <- colnames(sim_outs) |>
        keep(\(x) x %in% names(struct_pars)) |>
        map(\(x) sprintf("%s = %.3g", pretty_params(x, TRUE),
                         sim_outs[[x]][[1]])) |>
        c(list(collapse = ";")) |>
        do.call(what = paste)


    plot_list <- output |>
        mutate(across(x:y, \(x) factor(pretty_params(x, TRUE),
                                       levels = pretty_params(params, TRUE)))) |>
        split(~ x + y, drop = TRUE) |>
        map(\(z) {
            z |>
                ggplot(aes(xvar, yvar, color = output)) +
                geom_point(size = 0.5) +
                # scale_colour_gradientn(colours = grDevices::terrain.colors(10),
                scale_color_viridis_c(yvar_desc[[.yvar]], limits = y_range) +
                scale_x_continuous(z$x[[1]], breaks = scales::pretty_breaks(n = 3)) +
                scale_y_continuous(z$y[[1]], breaks = scales::pretty_breaks(n = 3)) +
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.title.x = element_markdown(family = "serif"),
                      axis.title.y = element_markdown(family = "serif"),
                      legend.background = element_rect(fill = "transparent", color = NA),
                      legend.key = element_rect(fill = "transparent", color = NA))
        })

    do.call(wrap_plots, c(plot_list, list(guides = "collect"))) +
        plot_annotation(title = .title,
                        theme = theme(plot.title = element_markdown())) &
        theme(legend.position = "top", ...)
}



# multi_scatter(sobol_summs[[2L]], np = 0L, .N = 2^8, axis.text = element_markdown(size = 6),
#               axis.title = element_markdown(family = "serif", size = 8))
# multi_scatter(sobol_summs[[2L]], np = 3L, .N = 2^8, axis.text = element_markdown(size = 6),
#               axis.title = element_markdown(family = "serif", size = 8))




set.seed(1998643658)
sobol_inds_np3 <- map(sobol_summs, \(ss) {
    Y <- ss[["outbreak_size"]][ss$n_pseudo > 0]
    ind <- sobol_indices(Y = Y, N = N, params = pretty_params(names(vary_pars)),
                         boot = TRUE, R = 2000)
    return(ind)
})

set.seed(512925036)
sobol_inds_diff <- map(diff_sobol_sims, \(ss) {
    Y <- ss[["outbreak_size"]]
    ind <- sobol_indices(Y = Y, N = N, params = pretty_params(names(vary_pars)),
                         boot = TRUE, R = 2000)
    return(ind)
})


sobol_inds_p1 <- plot(sobol_inds_np3[[2]]) +
    labs(y = "Sobol' index (outbreak size with<br>three *Pseudomonas*)") +
    theme(axis.text.y = element_markdown(),
          axis.title.x = element_markdown()) +
    coord_flip()
sobol_inds_p2 <- plot(sobol_inds_diff[[2]]) +
    labs(y = "Sobol' index (effect of *Pseudomonas*<br>on outbreak size)") +
    theme(axis.text.y = element_markdown(),
          axis.title.x = element_markdown()) +
    coord_flip()

sobol_inds_p <- (sobol_inds_p1 | sobol_inds_p2) +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "top")

save_plot("_plots/sobol-inds.pdf", sobol_inds_p, width = 8, height = 6)






