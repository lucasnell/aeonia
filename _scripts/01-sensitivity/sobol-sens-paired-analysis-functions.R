#'
#' Plotting + simulation functions for small-scale sensitivity via Sobol indices
#'
#'




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



# Paired plot functions ----

#' Add wasp density and attack survivals to a dataframe.
#' The attack survivals will only be added if the dataframe has
#' output by stage.
#' It's assumed that the dataframe is already filtered by rep.
org_timeseries_vars <- function(x) {
    # x = stg_max_diff_sims |> filter(rep == 1)
    # rm(x, zeta, a, h, k, R, attack_pars, out_survs)
    zeta <- x$zeta[[1]]
    attack_pars <- list()
    for (p in c("a", "h", "k", "R")) {
        if (p %in% colnames(x)) attack_pars[[p]] <- x[[p]][[1]]
        else attack_pars[[p]] <- pop_info[[p]]
    }
    a <- attack_pars$a
    h <- attack_pars$h
    k <- attack_pars$k
    R <- attack_pars$R
    out_survs <- "aphids_juv" %in% colnames(x)
    if (out_survs) {
        x <- x |>
            mutate(aphids = aphids_juv + aphids_adu,
                   alates = alates_juv + alates_adu)
    }
    x |>
        mutate(plant = interaction(x, y),
               n_pseudo = factor(n_pseudo)) |>
        select(-x, -y) |>
        select(n_pseudo, rep, plant, time, starts_with("aphids"),
               starts_with("alates"), parasitized:wasps) |>
        split(~ n_pseudo + rep + time, drop = TRUE) |>
        map(\(df_t) {
            # rm(df_t, Y, q, z_i, hat_z_i, x_i, Amat, Nmat)
            Y <- df_t$wasps[!is.na(df_t$wasps)]
            stopifnot(length(Y) == 1)
            df_t <- df_t |> filter(!is.na(plant))
            q <- nrow(df_t)
            z_i <- df_t |>
                mutate(z = aphids + alates + parasitized) |>
                getElement("z")
            hat_z_i <- z_i / sum(z_i)
            df_t$wasps <- Y * {(1 - zeta) / q + zeta * hat_z_i}
            if (out_survs) {
                x_i <- df_t |>
                    mutate(x = aphids + alates) |>
                    getElement("x")
                # Attack survivals by plant and stage:
                Amat <- map(R, \(R_i) (1 + R_i * a * df_t$wasps /
                                           (k * (h * x_i + 1)))^(-k)) |>
                    do.call(what  = cbind)
                # Relative aphid abundance by stage:
                Nmat <- df_t |>
                    select(aphids_juv, aphids_adu, alates_juv, alates_adu) |>
                    mutate(across(everything(),
                                  \(x) {
                                      x / (aphids_juv + aphids_adu +
                                               alates_juv + alates_adu)
                                  })) |>
                    as.matrix()
                # Attack survivals weighted by abundance:
                df_t$A <- rowSums(Amat * Nmat)
                df_t <- df_t |>
                    select(-ends_with("_juv"), -ends_with("_adu"))
            }
            return(df_t)
        }) |>
        list_rbind()
}

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

    # sim_df_p = stg_max_diff_sims; sim_df_np = stg_max_diff_sims0; .rep = NULL;
    # .alate_max = NULL; .aphid_max = NULL; .wasp_max = NULL
    # .title = waiver(); .tag = waiver()
    # rm(sim_df_p, sim_df_np, .rep, .alate_max, .aphid_max, .wasp_max, .title)
    # rm(all_sims, trans, itrans, aphid_breaks)
    # rm(aphid_labels, aphid_ylab, plot_list)

    stopifnot(identical(colnames(sim_df_p), colnames(sim_df_np)))
    stopifnot("n_pseudo" %in% colnames(sim_df_p))
    stopifnot("x" %in% colnames(sim_df_p) && "y" %in% colnames(sim_df_p))
    if (is.null(.rep)) .rep <- 1L
    stopifnot(length(.rep) == 1 && is.numeric(.rep))
    stopifnot(all(.rep %in% sim_df_p$rep) && all(.rep %in% sim_df_np$rep))

    stopifnot(!is.null(sim_df_p[["zeta"]]) && !is.null(sim_df_np[["zeta"]]))

    all_sims <- list(sim_df_p, sim_df_np) |>
        map(\(x) x |>
                filter(rep %in% .rep) |>
                org_timeseries_vars() |>
                mutate(aphids = aphids + parasitized) |>
                select(-parasitized, -mummies) |>
                pivot_longer(aphids:last_col(), names_to = "species",
                             values_to = "density"))

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

    # convert from wasps or alates --> aphids:
    trans <- c("wasps", "alates") |>
        set_names() |>
        map(\(n) {
            denom <- paste0('dens_maxes[["', n, '"]]')
            eval(parse(text = paste0("function(x) { return(" ,
                                     "x * dens_maxes$aphids / ", denom, ")}")))
        })
    # convert from aphids --> wasps or alates:
    itrans <- c("wasps", "alates") |>
        set_names() |>
        map(\(n) {
            numer <- paste0('dens_maxes[["', n, '"]]')
            eval(parse(text = paste0("function(x) { return(" ,
                                     "x *", numer, "/ dens_maxes$aphids)}")))
        })

    aphid_breaks <- scales::breaks_extended(n = 4)(c(0, dens_maxes$aphids))
    aphid_labels <- sprintf(paste0("%s (<span style=\"color: ", spp_pal[["alates"]],
                                   ";\">%.1f</span>)"), aphid_breaks,
                            itrans$alates(aphid_breaks))
    aphid_ylab <- paste0("Aphid density (<span style=\"color: ", spp_pal[["alates"]],
                        ";\">alate density</span>)")


    plot_list <- map(1:length(all_sims), \(i) {
        sims <- all_sims[[i]]
        tag__ <- .tag
        if (i > 1) tag__ <- waiver()
        sims |>
            mutate(density = case_when(species == "aphids" ~ density,
                                       species == "alates" ~ trans$alates(density),
                                       species == "wasps" ~ trans$wasps(density),
                                       species == "A" ~ density * dens_maxes$aphids,
                                       .default = NA)) |>
            ggplot(aes(time, density)) +
            geom_line(aes(color = species), linewidth = 1) +
            facet_wrap( ~ plant, nrow = 3) +
            scale_y_continuous(aphid_ylab,
                               breaks = aphid_breaks,
                               labels = aphid_labels,
                               sec.axis = sec_axis(itrans$wasps, "Wasp density")) +
            scale_color_manual(values = c(spp_pal, A = "gray70"), guide = "none") +
            labs(title = sprintf("%s *Pseudomonas* patches", sims$n_pseudo[[1]]),
                 x = "Time (days)",
                 tag = tag__) +
            coord_cartesian(ylim = c(0, dens_maxes$aphids)) +
            theme(plot.title = element_markdown(hjust = 0.5),
                  legend.title = element_markdown(),
                  plot.tag = element_markdown(size = 16),
                  plot.tag.location = "plot",
                  panel.grid.major.y = element_line(color = "gray80"),
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
    stopifnot(is.numeric(.rep) && length(.rep) == 1)
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
# Note: only takes output summarized by plant, and with stages output
paired_attack_plots <- function(sim_df_p,
                                sim_df_np,
                                .rep = NULL,
                                .y_min = NA,
                                .title = waiver(),
                                .tag = waiver()) {

    # sim_df_p = stg_max_diff_sims; sim_df_np = stg_max_diff_sims0
    # .rep = NULL; .y_min = NA; .title = waiver(); .tag = waiver()
    # rm(sim_df_p, sim_df_np, .rep, .y_min, .title, .tag, all_sims)

    stopifnot(identical(colnames(sim_df_p), colnames(sim_df_np)))
    stopifnot("n_pseudo" %in% colnames(sim_df_p))
    stopifnot("aphids_juv" %in% colnames(sim_df_p))
    stopifnot("x" %in% colnames(sim_df_p) && "y" %in% colnames(sim_df_p))
    if (is.null(.rep)) .rep <- 1L
    stopifnot(is.numeric(.rep))
    stopifnot(all(.rep %in% sim_df_p$rep) && all(.rep %in% sim_df_np$rep))

    stopifnot(!is.null(sim_df_p[["zeta"]]) && !is.null(sim_df_np[["zeta"]]))

    .alpha <- ifelse(length(.rep) == 1, 1, 0.25)


    all_sims <- bind_rows(sim_df_p, sim_df_np) |>
        filter(rep %in% .rep) |>
        org_timeseries_vars()

    if (!is.na(.y_min) && .y_min > min(all_sims$A)) {
        stop("\nIncrease .y_min to <= ", min(all_sims$A), call. = FALSE)
    }

    all_sims |>
        mutate(id = interaction(n_pseudo, rep, plant, drop = TRUE)) |>
        ggplot(aes(time, A)) +
        geom_line(aes(group = id, color = n_pseudo), linewidth = 1,
                  alpha = .alpha) +
        labs(title = .title,
             x = "Time (days)",
             y = "Parasitoid attack survival",
             tag = .tag) +
        facet_wrap( ~ plant, nrow = 3) +
        scale_color_manual(pretty_params("n_pseudo", TRUE),
                           values = np_pal) +
        coord_cartesian(ylim = c(.y_min, 1)) +
        # ggplot(aes((A), after_stat(density))) +
        # geom_histogram(aes(fill = n_pseudo), bins = 50) +
        # scale_fill_manual(values = c("#1E90FF", "gray60"), guide = "none") +
        # scale_y_sqrt() +
        # labs(x = "Parasitoid attack survival", y = "Density") +
        # facet_wrap( ~ n_pseudo) +
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
    # rm(.combo, pars, sim_list)

    pars <- list(Y0 = 1:20,
                 mean_N = 1:20 * 10,
                 sd_N = 0:20 * 5,
                 K = 1:25 * 1000,
                 # virus_attract = seq(1, 5, 0.25),
                 # pseudo_repel = seq(1, 5, 0.25),
                 virus_attract = 0:19 / 2 + 1,
                 pseudo_repel = 0:19 / 2 + 1,
                 pseudo_surv = seq(0.85, 1, 0.01),
                 zeta = seq(0, 1, 0.05),
                 spat_config = 0:4) |>
        map(\(x) round(x, 2))

    sim_list <- names(pars) |>
        map(\(x_name) {
            # x_name = "K"
            # rm(x_name, df_x)
            df_x <- pars[[x_name]] |>
                map(\(x) {
                    # x = pars[[x_name]][[1]]
                    # rm(x, args, .sims, .sims0, out)
                    args <- sobol_summs |>
                        filter(combo == .combo) |>
                        filter(alate_dens == 1, n_pseudo > 0) |>
                        select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
                        as.list()
                    args[[x_name]] <- x
                    args <- c(args, list(n_sims = 1000, summ = "all"))
                    .sims <- do.call(one_combo, args) |>
                        select(rep:last_col()) |>
                        select(-rep) |>
                        summarize(across(everything(), mean)) |>
                        mutate(n_pseudo = args[["n_pseudo"]])
                    args[["n_pseudo"]] <- 0L
                    .sims0 <- do.call(one_combo, args) |>
                        select(rep:last_col()) |>
                        select(-rep) |>
                        summarize(across(everything(), mean)) |>
                        mutate(n_pseudo = args[["n_pseudo"]])
                    out <- bind_rows(.sims, .sims0) |>
                        mutate(!! x_name := x,
                               n_pseudo = factor(n_pseudo)) |>
                        select(!! x_name, n_pseudo, everything())
                    return(out)
                }) |>
                list_rbind() |>
                mutate(combo = .combo) |>
                select(combo, everything())
            return(df_x)
        })

    return(sim_list)
}


one_combo_2par_manip_simmer <- function(.combo, .pars, progress_ = .prog_args) {

    # .combo = 5624; .pars = c("Y0", "mean_N"); progress_ = .prog_args
    # rm(.combo, .pars, progress_, par_perms, args, sim_list, out)

    stopifnot(length(.pars) == 2L && is.character(.pars))
    stopifnot(all(.pars %in% names(vary_pars)))

    par_perms <- list(Y0 = 1:20,
                 mean_N = 1:20 * 10,
                 sd_N = 0:20 * 5,
                 K = 1:25 * 1000,
                 # virus_attract = seq(1, 5, 0.25),
                 # pseudo_repel = seq(1, 5, 0.25),
                 virus_attract = 0:19 / 2 + 1,
                 pseudo_repel = 0:19 / 2 + 1,
                 pseudo_surv = seq(0.85, 1, 0.01),
                 zeta = seq(0, 1, 0.05),
                 spat_config = 0:4) |>
        map(\(x) round(x, 2)) |>
        base::`[`(.pars) |>
        do.call(what = crossing)

    args <- sobol_summs |>
        filter(combo == .combo) |>
        filter(alate_dens == 1, n_pseudo > 0) |>
        select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
        mutate(n_sims = 1000, summ = "all") |>
        as.list()

    sim_list <- 1:nrow(par_perms) |>
        map(\(i) {
            # i = 1L
            # rm(i)
            .args <- args
            for (p in .pars) .args[[p]] <- par_perms[[p]][[i]]
            .sims <- do.call(one_combo, .args) |>
                getElement("outbreak_size") |>
                mean()
            np <- .args[["n_pseudo"]]
            .args[["n_pseudo"]] <- 0L
            .sims0 <- do.call(one_combo, .args) |>
                getElement("outbreak_size") |>
                mean()
            return(tibble(n_pseudo = factor(c(0L, np)),
                          outbreak_size = c(.sims0, .sims)))

        }, .progress = progress_)


    out <- par_perms |>
        mutate(sims = sim_list) |>
        unnest(sims) |>
        mutate(combo = .combo) |>
        select(combo, everything())

    return(out)
}



one_combo_par_manip_plotter <- function(sim_list,
                                        .yvar = "outbreak_size",
                                        .title = waiver(),
                                        .tag = waiver()) {

    .combo <- sim_list[[1]][["combo"]][[1]]

    plot_list <- 1:length(sim_list) |>
        map(\(i) {
            # x_name = "K"
            # rm(x_name, df_og, df_x, df_max_diff)
            tag__ <- .tag
            if (i > 1) tag__ <- waiver()
            df_x <- sim_list[[i]]
            x_name <- colnames(df_x)[2]
            df_og <- sobol_summs |>
                filter(combo == .combo) |>
                filter(alate_dens == 1) |>
                rename(y = !!.yvar) |>
                select(n_pseudo, all_of(x_name), y) |>
                mutate(n_pseudo = factor(n_pseudo))
            df_max_diff <- df_x |>
                rename(y = !!.yvar) |>
                group_by(across(all_of(x_name))) |>
                summarize(dos = y[as.integer(n_pseudo) == 2L] -
                              y[as.integer(n_pseudo) == 1L],
                          .groups = "drop") |>
                filter(abs(dos) == max(abs(dos))) |>
                slice(1) |>
                getElement(x_name)
            p <- df_x |>
                rename(y = !!.yvar) |>
                ggplot(aes(.data[[x_name]], y, color = n_pseudo))
            if (.yvar == "outbreak_size") {
                p <- p +
                    geom_hline(yintercept = c(1, 9), color = "gray70") +
                    scale_y_continuous(breaks = c(1, 5, 9)) +
                    coord_cartesian(ylim = c(1, 9))
            }
            p <- p +
                geom_vline(xintercept = df_max_diff, color = "gray70",
                           linetype = "22") +
                geom_point() +
                geom_line() +
                geom_point(data = df_og, size = 3, shape = 8) +
                labs(x = pretty_params(x_name) |> first_cap(),
                     y = yvar_desc[[.yvar]] |> first_cap(),
                     tag = tag__) +
                scale_color_manual(pretty_params("n_pseudo", TRUE),
                                   values = np_pal) +
                guides(color = guide_legend(override.aes = list(shape = 19))) +
                theme(axis.title.x = element_markdown(),
                      axis.title.y = element_markdown(),
                      legend.title = element_markdown(),
                      plot.tag = element_markdown(size = 16),
                      plot.tag.location = "margin")
            return(p)
        })

    plot_list |>
        do.call(what = wrap_plots) +
        plot_layout(guides = "collect", axes = "collect") +
        plot_annotation(title = .title)
}


one_combo_2par_manip_plotter <- function(x, .contour = FALSE) {

    # x = combo_2par_sims[[1]]
    # rm(x, xvar, yvar, x0, y0, pf, p1, p2)

    xvar <- colnames(x)[colnames(x) != "combo"][1]
    yvar <- colnames(x)[colnames(x) != "combo"][2]

    x0 <- sobol_summs |>
        filter(combo == x$combo[[1]], alate_dens == 1, n_pseudo > 0) |>
        getElement(xvar)
    y0 <- sobol_summs |>
        filter(combo == x$combo[[1]], alate_dens == 1, n_pseudo > 0) |>
        getElement(yvar)

    # function to add shared plot parts:
    pf <- function(d, contour) {
        pp <- d |>
            ggplot(aes(.data[[xvar]], .data[[yvar]])) +
            geom_raster(aes(fill = outbreak_size))
        if (contour) {
            pp <- pp +
                geom_contour(aes(z = outbreak_size), color = "white")
        }
        pp <- pp +
            geom_hline(yintercept = y0, linetype = "22", color = "white",
                       linewidth = 1) +
            geom_vline(xintercept = x0, linetype = "22", color = "white",
                       linewidth = 1) +
            labs(x = pretty_params(xvar) |> first_cap(),
                 y = pretty_params(yvar) |> first_cap()) +
            theme(axis.title.y = element_markdown(),
                  axis.title.x = element_markdown(),
                  strip.text = element_markdown(),
                  legend.title = element_markdown())
        return(pp)
    }
    # Separate by n_pseudo:
    p1 <- x |>
        mutate(n_pseudo = factor(paste(n_pseudo), levels = levels(n_pseudo),
                                 labels = sprintf("n<sub>P</sub> = %s",
                                                  levels(n_pseudo)))) |>
        pf(contour = .contour) +
        facet_wrap(~ n_pseudo, nrow = 1) +
        scale_fill_viridis_c("Outbreak<br>size")
    # Effect of n_pseudo:
    p2 <- x |>
        group_by(across(all_of(c(xvar, yvar))))  |>
        summarize(outbreak_size = outbreak_size[n_pseudo != "0"] -
                      outbreak_size[n_pseudo == "0"],
                  .groups = "drop") |>
        mutate(outbreak_size = round(outbreak_size, 3)) |>
        pf(contour = .contour) +
        scale_fill_viridis_c("Effect of<br>*Pseudomonas* on<br>outbreak size",
                             option = "plasma")

    p1 + p2 +
        plot_layout(nrow = 1, widths = c(2, 1), axes = "collect")

}






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

