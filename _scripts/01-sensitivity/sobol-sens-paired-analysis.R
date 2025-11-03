
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
                    stat_summary(fun = "mean", size = 4, geom = "point") +
                    theme(legend.position = "none")
            } else {
                plt <- plt +
                    stat_smooth(method = "gam",
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




# LEFT OFF #1 ----
#' For parameter combos that seemed to result in a negative
#' effect of Pseudomonas, does alate ~ density affect outcomes?
# Test the top 100 parameter combinations:
diff_test_combos <- diff_sobol_summs |>
    filter(alate_dens == 1) |>
    arrange(desc(outbreak_size)) |>
    getElement("combo") |>
    head(n = 100)




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

# save_plot("_plots/alate-dens.pdf", alate_dens_p, width = 6, height = 6)





# LEFT OFF #2 ----
#' What parameter values are associated with Pseudomonas being bad for plants?

# Paired time series for with and without Pseudomonas
# Note: only takes output summarized by plant, will summarize in the function
#       if `by_plant  = FALSE`
timeseries <- function(sim_df_p,
                       sim_df_np,
                       by_plant = FALSE,
                       .rep = NULL,
                       .title = waiver()) {

    # sim_df_p = sim_df; sim_df_np = sim0_df; .rep = NULL; by_plant = TRUE; .title = waiver()
    # rm(sim_df_p, sim_df_np, by_plant, .rep, .title, all_sims, aphids, wasps, virus, max_virus, max_aphids, max_wasps, itrans, trans, itrans2, trans2, aphid_breaks, aphid_labels, aphid_ylab, col_pal, p)

    stopifnot(identical(colnames(sim_df_p), colnames(sim_df_np)))
    stopifnot("n_pseudo" %in% colnames(sim_df_p))
    stopifnot("x" %in% colnames(sim_df_p) && "y" %in% colnames(sim_df_p))
    if (is.null(.rep) && !by_plant) .rep <- unique(sim_df_p$rep)
    if (is.null(.rep) && by_plant) .rep <- sample(unique(sim_df_p$rep), 1)
    stopifnot(all(.rep %in% sim_df_p$rep) && all(.rep %in% sim_df_np$rep))
    stopifnot(!by_plant || length(.rep) == 1)

    all_sims <- bind_rows(sim_df_p, sim_df_np) |>
        filter(rep %in% .rep) |>
        mutate(x = ifelse(x == 0, NA, x),
               y = ifelse(y == 0, NA, y),
               plant = interaction(x, y),
               rep = factor(rep),
               n_pseudo = factor(n_pseudo),
               id = interaction(n_pseudo, rep)) |>
        mutate(aphids = aphids + parasitized) |>
        select(id, n_pseudo, rep, plant, time, virus, aphids, alates, wasps)
    wasps <- all_sims |>
        filter(!is.na(wasps)) |>
        mutate(id = interaction(n_pseudo, rep)) |>
        select(id, n_pseudo, rep, time, wasps)
    aphids <- all_sims |>
        filter(!is.na(plant)) |>
        select(id, n_pseudo:time, aphids)
    if (by_plant) {
        virus <- all_sims |>
            filter(!is.na(plant)) |>
            group_by(n_pseudo, rep, plant) |>
            filter(virus == 1) |>
            filter(time == min(time)) |>
            ungroup() |>
            select(id:time, virus)
        max_virus <- 1
        alates <- all_sims |>
            filter(!is.na(plant)) |>
            select(id, n_pseudo:time, alates)
    } else {
        virus <- all_sims |>
            filter(is.na(plant)) |>
            select(id, n_pseudo, rep, time, virus)
        # Max virus is the total number of plants:
        max_virus <- bind_rows(sim_df_p, sim_df_np) |>
            filter(x > 0) |>
            distinct(x, y) |>
            nrow()
        aphids <- aphids |>
            group_by(id, n_pseudo, rep, time) |>
            summarize(aphids = sum(aphids), .groups = "drop")
        alates <- all_sims |>
            filter(is.na(plant)) |>
            select(id, n_pseudo, rep, time, alates)
    }

    max_aphids <- max(aphids$aphids, na.rm = TRUE)
    max_wasps <- max(wasps$wasps, na.rm = TRUE)
    max_alates <- max(alates$alates, na.rm = TRUE)
    # convert from wasps --> aphids:
    trans <- \(x) x * max_aphids / max_wasps
    # convert from aphids --> wasps:
    itrans <- \(x) x * max_wasps / max_aphids
    # convert from alates --> aphids:
    trans2 <- \(x) x * max_aphids / max_alates
    # convert from aphids --> alates:
    itrans2 <- \(x) x * max_alates / max_aphids

    col_pal <- viridisLite::plasma(4) |>
        set_names(c("aphids", "alates", "wasps", "virus")) |>
        as.list()

    aphid_breaks <- scales::breaks_extended(n = 4)(aphids$aphids)
    aphid_labels <- sprintf(paste0("%s (<span style=\"color: ", col_pal$alates,
                                   ";\">%.1f</span>)"), aphid_breaks,
                            itrans2(aphid_breaks))
    aphid_ylab <- paste0("Aphid density (<span style=\"color: ", col_pal$alates,
                        ";\">alate density</span>)")


    p <- aphids |>
        ggplot(aes(time, alpha = n_pseudo, linewidth = n_pseudo)) +
        geom_line(aes(y = aphids, group = id),
                  color = col_pal$aphids) +
        geom_line(data = alates, aes(y = alates * max_aphids / max_alates, group = id),
                  color = col_pal$alates) +
        geom_line(data = wasps, aes(y = trans(wasps), group = id),
                  color = col_pal$wasps) +
        scale_y_continuous(aphid_ylab,
                           breaks = aphid_breaks,
                           labels = aphid_labels,
                           sec.axis = sec_axis(itrans, "Wasp density")) +
        # scale_linetype_manual(pretty_params("n_pseudo", TRUE),
        #                       values = c("22", "solid")) +
        scale_alpha_manual(pretty_params("n_pseudo", TRUE), values = c(0.3, 1)) +
        scale_linewidth_manual(pretty_params("n_pseudo", TRUE), values = c(1.25, 0.75)) +
        guides(alpha = guide_legend(override.aes =
                                        list(color = alpha("gray", c(0.3, 1))))) +
        labs(title = .title, x = "Time (days)") +
        theme(plot.title = element_markdown(hjust = 0.5),
              legend.title = element_markdown(),
              axis.title.y.left = element_markdown(color = col_pal$aphids,
                                                   face = "bold"),
              axis.title.y.right = element_markdown(color = col_pal$wasps,
                                                    face = "bold"),
              axis.text.y.left = element_markdown(color = col_pal$aphids),
              axis.text.y.right = element_markdown(color = col_pal$wasps))

    if (by_plant) {
        p <- p +
            # geom_vline(data = virus, aes(xintercept = time), color = "#EC008C") +
            facet_wrap( ~ plant, nrow = 3) +
            scale_color_viridis_d(begin = 0.1, end = 0.9, guide = "none")
    } else {
        p <- p +
            geom_line(data = virus, aes(y = virus * max_aphids / max_virus),
                      color = col_pal$virus) +
            facet_wrap(~ rep, labeller = label_both)
    }

    return(p)

}

calc_metrics <- function(sims) {
    sims |>
        filter(!is.na(wasps)) |>
        group_by(rep) |>
        summarize(outbreak_size = max(virus),
                  log_aphids = mean(log10(1 + aphids + parasitized)),
                  log_alates = mean(log10(1 + alates)),
                  log_wasps = mean(log10(1 + wasps)))
}


diff_sobol_summs |>
    filter(alate_dens == 1) |>
    arrange(desc(outbreak_size)) |>
    select(combo, all_of(names(vary_pars)), outbreak_size)

i = 7112L

set.seed(1032439295)
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

timeseries(sim_df, sim0_df, .title = sprintf("combo %i", i))
timeseries(sim_df, sim0_df, TRUE, 4, .title = sprintf("combo %i", i))



diff_sobol_summs |>
    filter(alate_dens == 1, spat_config == 1) |>
    ggplot(aes(zeta, K)) +
    geom_point(aes(color = outbreak_size)) +
    scale_color_viridis_c()

diff_sobol_summs |>
    filter(alate_dens == 1, spat_config == 1) |>
    ggplot(aes(zeta, outbreak_size)) +
    geom_point(alpha = 0.1) +
    stat_smooth(method = "gam",
                formula = y ~ s(x, bs = "cs"),
                se = TRUE, linewidth = 1)




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






