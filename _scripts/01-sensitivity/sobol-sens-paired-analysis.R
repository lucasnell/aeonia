
#'
#' Small-scale sensitivity via Sobol indices
#'
#'


source("_scripts/01-sensitivity/00-preamble.R")
source("_scripts/01-sensitivity/sobol-preamble.R")


#' Directory containing output from `sobol-sens-paired.R`:
sobol_dir <- "~/_globus"


#' Output from `sobol-sens-paired.R`:
sobol_sims <- paste0(sobol_dir, "/sobol-sims-paired.rds") |>
    read_rds()




# Summarize each set of simulations:
# Takes ~13 sec  (multithreading doesn't help)
sobol_summs <- map(sobol_sims, \(ss) {
    ss |>
        imap(\(sim_set, i) {
            # Number of pseudo patches in sim with them:
            np <- sim_set$pseudo[["n_pseudo"]][[1]]
            out_df <- sim_set[[1]][1:2, c(names(struct_pars), names(vary_pars))]
            out_df[["combo"]] <- factor(i, levels = 1:length(ss))
            out_df[["n_pseudo"]] <- c(0L, np)
            for (y in yvars) {
                y_p <- sim_set$pseudo[[y]]
                y_np <- sim_set$no_pseudo[[y]]
                if (y != "infect_time" && any(is.na(c(y_p, y_np)))) {
                    stop(y, " has NA values")
                }
                out_df[[y]] <- c(mean(y_np, na.rm = TRUE),
                                 mean(y_p, na.rm = TRUE))
            }
            y_p <- sim_set$pseudo[["outbreak_size"]]
            y_np <- sim_set$no_pseudo[["outbreak_size"]]
            # now do prob. outbreak happened:
            out_df[["p_outbreak"]] <-c(mean(y_np > 1), mean(y_p > 1))
            # and outbreak size when there was one:
            out_df[["outbreak_size2"]] <-c(mean(y_np[y_np > 1]), mean(y_p[y_p > 1]))
            # Lastly, SD(outbreak size):
            out_df[["sd_outbreak_size"]] <-c(sd(y_np), sd(y_p))
            return(out_df)

        }) |>
        list_rbind()
})

# Same thing but looking at differences between with and without Pseudo:
diff_sobol_summs <- map(sobol_sims, \(ss) {
    ss |>
        map(\(sim_set) {
            out_df <- sim_set[[1]][1, c(names(struct_pars), names(vary_pars))]
            out_df[["n_pseudo"]] <- NULL
            for (y in yvars) {
                y_p <- sim_set$pseudo[[y]]
                y_np <- sim_set$no_pseudo[[y]]
                if (y != "infect_time" && any(is.na(c(y_p, y_np)))) {
                    stop(y, " has NA values")
                }
                # using rounding bc weird, very small values show up here (~1e-16)
                # if I don't. Plus, I know that the difference can't have more than
                # 2 decimal digits bc n_sims = 100.
                out_df[[y]] <- round(mean(y_p, na.rm = TRUE) -
                                         mean(y_np, na.rm = TRUE), 2)
            }
            # now do prob. outbreak happened:
            out_df[["p_outbreak"]] <- mean(y_p > 1) - mean(y_np > 1)
            # and outbreak size when there was one:
            out_df[["outbreak_size2"]] <- mean(y_p[y_p > 1]) - mean(y_np[y_np > 1])

            return(out_df)

        }) |>
        list_rbind()
})







scatter <- function(sim_outs, yvar = "outbreak_size", .facet_nrow = NULL, .title = NULL) {

    # sim_outs = sobol_summs[[2]]; yvar = "sd_outbreak_size"; .facet_nrow = NULL; .title = NULL
    # rm(sim_outs, yvar, .facet_nrow, .title, .vary_pars, .df, has_n_pseudo, np_name, .ylab, .y_range, .y_breaks, plot_list)

    # if (yvar == "infect_time") sim_outs <- sim_outs |>
    #     mutate(infect_time = ifelse(is.na(infect_time), 101, infect_time))


    if (is.null(.title)) {
        .title <- colnames(sim_outs) |>
            keep(\(x) x %in% names(struct_pars) & x != "n_pseudo") |>
            map(\(x) sprintf("%s = %.3g", pretty_params(x, TRUE),
                             sim_outs[[x]][[1]])) |>
            c(list(collapse = ";")) |>
            do.call(what = paste)
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
                     xvars = c("log_aphids", "log_alates"),
                     yvar = "outbreak_size",
                     .title = NULL) {

    if (is.null(.title)) {
        .title <- colnames(sim_outs) |>
            keep(\(x) x %in% names(struct_pars) & x != "n_pseudo") |>
            map(\(x) sprintf("%s = %.3g", pretty_params(x, TRUE),
                             sim_outs[[x]][[1]])) |>
            c(list(collapse = ";")) |>
            do.call(what = paste)
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


sum(diff_sobol_summs[[2L]][["outbreak_size"]] > 0)
sum(diff_sobol_summs[[2L]][["outbreak_size"]] == 0)
sum(diff_sobol_summs[[2L]][["outbreak_size"]] < 0)




scatter(sobol_summs[[2L]], "outbreak_size")
scatter(sobol_summs[[2L]], "sd_outbreak_size")

sobol_summs[[2L]] |>
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
    theme(axis.title = element_markdown(),
          legend.title = element_markdown())




# LEFT OFF #1 ----
#' For parameter combos that seemed to result in a negative
#' effect of Pseudomonas, do they continue to show this pattern with greater
#' numbers of simulations?
# Test the top 100 parameter combinations:
diff_test_df <- diff_sobol_summs[[2L]] |>
    arrange(desc(outbreak_size)) |>
    select(all_of(names(vary_pars)), outbreak_size) |>
    slice_head(n = 100)

#' Output from `sobol-sens-paired-test.R`:
sobol_test_sims <- paste0(sobol_dir, "/sobol-test-sims-paired.rds") |>
    read_rds()


sobol_test_df <- imap(sobol_test_sims,
                     \(sim_set,  i) {
                         out_df <- sim_set[[1]][1:2, names(vary_pars)]
                         out_df[["n_pseudo"]] <- c(3L, 0L)
                         y_p <- sim_set$pseudo[["outbreak_size"]]
                         y_np <- sim_set$no_pseudo[["outbreak_size"]]
                         if (any(is.na(c(y_p, y_np)))) {
                             stop("outbreak_size has NA values")
                         }
                         out_df[["outbreak_size"]] <- round(c(mean(y_p),
                                                              mean(y_np)), 3)
                         out_df[["sd_outbreak_size"]] <- round(c(sd(y_p),
                                                                 sd(y_np)), 3)
                         # now do prob. outbreak happened:
                         out_df[["p_outbreak"]] <- round(c(mean(y_p > 1),
                                                           mean(y_np > 1)), 3)
                         out_df[["combo"]] <- factor(i, levels = 1:nrow(diff_test_df))
                         return(out_df)
                     }) |>
    list_rbind() |>
    select(combo, everything())

# How do these compare to previous simulations?
sobol_test_df |>
    group_by(combo) |>
    summarize(outbreak_size = outbreak_size[n_pseudo > 0L] -
                  outbreak_size[n_pseudo == 0L]) |>
    mutate(outbreak_size0 = diff_test_df$outbreak_size) |>
    ggplot(aes(outbreak_size, outbreak_size0)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "22") +
    geom_point() +
    coord_equal(xlim = c(0, 3), ylim = c(0, 3))
# Pretty good!



sobol_test_df |>
    mutate(n_pseudo = factor(n_pseudo)) |>
    ggplot(aes(outbreak_size, sd_outbreak_size, color = n_pseudo)) +
    geom_point() +
    scale_color_manual(pretty_params("n_pseudo", TRUE),
                       values = c("goldenrod", "dodgerblue"))

sobol_test_df |>
    group_by(combo) |>
    summarize(outbreak_size = outbreak_size[n_pseudo > 0L] -
                  outbreak_size[n_pseudo == 0L],
              sd_outbreak_size = sd_outbreak_size[n_pseudo > 0L] -
                  sd_outbreak_size[n_pseudo == 0L]) |>
    ggplot(aes(outbreak_size, sd_outbreak_size)) +
    geom_point()


# For these top 100 combos that result in Pseudomonas being bad for plants,
# is the aphid density ~ alate production feedback necessary?

# Takes ~20 sec
set.seed(743235938)
sobol_test0_sims <- diff_test_df |>
    select(names(vary_pars)) |>
    pmap(\(Y0, mean_N, sd_N, K, virus_attract, pseudo_repel, pseudo_surv,
           zeta, spat_config) {

        args <- list(Y0 = Y0, mean_N = mean_N, sd_N = sd_N, K = K,
                     virus_attract = virus_attract, pseudo_repel = pseudo_repel,
                     pseudo_surv = pseudo_surv, zeta = zeta,
                     spat_config = spat_config,
                     n_sims = 1000L,
                     alate_slope = 0, alate_max = 0.05, n_pseudo = 3L)

        sim <- do.call(one_combo, args)

        args[["n_pseudo"]] <- 0L
        sim0 <- do.call(one_combo, args)

        out <- list(pseudo = sim, no_pseudo = sim0)
    }, .progress = .prog_args)


sobol_test0_df <- imap(sobol_test0_sims,
     \(sim_set, i) {
         out_df <- sim_set[[1]][1:2, names(vary_pars)]
         out_df[["n_pseudo"]] <- c(sim_set$pseudo[["n_pseudo"]][[1]], 0L)
         y_p <- sim_set$pseudo[["outbreak_size"]]
         y_np <- sim_set$no_pseudo[["outbreak_size"]]
         if (any(is.na(c(y_p, y_np)))) {
             stop("outbreak_size has NA values")
         }
         out_df[["outbreak_size"]] <- round(c(mean(y_p),
                                              mean(y_np)), 3)
         out_df[["sd_outbreak_size"]] <- round(c(sd(y_p),
                                                 sd(y_np)), 3)
         # now do prob. outbreak happened:
         out_df[["p_outbreak"]] <- round(c(mean(y_p > 1),
                                           mean(y_np > 1)), 3)
         out_df[["combo"]] <- factor(i, levels = 1:nrow(diff_test_df))
         return(out_df)
     }) |>
    list_rbind() |>
    select(combo, everything())




bind_cols(sobol_test_df |>
              select(combo, n_pseudo, outbreak_size),
          sobol_test0_df |>
              select(outbreak_size) |>
              rename(outbreak_size0 = outbreak_size)) |>
    ggplot(aes(outbreak_size, outbreak_size0)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "22") +
    geom_point() +
    coord_equal() +
    labs(x = "Outbreak size with alates ~ density",
         y = "Outbreak size constant alates")


bind_cols(sobol_test_df |>
              select(combo, n_pseudo, outbreak_size),
          sobol_test0_df |>
              select(outbreak_size) |>
              rename(outbreak_size0 = outbreak_size)) |>
    group_by(combo) |>
    summarize(outbreak_size = outbreak_size[n_pseudo > 0L] -
                  outbreak_size[n_pseudo == 0L],
              outbreak_size0 = outbreak_size0[n_pseudo > 0L] -
                  outbreak_size0[n_pseudo == 0L]) |>
    ggplot(aes(outbreak_size, outbreak_size0)) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "22") +
    geom_point() +
    coord_equal() +
    labs(x = "Effect of *Pseudomonas* on outbreak size with alates ~ density",
         y = "Effect of *Pseudomonas* on outbreak size constant alates") +
    theme(axis.title.x = element_markdown(),
          axis.title.y = element_markdown())









# LEFT OFF #2 ----
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







# =============================================================================*
# =============================================================================*
# Random Forest approach ----
# =============================================================================*
# =============================================================================*


library(randomForest)
library(iml)

# Dataset for outbreak size:
set.seed(99240974) # used for shuffling dataset
ob_sim_df <- sobol_summs[[2]] |>
    filter(n_pseudo > 0) |>
    select(outbreak_size, all_of(c(names(vary_pars)))) |>
    slice_sample(prop = 1)

# Dataset for effect of pseudo on outbreak size:
set.seed(1976042860) # used for shuffling dataset
dob_sim_df <- sobol_sims[[2]] |>
    map(\(sim_set) {
        out_df <- sim_set[[1]][1, names(vary_pars)]
        y_p <- sim_set$pseudo[["outbreak_size"]]
        y_np <- sim_set$no_pseudo[["outbreak_size"]]
        out_df[["outbreak_size"]] <- mean(y_p, na.rm = TRUE) - mean(y_np, na.rm = TRUE)
        return(out_df)
    }) |>
    list_rbind() |>
    select(outbreak_size, everything()) |>
    slice_sample(prop = 1)

n_train <- as.integer(round(nrow(ob_sim_df) * 0.75))
n_test <- nrow(ob_sim_df) - n_train

ob_train_df <- ob_sim_df[1:n_train,]
dob_train_df <- dob_sim_df[1:n_train,]
ob_test_df <- ob_sim_df[(n_train+1):nrow(ob_sim_df),]
dob_test_df <- dob_sim_df[(n_train+1):nrow(dob_sim_df),]





# # Which value of `mtry` results in most explained variance?
# # Each step takes ~3 min
# set.seed(45181283)
# ob_m_df <- tibble(m = 2:length(vary_pars)) |>
#     mutate(v = map_dbl(m, \(.m) {
#         rf <- randomForest(ob_sim_df[1:n_train,names(vary_pars)],
#                            ob_sim_df[["outbreak_size"]][1:n_train],
#                            corr.bias = TRUE, mtry = .m)
#         return(tail(rf$rsq, 1))
#     }))
# set.seed(1547476225)
# dob_m_df <- tibble(m = 2:length(vary_pars)) |>
#     mutate(v = map_dbl(m, \(.m) {
#         rf <- randomForest(dob_sim_df[1:n_train,names(vary_pars)],
#                            dob_sim_df[["outbreak_size"]][1:n_train],
#                            corr.bias = TRUE, mtry = .m)
#         return(tail(rf$rsq, 1))
#     }))
# ob_m_df |>
#     ggplot(aes(m, v)) +
#     geom_line() +
#     geom_point() +
#     scale_x_continuous(breaks = 2:10) +
#     theme(panel.grid.major = element_line(color = "gray80"))
# dob_m_df |>
#     ggplot(aes(m, v)) +
#     geom_line() +
#     geom_point() +
#     scale_x_continuous(breaks = 2:10) +
#     theme(panel.grid.major = element_line(color = "gray80"))
# # The best for both was mtry = 5



if (!file.exists("_scripts/interm-data/randomforest-outbreak_size.rds") |
    !file.exists("_scripts/interm-data/randomforest-diff_outbreak_size.rds")) {
    set.seed(359013909)
    ob_rf <- randomForest(ob_train_df[,names(vary_pars)],
                          ob_train_df[["outbreak_size"]],
                          importance = TRUE, corr.bias = TRUE, mtry = 5)
    set.seed(84490683)
    dob_rf <- randomForest(dob_train_df[,names(vary_pars)],
                           dob_train_df[["outbreak_size"]],
                           importance = TRUE, corr.bias = TRUE, mtry = 5)
    write_rds(ob_rf, "_scripts/interm-data/randomforest-outbreak_size.rds", compress = "gz")
    write_rds(dob_rf, "_scripts/interm-data/randomforest-diff_outbreak_size.rds", compress = "gz")
} else {
    ob_rf <- read_rds("_scripts/interm-data/randomforest-outbreak_size.rds")
    dob_rf <- read_rds("_scripts/interm-data/randomforest-diff_outbreak_size.rds")
}



# variance explained:
(ob_rf_rsq <- tail(ob_rf$rsq, 1))
# [1] 0.9556963

# variance explained:
(dob_rf_rsq <- tail(dob_rf$rsq, 1))
# [1] 0.8015296



# randomForest:::partialPlot.randomForest



# Calculate data for partial dependence plots (takes a while to run, so this
# should be done separately from plots)
partial_calc <- function(x, pred_data, .prog_args = FALSE) {

    # x = ob_rf; pred_data = ob_train_df
    # rm(x, pred_data, x_names, n)

    n <- nrow(pred_data)
    x$forest$xlevels |>
        names() |>
        map(\(x_name) {
            # x_name = x_names[[1]]
            # rm(x_name, xv, n_pt, x_pt, y_pt, x_data, i)
            xv <- pred_data[[x_name]]
            n_pt <- min(length(unique(xv)), 51)
            x_pt <- seq(min(xv), max(xv), length = n_pt)
            y_pt <- numeric(n_pt)
            x_data <- pred_data
            for (i in 1:n_pt) {
                x_data[, x_name] <- rep(x_pt[i], n)
                y_pt[i] <- mean(predict(x, x_data), na.rm = TRUE)
            }
            return(tibble(param = x_name, x = x_pt, y = y_pt))
        }, .progress = .prog_args) |>
        list_rbind()

}


if (!file.exists("_scripts/interm-data/randomforestParts-outbreak_size.rds") |
    !file.exists("_scripts/interm-data/randomforestParts-diff_outbreak_size.rds")) {
    # Each step takes ~ 5 min
    ob_rf_parts <- partial_calc(ob_rf, ob_train_df, .prog_args)
    dob_rf_parts <- partial_calc(dob_rf, dob_train_df, .prog_args)
    write_rds(ob_rf_parts, "_scripts/interm-data/randomforestParts-outbreak_size.rds", compress = "gz")
    write_rds(dob_rf_parts, "_scripts/interm-data/randomforestParts-diff_outbreak_size.rds", compress = "gz")
} else {
    ob_rf_parts <- read_rds("_scripts/interm-data/randomforestParts-outbreak_size.rds")
    dob_rf_parts <- read_rds("_scripts/interm-data/randomforestParts-diff_outbreak_size.rds")
}





partial_plot <- function(x, .ylab) {
    # x = ob_rf_parts; .ylab = "Outbreak size"
    # rm(x, .ylab)
    x |>
        mutate(param = pretty_params(param) |>
                   factor(levels = pretty_params(names(vary_pars)))) |>
        ggplot(aes(x, y)) +
        geom_line() +
        geom_point() +
        facet_wrap(~ param, scales = "free_x", strip.position = "bottom") +
        labs(y = .ylab) +
        theme(plot.title = element_markdown(),
              axis.title.y = element_markdown(),
              strip.text = element_markdown(),
              strip.placement = "outside",
              axis.title.x = element_blank())
}




# Plot predicted vs observed:
pred_plot <- function(x, test_x, test_y, .title = "") {
    # x = ob_rf
    # test_x = ob_test_df[,names(vary_pars)]
    # test_y = ob_test_df[["outbreak_size"]]
    # .title = ""
    # rm(x, test_x, test_y, .title, pred_rf, rf_r2)
    stopifnot(inherits(x, "randomForest"))
    stopifnot(inherits(test_x, "data.frame"))
    stopifnot(inherits(test_y, "numeric") || inherits(test_y, "factor"))
    stopifnot(length(test_y) == nrow(test_y))
    stopifnot(inherits(.title, "character") && length(.title) == 1L)
    pred_rf <- predict(x, test_x)
    rf_r2 <- sprintf("r^2 == %.3f", cor(pred_rf, test_y)^2)
    tibble(Predicted = pred_rf, Observed = test_y) |>
        ggplot(aes(Observed, Predicted)) +
        geom_point(shape = 1) +
        ggtitle(.title) +
        geom_text(data = tibble(lab = rf_r2, x = min(pred_rf), y = max(pred_rf)),
                  aes(x, y, label = lab), hjust = 0, vjust = 1, parse = TRUE) +
        geom_abline(slope = 1, intercept = 0, linetype = 2, color = "red") +
        coord_equal() +
        theme(plot.title = element_markdown())
}

# Plot variable importance:
imp_plot <- function(x, .title = "") {
    # x = ob_rf; .title = ""
    # rm(x, .title)
    stopifnot(inherits(x, "randomForest"))
    stopifnot(inherits(.title, "character") && length(.title) == 1L)
    stopifnot("%IncMSE" %in% colnames(x$importance))
    importance(x) |>
        as.data.frame() |>
        rownames_to_column("par") |>
        select(par, `%IncMSE`) |>
        rename(inc_mse = `%IncMSE`) |>
        mutate(par = pretty_params(par),
               par = fct_reorder(par, inc_mse, .na_rm = TRUE)) |>
        ggplot(aes(inc_mse, par)) +
        geom_vline(xintercept = 0, linewidth = 1, color = "gray70") +
        geom_point(size = 3) +
        geom_segment(aes(xend = 0, yend = par), linewidth = 1) +
        ggtitle(.title) +
        xlab("Mean increase in MSE") +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_markdown(size = 11, color = "black"),
              plot.title = element_markdown())
}


# Plot variable interaction strength:
inter_plot <- function(x, .title = "") {
    # x = ob_rf_int; .title = ""
    # rm(x, .title, x_df)
    stopifnot(inherits(x, "Interaction"))
    stopifnot(inherits(.title, "character") && length(.title) == 1L)

    x_df <- x[["results"]] |> as_tibble()

    p <- x_df |>
        filter(!is.na(.interaction)) |>
        mutate(.feature = pretty_params(.feature)) |>
        mutate(.feature = fct_reorder(.feature, .interaction)) |>
        ggplot(aes(.interaction, .feature)) +
        geom_vline(xintercept = 0, linewidth = 1, color = "gray70") +
        geom_point(size = 3) +
        geom_segment(aes(xend = 0, yend = .feature), linewidth = 1) +
        ggtitle(.title) +
        labs(x = "Overall interaction strength") +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_markdown(size = 11, color = "black"),
              plot.title = element_markdown())
    if (".class" %in% colnames(x_df)) {
        p <- p + facet_wrap(~ .class)
    }
    return(p)
}


imp_inter_plot <- function(rf, rf_int, .title = "") {

    # rf = ob_rf; rf_int = ob_rf_int; .title = ""
    # rm(rf, rf_int, .title, .df, mse_max, int_max, trans, itrans, .mse_col, .int_col)
    stopifnot(inherits(rf, "randomForest"))
    stopifnot("%IncMSE" %in% colnames(rf$importance))
    stopifnot(inherits(rf_int, "Interaction"))
    stopifnot(inherits(.title, "character") && length(.title) == 1L)

    .df <- importance(rf) |>
        as.data.frame() |>
        rownames_to_column("par") |>
        select(par, `%IncMSE`) |>
        rename(inc_mse = `%IncMSE`) |>
        left_join(rf_int[["results"]] |> as_tibble() |> rename(par = .feature, inter = .interaction),
                  by = "par") |>
        pivot_longer(-par, names_to = "measure")

    mse_max <- max(.df$value[.df$measure == "inc_mse"])
    int_max <- max(.df$value[.df$measure == "inter"])
    trans <- \(x) x * mse_max / int_max
    itrans <- \(x) x * int_max / mse_max

    .mse_col <- "dodgerblue"
    .int_col <- "goldenrod"

    .df |>
        mutate(par = pretty_params(par) |>
                   factor(levels = rev(pretty_params(names(vary_pars)))),
               value = ifelse(measure == "inter", trans(value), value)) |>
        ggplot(aes(value, par, fill = measure)) +
        geom_vline(xintercept = 0, linewidth = 1, color = "gray70") +
        geom_col(position = position_dodge(0.5), color = NA, width = 0.4) +
        ggtitle(.title) +
        scale_x_continuous("Mean increase in MSE",
                           sec.axis = sec_axis(itrans, "Overall interaction strength")) +
        scale_fill_manual(values = c(inc_mse = .mse_col, inter = .int_col), guide = "none") +
        theme(axis.title.y = element_blank(),
              axis.text.y = element_markdown(size = 11, color = "black"),
              axis.title.x.top = element_markdown(color = .int_col, face = "bold"),
              axis.text.x.top = element_markdown(color = .int_col),
              axis.ticks.x.top = element_line(color = .int_col),
              axis.title.x.bottom = element_markdown(color = .mse_col, face = "bold"),
              axis.text.x.bottom = element_markdown(color = .mse_col),
              axis.ticks.x.bottom = element_line(color = .mse_col),
              plot.title = element_markdown())

}




#' Create `iml` Predictor objects.
ob_rf_pred <- Predictor$new(ob_rf, data = ob_train_df[,names(vary_pars)],
                            y = ob_train_df[["outbreak_size"]])
dob_rf_pred <- Predictor$new(dob_rf, data = dob_train_df[,names(vary_pars)],
                             y = dob_train_df[["outbreak_size"]])
#' Create `iml` Interaction objects.
#' This function creates an Interaction option and temporarily allows
#' the future package to access large objects:
make_inter <- \(rf_obj) {
    oopts <- options(future.globals.maxSize = 1.5e9)  ## 1.5 GB
    on.exit(options(oopts))
    Interaction$new(rf_obj)
}
if (!file.exists("_scripts/interm-data/randomforestInt-outbreak_size.rds") |
    !file.exists("_scripts/interm-data/randomforestInt-diff_outbreak_size.rds")) {
    # Make futures garbage collect to save memory:
    plan(multisession, workers = options()[["mc.cores"]], gc = TRUE)
    # Takes ~1.5 min each
    ob_rf_int <- make_inter(ob_rf_pred)
    dob_rf_int <- make_inter(dob_rf_pred)
    write_rds(ob_rf_int, "_scripts/interm-data/randomforestInt-outbreak_size.rds", compress = "gz")
    write_rds(dob_rf_int, "_scripts/interm-data/randomforestInt-diff_outbreak_size.rds", compress = "gz")
    # Go back to default configuration:
    plan(multisession, workers = options()[["mc.cores"]])
} else {
    ob_rf_int <- read_rds("_scripts/interm-data/randomforestInt-outbreak_size.rds")
    dob_rf_int <- read_rds("_scripts/interm-data/randomforestInt-diff_outbreak_size.rds")
}





(partial_plot(ob_rf_parts, "Outbreak size") |
        partial_plot(dob_rf_parts, "Effect of *Pseudomonas* on outbreak size"))






rf_pred_p <- (pred_plot(ob_rf, ob_test_df[,names(vary_pars)],
           ob_test_df[["outbreak_size"]],
           "Outbreak size") |
        pred_plot(dob_rf, dob_test_df[,names(vary_pars)],
                  dob_test_df[["outbreak_size"]],
                  "Effect of *Pseudomonas* on outbreak size")) +
    plot_layout(nrow = 1) +
    plot_annotation(title = "Observed vs predicted",
                    theme = theme(plot.title = element_markdown(size = 18)))

save_plot("_plots/rf-pred.pdf", rf_pred_p, width = 8, height = 5)



# (imp_plot(ob_rf, "Outbreak size") | imp_plot(dob_rf, "Effect of *Pseudomonas* on outbreak size")) +
#     plot_layout(nrow = 1) +
#     plot_annotation(title = "Variable importance",
#                     theme = theme(plot.title = element_markdown(size = 18)))
#
# (inter_plot(ob_rf_int, "Outbreak size") | inter_plot(dob_rf_int, "Effect of *Pseudomonas* on outbreak size")) +
#     plot_layout(nrow = 1) +
#     plot_annotation(title = "Variable interaction strength",
#                     theme = theme(plot.title = element_markdown(size = 18)))


rf_imp_inter_p <- (imp_inter_plot(ob_rf, ob_rf_int, .title = "Outbreak size") |
        imp_inter_plot(dob_rf, dob_rf_int, .title = "Effect of *Pseudomonas* on<br>outbreak size")) +
    plot_layout(nrow = 1) +
    plot_annotation(title = "Variable importance and interaction strength",
                    theme = theme(plot.title = element_markdown(size = 18)))

save_plot("_plots/rf-imp_inter.pdf", rf_imp_inter_p, width = 10, height = 5)





