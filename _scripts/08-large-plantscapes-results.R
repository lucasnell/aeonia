#'
#' Plots for larger landscape simulations
#'

source("_scripts/00-preamble.R")


.overwrite <- FALSE


#' Because below returns an empty tibble, it shows that when we set
#' `wt_vp = 1e-6` it does indeed result in no overlap between
#' *Pseudomonas* and virus.
#' Unsurprisingly, when we set `wt_vp = 100` and manually have them overlap,
#' there is always overlap.
#'
# list.files("_scripts/interm-data", "large-plantscapes-.?.?.rds",
#            full.names = TRUE) |>
#     map(read_rds) |>
#     list_rbind() |>
#     select(n_pseudo:wt_pp, landscape) |>
#     filter(!is.na(wt_pp)) |>
#     # Predicted number of instances of both Pseudomonas and virus
#     # across 100 sims:
#     mutate(n_both_pred = ifelse(wt_vp < 1, 0L, 100L)) |>
#     # Observed:
#     mutate(n_both_obs = map_int(landscape, \(x) sum(x == 3L))) |>
#     select(n_pseudo:wt_pp, starts_with("n_both")) |>
#     filter(n_both_pred != n_both_obs)






# from 06-large-plantscapes.sh:
# Takes ~5 sec
set.seed(2120927824) # for bootstrapping
sim_df <- list.files("_scripts/interm-data", "large-plantscapes-.?.?.rds",
                     full.names = TRUE) |>
    map(\(x) {
        read_rds(x) |>
            # Remove this vector immediately bc it's quite large
            select(-landscape)
    }) |>
    list_rbind() |>
    mutate(outbreak_size = map_dbl(sim, \(x) mean(x$n_infected[x$n_infected > 1])),
           log_outbreak_size = map_dbl(sim, \(x) mean(log10(x$n_infected[x$n_infected > 1]))),
           p_emerge = map_dbl(sim, \(x) mean(x$n_infected > 1)),
           n_infected = map_dbl(sim, \(x) mean(x$n_infected)),
           boots = map(sim, \(x) ci_booter(x$n_infected, "all"))) |>
    mutate(wt_vp = ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") |>
               factor(levels = c("off *Pseudo.*", "on *Pseudo.*")),
           wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
               factor(levels = c("uniform", "clustered")),
           outbreaks = factor(outbreaks, levels = c("small", "big_zh", "big_zl")),
           wasp_resp = factor(wasp_resp, levels = levels(wasp_resp_fct)))




#' Basic check to make sure that the two scenarios still do what we expect
#' them to: Pseudomonas decreases outbreak size in "strong" and increases in "weak"

sim_df |>
    filter(sd_N == 0 & virus_attract == 1 & pseudo_repel == 1 &
               ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform" & n_pseudo == 7e3) |
                    n_pseudo == 0)) |>
    select(wasp_resp, outbreaks, n_pseudo, outbreak_size) |>
    arrange(wasp_resp, outbreaks, n_pseudo)
#    wasp_resp outbreaks n_pseudo outbreak_size
#    <fct>     <fct>        <dbl>         <dbl>
#  1 weak      small            0          3.03
#  2 weak      small         7000          4.40
#  3 weak      big_zh           0       6338.
#  4 weak      big_zh        7000       9337.
#  5 weak      big_zl           0       2370.
#  6 weak      big_zl        7000       7363.
#  7 strong    small            0          2.67
#  8 strong    small         7000          2.31
#  9 strong    big_zh           0       6355.
# 10 strong    big_zh        7000       2412.
# 11 strong    big_zl           0       2220.
# 12 strong    big_zl        7000        606.




#
# Add np = 0 for each landscape type.
# Used to construct input for `baseline_plotter`.
#
add_no_pseudo_points <- function(data_df) {

    stopifnot(all(c("wt_vp", "wt_pp") %in% colnames(data_df)))
    split_cols <- c("wasp_resp", "outbreaks", "sd_N", "virus_attract", "pseudo_repel") |>
        keep(\(x) x %in% colnames(data_df))
    stopifnot(length(split_cols) >= 1)
    split_form <- paste("~", paste(split_cols, collapse = "+")) |> as.formula()

    data_df  |>
        split(split_form, drop = TRUE) |>
        map(\(x) {
            row_n0 <- filter(x, is.na(wt_vp)) |> select(-wt_vp, -wt_pp)
            stopifnot(nrow(row_n0) == 1)
            rows_np <- x |>
                filter(!is.na(wt_vp))
            rows_np |>
                distinct(wt_vp, wt_pp) |>
                mutate(obs = map(1:n(), \(i) row_n0)) |>
                unnest(obs) |>
                bind_rows(rows_np)
        }) |>
        list_rbind()

}



baseline_plotter <- function(filt_df, outcomes = "all",
                             col_fct = NULL,
                             color_vals = NULL,
                             inters = FALSE,
                             fix_y = TRUE,
                             obs_breaks = NULL,
                             incl_vals = FALSE,
                             return_list = FALSE,
                             .title = list(waiver()),
                             .subtitle = list(waiver()),
                             .tag = list(waiver())) {

    # filt_df <- sim_df |>
    #     filter(outbreaks == "small" & sd_N == 0 & virus_attract == 1 & pseudo_repel == 1 &
    #                ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform") |
    #                     n_pseudo == 0))
    # outcomes = "n_infected"
    # # color_vals = c("black", "dodgerblue", "firebrick2")
    # # col_fct = c("sd_N", "wt_pp")
    # color_vals = NULL
    # col_fct = NULL
    # inters = FALSE
    # fix_y = TRUE; obs_breaks = NULL; incl_vals = TRUE; return_list = FALSE
    # .title = list(waiver()); .subtitle = list(waiver()); .tag = list(waiver())
    # rm(filt_df, outcomes, col_fct, color_vals, inters, fix_y, obs_breaks)
    # rm(incl_vals, return_list, .title, .subtitle, .tag, default_vals)
    # rm(x, unq_x, panel_dims, dsn, dsn_widths, col_scale, line_geom, col_title, n_diff_lvls)

    default_vals <- list(wt_vp = "off *Pseudo.*",
                         wt_pp = "uniform",
                         virus_attract = 1,
                         pseudo_repel = 1,
                         outbreaks = "small",
                         sd_N = 0)

    if (is.null(obs_breaks)) obs_breaks <- scales::breaks_extended(n = 3)

    # Make sure other variables have no extra levels to them:
    for (x in c("outbreaks", "sd_N", "virus_attract", "pseudo_repel",
                "wt_vp", "wt_pp")) {
        if (!is.null(col_fct) && x %in% col_fct) next
        if (x %in% colnames(filt_df)) {
            unq_x <- unique(filt_df[[x]])
            if (length(unq_x[!is.na(unq_x)]) != 1L)
                stop(sprintf("Must have exactly one unique value of `%s`", x))
        }
    }

    stopifnot(all(c("wasp_resp", "n_pseudo", "p_emerge", "outbreak_size", "n_infected") %in%
                      colnames(filt_df)))

    stopifnot(length(outcomes) == 1L && outcomes %in% c("all", "p_emerge", "outbreak_size", "n_infected"))
    if (outcomes == "all") outcomes <- c("p_emerge", "outbreak_size")

    # Plot panel rows, columns, respectively:
    panel_dims <- c(length(outcomes), length(unique(filt_df$wasp_resp)))
    if (!prod(panel_dims) %in% c(1L,2L,4L))
        stop("\n# of unique outcomes-wasp_resp combos must be 1, 2, or 4")
    dsn <- map(1:panel_dims[1], \(i)
               paste(LETTERS[((i-1L)*panel_dims[2]+1L):(i*panel_dims[2])],
                     collapse = "#")) |>
        c(sep = "\n") |> do.call(what = paste)
    dsn_widths <- head(rep(c(1, 0.05), panel_dims[2]), -1)

    # -----------------------------*
    # -----------------------------*
    col_scale <- NULL
    line_geom <- list(geom_line(linewidth = 1),
                      geom_ribbon(aes(ymin = lower, ymax = upper),
                                  fill = "black", alpha = 0.25, color = NA))
    if (!is.null(col_fct)) {

        stopifnot(all(col_fct %in% colnames(filt_df)))

        if (length(col_fct) == 1L) {
            if (!is.factor(filt_df[[col_fct]]))
                filt_df[[col_fct]] <- factor(filt_df[[col_fct]])
            col_title <- pretty_params(col_fct, cap1 = TRUE) |>
                (\(x) {
                    xs <- str_split(x, "\\(")[[1]]
                    xs[1] <- str_replace_all(xs[1], " ", "<br>")
                    str_c(xs, collapse = "(")
                })()
            if (!is.null(color_vals)) {
                stopifnot(length(color_vals) == length(levels(filt_df[[col_fct]])))
                col_scale <- scale_color_manual(col_title, values = color_vals,
                                                aesthetics = c("colour","fill"))
            } else {
                col_scale <- scale_color_scico_d(col_title, end = 0.8,
                                                 palette = "hawaii",
                                                 aesthetics = c("colour","fill"))
            }
            line_geom <- list(geom_line(aes(color = .data[[col_fct]]), linewidth = 1),
                              geom_ribbon(aes(ymin = lower, ymax = upper,
                                              fill = .data[[col_fct]]),
                                          alpha = 0.25, color = NA))

        } else {

            if (!all(col_fct %in% names(default_vals))) {
                stop("col_fct must be one of ", paste(names(default_vals),
                                                      collapse = ", "))
            }


            if (inters) {

                col_title <- "Param(s).<br>that<br>differ:"
                for (x in col_fct) {
                    if (length(unique(filt_df[[x]])) > 2L) stop(x, " has > 2 levels")
                    if (!is.factor(filt_df[[x]])) filt_df[[x]] <- factor(filt_df[[x]])
                    if (! default_vals[[x]] %in% levels(filt_df[[x]])) {
                        stop("levels for ", x, " should contain '",
                             default_vals[[x]], "'")
                    }
                }
                filt_df[["col_fct"]] <- map_chr(1:nrow(filt_df), \(i) {
                    z <- keep(col_fct, \(x) filt_df[[x]][[i]] != default_vals[[x]])
                    if (length(z) == 0) return("none")
                    paste(pretty_params(z, TRUE, serif = TRUE), collapse = "&amp;")
                })

            } else {

                col_title <- "Differences:"

                n_diff_lvls <- rep(0L, nrow(filt_df))
                for (x in col_fct) {
                    if (length(unique(filt_df[[x]])) > 2L) stop(x, " has > 2 levels")
                    if (!is.factor(filt_df[[x]])) filt_df[[x]] <- factor(filt_df[[x]])
                    if (! default_vals[[x]] %in% levels(filt_df[[x]])) {
                        stop("levels for ", x, " should contain '",
                             default_vals[[x]], "'")
                    }
                    n_diff_lvls <- n_diff_lvls + as.integer(filt_df[[x]] != default_vals[[x]])
                }
                filt_df <- filt_df[n_diff_lvls <= 1L,]
                filt_df[["col_fct"]] <- map_chr(1:nrow(filt_df), \(i) {
                    z <- keep(col_fct, \(x) filt_df[[x]][[i]] != default_vals[[x]])
                    if (length(z) == 0) return("none")
                    if (incl_vals) {
                        val <- filt_df[[z]][[i]]
                        if (z == "wt_pp") val <- 3
                        out <- paste(pretty_params(z, TRUE), "=", val) |>
                            serify(prefix = "", suffix = "")
                    } else {
                        out <- pretty_params(z, TRUE, serif = TRUE)
                    }
                    return(out)
                })

            }

            filt_df[["col_fct"]] <- filt_df[["col_fct"]] |>
                (\(x) {
                    x_unq <- unique(x)[unique(x) != "none"]
                    factor(x, levels = c("none", x_unq))
                })()

            if (!is.null(color_vals)) {
                stopifnot(length(color_vals) == length(levels(filt_df[["col_fct"]])))
                col_scale <- scale_color_manual(col_title, values = color_vals,
                                                aesthetics = c("colour","fill"))
            } else {
                col_scale <- scale_color_scico_d(col_title, end = 0.8,
                                                 palette = "hawaii",
                                                 aesthetics = c("colour","fill"))
            }
            line_geom <- list(geom_line(aes(color = col_fct, linewidth = col_fct)),
                              geom_ribbon(aes(ymin = lower, ymax = upper,
                                              fill = col_fct),
                                          alpha = 0.25, color = NA),
                              scale_linewidth_manual(
                                  values = c(1, rep(0.75, length(levels(filt_df[["col_fct"]]))-1)),
                                                     guide = "none"))

            col_fct <- "col_fct"

        }
    }

    if (length(.title) == 1) .title <- rep(.title, prod(panel_dims))
    stopifnot(length(.title) == prod(panel_dims))
    if (length(.subtitle) == 1) .subtitle <- rep(.subtitle, prod(panel_dims))
    stopifnot(length(.subtitle) == prod(panel_dims))
    if (length(.tag) == 1) .tag <- rep(.tag, prod(panel_dims))
    stopifnot(length(.tag) == prod(panel_dims))


    p_list <- filt_df |>
        select(wasp_resp, n_pseudo, all_of(col_fct), p_emerge, outbreak_size,
               n_infected, boots) |>
        pivot_longer(p_emerge:n_infected, names_to = "outcome") |>
        mutate(outcome = factor(outcome, levels = c("p_emerge", "outbreak_size",
                                                    "n_infected"))) |>
        filter(outcome %in% outcomes) |>
        mutate(lower = map2_dbl(boots, outcome, \(b, o) b[[paste(o)]][["Lower"]]),
               upper = map2_dbl(boots, outcome, \(b, o) b[[paste(o)]][["Upper"]])) |>
        (\(dd) {
            max_obs <- NULL
            if (any(dd$outcome != "p_emerge")) max_obs <- max(dd$value, dd$upper)
            dds <- dd |>
                split(~ wasp_resp + outcome, drop = TRUE)
            names(.title) <- names(dds)
            names(.subtitle) <- names(dds)
            names(.tag) <- names(dds)
            dds |>
                imap(\(ddd, n) {
                    yvar <- paste(ddd$outcome[[1]])
                    y_lines <- switch(yvar, p_emerge = c(0, 1),
                                      outbreak_size = 2, n_infected = 1)
                    y_lims <- switch(yvar, p_emerge = c(0, 1.1),
                                     outbreak_size = c(2, max_obs),
                                     n_infected = c(1, max_obs))
                    free_y <- !fix_y && yvar != "p_emerge"
                    if (free_y) y_lims[2] <- NA
                    yb <- switch(yvar, p_emerge = 0:2/2,
                                 outbreak_size = obs_breaks,
                                 n_infected = obs_breaks)
                    y_lab <- yvar_desc[[yvar]] |>
                        (\(x) ifelse(yvar != "p_emerge", paste("mean", x), x))() |>
                        first_cap()
                    p <- ddd |>
                        ggplot(aes(n_pseudo / 10e3 * 100, value)) +
                        geom_hline(yintercept = y_lines, color = "gray70") +
                        line_geom +
                        labs(x = "Percent *Pseudomonas* plants", y = y_lab,
                             title = .title[[n]], subtitle = .subtitle[[n]],
                             tag = .tag[[n]]) +
                        scale_y_continuous(breaks = yb) +
                        scale_x_continuous(breaks = n_pseudo_lvls / 10e3 * 100) +
                        col_scale +
                        coord_cartesian(ylim = y_lims) +
                        theme(plot.margin = margin(0,0,0,0))
                    if (fix_y && panel_dims[2] > 1L &&
                        ddd$wasp_resp[[1]] == levels(wasp_resp_fct)[[2]]) {
                        p <- p + theme(axis.text.y = element_blank())
                    }
                    if (panel_dims[1] > 1L && yvar == "p_emerge") {
                        p <- p + theme(axis.text.x = element_blank())
                    }
                    return(p)
                })
        })()

    if (return_list) return(p_list)

    p_list |>
        wrap_plots() +
        plot_layout(design = dsn, widths = dsn_widths,
                    guides = "collect", axis_titles = "collect")

}



# =============================================================================*
# =============================================================================*
# Basic plot ----
# =============================================================================*
# =============================================================================*

large_outcomes_p <- sim_df |>
    filter(outbreaks == "small" & sd_N == 0 & virus_attract == 1 & pseudo_repel == 1 &
               ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform") |
                    n_pseudo == 0)) |>
    baseline_plotter(outcomes = "n_infected") &
    illustrator_theme



if (.overwrite) {
    save_plot("_plots/large-plantscapes-outcomes.pdf", large_outcomes_p,
              width = 4.88, height = 1.5)
}



# =============================================================================*
# =============================================================================*
# Manips ----
# =============================================================================*
# =============================================================================*

#' - `wt_vp`
#' - `pseudo_repel`
#' - `virus_attract`


# Effects of
# (1) SD in initial aphid abundances and
# (2) clustering vs uniform Pseudomonas:
large_manip_sd_clust_p <- sim_df |>
     filter(outbreaks == "small" & virus_attract == 1 & pseudo_repel == 1 &
                ((wt_vp == "off *Pseudo.*") | n_pseudo == 0)) |>
    add_no_pseudo_points() |>
    baseline_plotter(outcomes = "all", col_fct = c("sd_N", "wt_pp"),
                     color_vals = c("black", "dodgerblue", "firebrick2"),
                     incl_vals = TRUE,
                     .title = list(scenario_title(levels(wasp_resp_fct)[1], TRUE, TRUE),
                                   scenario_title(levels(wasp_resp_fct)[2], TRUE, TRUE),
                                   waiver(), waiver())) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"),
          plot.tag.location = c("panel", "plot", "margin")[1],
          plot.tag.position = c(0.075, 0.92))

# large_manip_sd_clust_p

if (.overwrite) {
    save_plot("_plots/large-plantscapes-baseline-sd_N-wt_pp.pdf",
              large_manip_sd_clust_p, width = 5, height = 5)
}




large_manip_plots <- list(wt_vp = sim_df |>
         filter(outbreaks == "small" & sd_N == 0 & virus_attract == 1 & pseudo_repel == 1 &
                    ((wt_pp == "uniform") | n_pseudo == 0)) |>
         add_no_pseudo_points(),
     pseudo_repel = sim_df |>
         filter(outbreaks == "small" & sd_N == 0 & virus_attract == 1 &
                    ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform") |
                         n_pseudo == 0)),
     virus_attract = sim_df |>
         filter(outbreaks == "small" & sd_N == 0 & pseudo_repel == 1 &
                    ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform") |
                         n_pseudo == 0))) |>
    imap(\(x, n) {
        p <- baseline_plotter(x, outcomes = "all", col_fct = n,
                         color_vals = c("black", par_pal[[n]])) &
            theme(plot.margin = margin(0,0,0,0))
        if (n != "virus_attract") p <- p & theme(axis.text.x = element_markdown(colour = NA))
        return(p)
    })

# wrap_plots(large_manip_plots, ncol = 1, guides = "collect", axes = "collect")



if (.overwrite) {
    for (n in names(large_manip_plots)) {
        p <- large_manip_plots[[n]] & illustrator_theme
        f <- sprintf("_plots/large-plantscapes-baseline-%s.pdf", n)
        save_plot(f, p, width = 4, height = 2.5)
    }; rm(n, p, f)
}


# ----------------------------------*
## Plots of interactions, if interested:
# ----------------------------------*

# sim_df |>
#     filter(outbreaks == "small" & sd_N == 0 & virus_attract == 1 &
#                (wt_pp == "uniform" | n_pseudo == 0)) |>
#     filter(wasp_resp == "weak") |>
#     add_no_pseudo_points() |>
#     select(n_pseudo, wasp_resp, wt_vp, virus_attract, pseudo_repel,
#            p_emerge, outbreak_size) |>
#     baseline_plotter(outcomes = "all",
#                      col_fct = c("wt_vp", "pseudo_repel"),
#                      color_vals = c("black", "gold", "red", "orange"),
#                      inters = TRUE)
# sim_df |>
#     filter(outbreaks == "small" & sd_N == 0 & virus_attract == 1 &
#                (wt_pp == "uniform" | n_pseudo == 0)) |>
#     filter(wasp_resp == "strong") |>
#     add_no_pseudo_points() |>
#     select(n_pseudo, wasp_resp, wt_vp, virus_attract, pseudo_repel,
#            p_emerge, outbreak_size) |>
#     baseline_plotter(outcomes = "all",
#                      col_fct = c("wt_vp", "pseudo_repel"),
#                      color_vals = c("black", "gold", "red", "orange"),
#                      inters = TRUE)





# =============================================================================*
# =============================================================================*
# pseudo_repel when virus starts on *Pseudomonas* ----
# =============================================================================*
# =============================================================================*

pr_onp_p <- sim_df |>
    filter(outbreaks == "small" & sd_N == 0 & virus_attract == 1 &
               ((wt_vp == "on *Pseudo.*" & wt_pp == "uniform") |
                    n_pseudo == 0)) |>
    baseline_plotter(outcomes = "all", col_fct = "pseudo_repel",
                     color_vals = c("black", par_pal[["pseudo_repel"]]),
                     .title = list(scenario_title(levels(wasp_resp_fct)[1], TRUE, TRUE),
                                   scenario_title(levels(wasp_resp_fct)[2], TRUE, TRUE),
                                   waiver(), waiver())) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"),
          plot.tag.location = c("panel", "plot", "margin")[1],
          plot.tag.position = c(0.075, 0.92))
# pr_onp_p

# Note that outbreak sizes above are very noisy because the emergence
# probabilities are very low

if (.overwrite) {
    save_plot("_plots/large-plantscapes-baseline-pseudo_repel-on_pseudo.pdf",
              pr_onp_p, width = 5, height = 5)
}





sim_df |>
    filter(outbreaks == "big_zl" & sd_N == 0 & virus_attract == 1 &
               ((wt_vp == "on *Pseudo.*" & wt_pp == "uniform") |
                    n_pseudo == 0)) |>
    baseline_plotter(outcomes = "all", col_fct = "pseudo_repel",
                     color_vals = c("black", par_pal[["pseudo_repel"]]))
sim_df |>
    filter(outbreaks == "big_zh" & sd_N == 0 & virus_attract == 1 &
               ((wt_vp == "on *Pseudo.*" & wt_pp == "uniform") |
                    n_pseudo == 0)) |>
    baseline_plotter(outcomes = "all", col_fct = "pseudo_repel",
                     color_vals = c("black", par_pal[["pseudo_repel"]]))




va_onp_p <- sim_df |>
    filter(outbreaks == "small" & sd_N == 0 & pseudo_repel == 1 &
               ((wt_vp == "on *Pseudo.*" & wt_pp == "uniform") |
                    n_pseudo == 0)) |>
    baseline_plotter(outcomes = "all", col_fct = "virus_attract",
                     color_vals = c("black", par_pal[["virus_attract"]]),
                     .title = list(scenario_title(levels(wasp_resp_fct)[1], TRUE, TRUE),
                                   scenario_title(levels(wasp_resp_fct)[2], TRUE, TRUE),
                                   waiver(), waiver())) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"),
          plot.tag.location = c("panel", "plot", "margin")[1],
          plot.tag.position = c(0.075, 0.92))

if (.overwrite) {
    save_plot("_plots/large-plantscapes-baseline-virus_attract-on_pseudo.pdf",
              va_onp_p, width = 5, height = 5)
}



# =============================================================================*
# =============================================================================*
# Big outbreaks ----
# =============================================================================*
# =============================================================================*



large_big_manip_plots <- c("big_zl", "big_zh") |>
    set_names() |>
    map(\(ob) {
        y0 <- ifelse(ob == "big_zh", "200", "400")
        .title <- serify("Big outbreaks, ", paste(pretty_params("Y0",TRUE), "=", y0), "")
        list(wt_vp = sim_df |>
                 filter(outbreaks == ob & sd_N == 0 & virus_attract == 1 &
                            pseudo_repel == 1 &
                            ((wt_pp == "uniform") | n_pseudo == 0)) |>
                 add_no_pseudo_points(),
             pseudo_repel = sim_df |>
                 filter(outbreaks == ob & sd_N == 0 & virus_attract == 1 &
                            ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform") |
                                 n_pseudo == 0)),
             virus_attract = sim_df |>
                 filter(outbreaks == ob & sd_N == 0 & pseudo_repel == 1 &
                            ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform") |
                                 n_pseudo == 0))) |>
            imap(\(x, n) {
                ..title <- list(waiver())
                if (n == "wt_vp") ..title <- map(levels(wasp_resp_fct),
                                                 \(wr) scenario_title(wr, TRUE, TRUE))
                pl <- baseline_plotter(x, outcomes = "outbreak_size", col_fct = n,
                                       color_vals = c("black", par_pal[[n]]),
                                       .title = ..title,
                                       fix_y = FALSE,
                                       obs_breaks = scales::breaks_extended(n = 5),
                                       return_list = TRUE)
                if (n != "virus_attract")
                    pl <- c(pl, rep(list(plot_spacer()), length(pl)))
                return(pl)
            }) |>
            do.call(what = c) |>
            wrap_plots(ncol = 2, guides = "collect", axes = "collect",
                       heights = c(1, 0.05, 1, 0.05, 1)) +
            plot_annotation(tag_levels = "A",
                            title = .title,
                            theme = theme(plot.title = element_markdown(
                                face = "bold", size = 16))) &
            theme(legend.position = "right",
                  plot.tag.location = c("panel", "plot", "margin")[1],
                  plot.tag.position = c(0.05, 1.05),
                  axis.title.y = element_markdown(size = 14),
                  axis.title.x = element_markdown(size = 14))
    })


# large_big_manip_plots[["big_zl"]]
# large_big_manip_plots[["big_zh"]]



if (.overwrite) {
    for (n in names(large_big_manip_plots)) {
        f <- sprintf("_plots/large-plantscapes-big-manips-Y0_%s.pdf",
                     ifelse(n == "big_zh", "200", "400"))
        save_plot(f, large_big_manip_plots[[n]], width = 6, height = 7)
    }; rm(n, f)
}


# =============================================================================*
# =============================================================================*
# w=1 ----
# =============================================================================*
# =============================================================================*



# from 07-large-plantscapes-w1.sh:
sim_w1_df <- list.files("_scripts/interm-data", "large-plantscapes-w1-.?.?.rds",
                        full.names = TRUE) |>
    map(\(x) {
        read_rds(x) |>
            # Remove this vector immediately bc it's quite large
            select(-landscape)
    }) |>
    list_rbind() |>
    mutate(outbreak_size = map_dbl(sim, \(x) mean(x$n_infected[x$n_infected > 1])),
           log_outbreak_size = map_dbl(sim, \(x) mean(log10(x$n_infected[x$n_infected > 1]))),
           p_emerge = map_dbl(sim, \(x) mean(x$n_infected > 1))) |>
    mutate(wt_vp = ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") |>
               factor(levels = c("off *Pseudo.*", "on *Pseudo.*")),
           wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
               factor(levels = c("uniform", "clustered")),
           wasp_resp = factor(wasp_resp, levels = levels(wasp_resp_fct)))



large_w1_manip_p <- list(wt_vp = sim_w1_df |>
                                  filter(virus_attract == 1 & pseudo_repel == 1 &
                                             ((wt_pp == "uniform") | n_pseudo == 0)) |>
                                  add_no_pseudo_points(),
                              pseudo_repel = sim_w1_df |>
                                  filter(virus_attract == 1 &
                                             ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform") |
                                                  n_pseudo == 0)),
                              virus_attract = sim_w1_df |>
                                  filter(pseudo_repel == 1 &
                                             ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform") |
                                                  n_pseudo == 0))) |>
    imap(\(x, n) {
        .pal <- c("black", par_pal[[n]])
        if (n == "virus_attract") .pal <- c(.pal, lighten(par_pal[[n]], 0.4))
        .title <- list(waiver())
        if (n == "wt_vp") .title <- map(levels(wasp_resp_fct),
                                        \(wr) scenario_title(wr, TRUE, TRUE))
        pl <- baseline_plotter(x, outcomes = "outbreak_size", col_fct = n,
                               color_vals = .pal, fix_y = FALSE,
                               obs_breaks = scales::breaks_extended(n = 5),
                               .title = .title, return_list = TRUE)
        if (n != "virus_attract")
            pl <- c(pl, rep(list(plot_spacer()), length(pl)))
        return(pl)
    }) |>
    do.call(what = c) |>
    wrap_plots(ncol = 2, guides = "collect", axes = "collect",
               heights = c(1, 0.05, 1, 0.05, 1)) +
    plot_annotation(tag_levels = "A",
                    title = serify("Big outbreaks, ", "*w* = 1", ""),
                    theme = theme(plot.title = element_markdown(
                        face = "bold", size = 16))) &
    theme(legend.position = "right",
          plot.tag.location = c("panel", "plot", "margin")[1],
          plot.tag.position = c(0.05, 1.05),
          axis.title.y = element_markdown(size = 14),
          axis.title.x = element_markdown(size = 14))

# large_w1_manip_p

if (.overwrite) {
    save_plot("_plots/large-plantscapes-big-manips-w1.pdf",
              large_w1_manip_p, width = 6, height = 7)
}





# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# OLD CODE ----
#




big_land_plotter <- function(yvar,
                             wasp_resp,
                             outbreaks = "small",
                             col_fct = "sd_N",
                             x_facet_fct = "land_type",
                             y_facet_fct = c("virus_attract", "pseudo_repel"),
                             shp_lty_fct = NULL,
                             fixed = NULL,
                             color_vals = NULL,
                             facet_scales = "fixed",
                             facet_ncol = NULL,
                             y_breaks = waiver(),
                             y_max = NULL,
                             .title = NULL,
                             .subtitle = waiver(),
                             w1 = FALSE,
                             v100 = FALSE) {

    # yvar = "outbreak_size"; wasp_resp = "strong"; outbreaks = "small"; col_fct = "pseudo_repel"
    # x_facet_fct = "virus_attract"; y_facet_fct = NULL; shp_lty_fct = "wt_vp"
    # fixed = list(wt_pp = "uniform", sd_N = 0)
    # color_vals = NULL; facet_scales = "fixed"; facet_ncol = NULL
    # y_breaks = waiver(); y_max = NULL; .title = NULL; .subtitle = waiver()
    # w1 = FALSE; v100 = FALSE
    # rm(yvar, wasp_resp, outbreaks, col_fct, x_facet_fct, y_facet_fct, shp_lty_fct, fixed)
    # rm(color_vals, facet_scales, facet_ncol, y_breaks, y_max, .title, .subtitle, w1, v100)
    # rm(pretty_greek_fct, make_land_type, facet_scales, dd, facet_fct)
    # rm(points_lines, facet_form, facets_coord, col_scale, y_ints, y_lab)

    pretty_greek_fct <- function(x, letter, suffix = "") {
        lvls <- sort(unique(x))
        factor(x, levels = lvls,
               labels = serify("", paste0("&", letter, ";", suffix),
                               sprintf(" = %.1f", lvls)))
    }
    make_land_type <- function(wt_vp, wt_pp, sep = "<br>") {
        lvls <- paste(rep(c("off *Pseudo.*", "on *Pseudo.*"), each = 2),
                      rep(c("uniform", "clustered"), 2), sep = "<br>")
        factor(paste0(wt_vp, sep, wt_pp), levels = lvls)
    }

    if (v100 && !w1) stop("v100 implies w1")
    if ((v100 || w1) && outbreaks == "small")
        stop("outbreaks are large when v100 or w1")

    if (is.null(.title)) {
        .title <- scenario_title(wasp_resp, TRUE, TRUE)
    }

    facet_scales <- match.arg(facet_scales, c("fixed", "free", "free_x", "free_y"))

    if (w1) {
        dd <- sim_w1_df |>
            filter(wasp_resp == .env$wasp_resp, virus_attract < 100) |>
            mutate(sd_N = 0)
    } else if (v100) {
        dd <- sim_w1_df |>
            filter(wasp_resp == .env$wasp_resp, virus_attract %in% c(1, 100)) |>
            mutate(sd_N = 0)
    } else dd <- sim_df |>
        filter(wasp_resp == .env$wasp_resp, outbreaks == .env$outbreaks) |>
        select(-outbreaks)

    if (nrow(dd) == 0) stop("Some combination of args resulted in no output from sim_df")

    dd <- dd |>
        mutate(outbreak_size = ifelse(is.na(outbreak_size), -Inf, outbreak_size)) |>
        select(-wasp_resp, -sim) |>
        add_no_pseudo_points() |>
        mutate(land_type = make_land_type(wt_vp, wt_pp),
               grp = interaction(sd_N, land_type, virus_attract,
                                 pseudo_repel, drop = TRUE))
    if (!is.null(fixed)) {
        stopifnot(is.list(fixed) && !is.null(names(fixed)) && all(names(fixed) != ""))
        stopifnot(all(names(fixed) %in% colnames(dd)))
        stopifnot(!any(names(fixed) %in% c(col_fct, x_facet_fct, y_facet_fct, shp_lty_fct)))
        for (i in 1:length(fixed)) {
            .col <- names(fixed)[[i]]
            .val <- fixed[[i]]
            dd <- dd[dd[[.col]] == .val, ]
        }; rm(i, .col, .val)
        dd <- dd |> select(-all_of(names(fixed)))
    }
    facet_fct <- c(x_facet_fct, y_facet_fct)

    if ("sd_N" %in% facet_fct) dd$sd_N <- pretty_greek_fct(dd$sd_N, "sigma",
                                                           "<sub>N</sub>")
    if ("pseudo_repel" %in% facet_fct) dd$pseudo_repel <- pretty_greek_fct(
        dd$pseudo_repel, "rho")
    if ("virus_attract" %in% facet_fct) dd$virus_attract <- pretty_greek_fct(
        dd$virus_attract, "nu")
    dd[[col_fct]] <- factor(dd[[col_fct]])
    if (!is.null(shp_lty_fct)) {
        dd[[shp_lty_fct]] <- factor(dd[[shp_lty_fct]])
        points_lines <- list(geom_point(aes(shape = .data[[shp_lty_fct]])),
                             geom_line(aes(group = grp, linetype = .data[[shp_lty_fct]])))
                             # geom_line(aes(group = grp, linewidth = .data[[shp_lty_fct]])))
        shp_lty_ttl <- pretty_params(shp_lty_fct, cap1 = TRUE) |>
            str_replace_all("\\s+", "<br>")
    } else {
        points_lines <- list(geom_point(), geom_line(aes(group = grp)))
        shp_lty_ttl <- NULL
    }

    if (length(facet_fct) > 0) {
        facet_form <- paste(paste(y_facet_fct, collapse = "+"), "~",
                            paste(x_facet_fct, collapse = "+")) |>
            (\(x) {
                if (grepl("~ $", x)) x <- paste0(x, ".")
                return(x)
            })() |>
            as.formula()
        if (facet_scales == "free" || facet_scales == "free_y") {
            facets_coord <- facet_wrap(facet_form, scales = facet_scales,
                                       ncol = facet_ncol)
        } else {
            if (is.null(y_max)) y_max <- max(sim_df[[yvar]])
            facets_coord <- list(facet_grid(facet_form),
                                 coord_cartesian(ylim = c(0, y_max)))
        }

    } else {
        facets_coord <- NULL
        if (facet_scales != "free" && facet_scales != "free_y") {
            if (is.null(y_max)) y_max <- max(sim_df[[yvar]])
            facets_coord <- coord_cartesian(ylim = c(0, y_max))
        }
    }

    if (!is.null(color_vals)) {
        col_scale <- scale_color_manual(pretty_params(col_fct, cap1 = TRUE,
                                                      serif = TRUE) |>
                                            str_replace(" ", "<br>"),
                                        values = color_vals)
    } else {
        col_scale <- scale_color_viridis_d(pretty_params(col_fct, cap1 = TRUE,
                                                         serif = TRUE) |>
                                               str_replace(" ", "<br>"),
                                           begin = 0.2, end = 0.8, option = "plasma")
    }

    y_ints <- 0
    if (yvar == "p_emerge") y_ints <- c(0, 1)
    if (yvar == "outbreak_size" && (outbreaks == "big" || w1)) y_ints <- c(0, 10e3)

    y_lab <- yvar_desc[[yvar]] |>
        (\(x) ifelse(yvar == "outbreak_size", paste("mean", x), x))() |>
        first_cap()

    dd |>
        ggplot(aes(n_pseudo / 10e3 * 100, .data[[yvar]], color = .data[[col_fct]])) +
        geom_hline(yintercept = y_ints, color = "gray70") +
        points_lines +
        col_scale +
        scale_linetype_manual(shp_lty_ttl, values = c("solid", "22")) +
        # scale_linewidth_manual(shp_lty_ttl, values = c(1, 0.5)) +
        scale_shape_manual(shp_lty_ttl, values = c(19, 17)) +
        scale_x_continuous(breaks = c(0, 10, 30, 50, 70, 90)) +
        scale_y_continuous(breaks = y_breaks) +
        facets_coord +
        labs(x = "Percent *Pseudomonas* plants",
             y = y_lab, title = .title, subtitle = .subtitle) +
        theme(legend.text = element_markdown(),
              plot.subtitle = element_markdown())
}



#'
#' Possible values for *_fct arguments:
#'
#' - `wt_vp`
#' - `wt_pp`
#' - `sd_N`
#' - `virus_attract`
#' - `pseudo_repel`
#'


#'
#' Effects of `sd_N` and `wt_pp`:
#'
# .y <- "outbreak_size"
# .t <- "strong"
# .o <- "small"
# big_land_plotter(yvar = .y, wasp_resp = .t, outbreaks = .o,
#                  col_fct = "wt_pp",
#                  shp_lty_fct = "sd_N",
#                  x_facet_fct = c("pseudo_repel", "wt_vp"),
#                  y_facet_fct = "virus_attract",
#                  facet_scales = "free", facet_ncol = 4L)


#'
#' Effects of others:
#'
# .y <- "outbreak_size"
# .t <- "strong"
# .o <- "big"
# big_land_plotter(yvar = .y, wasp_resp = .t, outbreaks = .o,
#                  col_fct = "pseudo_repel",
#                  shp_lty_fct = "wt_vp",
#                  x_facet_fct = NULL,
#                  y_facet_fct = "virus_attract",
#                  fixed = list(wt_pp = "uniform", sd_N = 0),
#                  color_vals = scico(2, end = 0.8, palette = "hawaii"),
#                  facet_scales = "free")







# =============================================================================*
# =============================================================================*
# Results ----
# =============================================================================*
# =============================================================================*
#'
#' # =============================================
# _ Outbreak size ----
#' # =============================================
#'
#' # ------------
#' When `wasp_resp = "weak"`
#' # ------------
#'
#' * `pseudo_repel` (rho):
#'     - `pseudo_repel` = 5 increases outbreak sizes when virus starts off
#'       *Pseudomonas*
#'     - `pseudo_repel` = 5 decreases outbreak sizes or has little effect
#'       when virus starts on *Pseudomonas*
#' * `virus_attract` (nu):
#'     - `virus_attract` > 1 always increases outbreak sizes, quite strongly
#'     - when `virus_attract` = 5, there is an especially large drop in the
#'       effect of *Pseudomonas* on outbreak sizes at 90% *Pseudomonas* plants
#'       compared to 70%
#'         - this is especially pronounced virus starts on *Pseudomonas*  and
#'           when landscape is uniform
#' * `sd_N`:
#'     - `sd_N` > 0 can cause increase in outbreak sizes at higher
#'       *Pseudomonas* densities; this effect is slightly more consistent
#'       in clustered landscapes
#' * `wt_vp` (Pseudomonas location in relation to virus):
#'     - outbreak sizes always lower when virus starts on *Pseudomonas*
#'     - this effect is strongest when `pseudo_repel` = 5.0
#'     - this effect is stronger when in a clustered landscape only
#'       when `virus_attract` = 1
#'     - can cause *Pseudomonas* to have little effect on outbreak size,
#'       or to even slightly decrease them
#' * `wt_pp` (Pseudomonas clustering):
#'     - clustered causes larger outbreaks at high *Pseudomonas* densities
#'     - usually (but not always) strongest when `sd_N` = 50
#'     - slightly stronger when `virus_attract` or `pseudo_repel` equals 1,
#'       but not when both equal 1
#'
#'
#'
#' # ------------
#' When `wasp_resp = "strong"`
#' # ------------
#'
#' * `pseudo_repel` (rho):
#'     - `pseudo_repel` = 5 increases outbreak sizes when virus starts off
#'       *Pseudomonas*
#'     - when virus starts on *Pseudomonas*, `pseudo_repel` = 5 decreases
#'       outbreak sizes at strong *Pseudomonas*  densities and has little effect
#'       at higher densities (>=50%)
#'     - both effects are more consistent when `pseudo_repel` = 5
#'     - neither effect alters much regarding the effect of *Pseudomonas* on
#'       outbreak size when measuring between `n_pseudo` = 0 and 9000
#' * `virus_attract` (nu):
#'     - greater value increases outbreak size significantly
#'     - this effect dissipates with increasing *Pseudomonas* densities because
#'       *Pseudomonas* reduces outbreak sizes overall, such that at very high
#'       *Pseudomonas* densities (90%), the effect of `virus_attract`
#'       is negligible
#' * `sd_N`:
#'     - `sd_N` > 1 can slightly reduce outbreak size, especially at lower
#'       *Pseudomonas* densities, so it slightly reduces the outbreak-inhibiting
#'       effect of *Pseudomonas*.
#'     - effect is weak
#' * `wt_vp` (Pseudomonas location in relation to virus):
#'     - outbreak size reduced when virus starts on *Pseudomonas*
#'     - this effect is most pronounced when `pseudo_repel` = 5
#'     - this effect obviously only happens when there is *Pseudomonas* on
#'       landscape, plus its effect levels off at very high *Pseudomonas*
#'       densities (90%) since at these densities, outbreak sizes are so small
#'       across the board
#'     - this results in a more convex outbreak ~ *Pseudomonas* curve instead
#'       of linear
#' * `wt_pp` (Pseudomonas clustering):
#'     - little to no effect
#'
#'
#'
#'
#' # =============================================
# _ Prob. emergence ----
#' # =============================================
#'
#' # ------------
#' When `wasp_resp = "weak"`
#' # ------------
#' * `pseudo_repel` (rho):
#'     - `pseudo_repel` = 5 decreases chance of outbreak when virus starts on
#'       *Pseudomonas*
#'     - `pseudo_repel` = 5 has little effect when virus starts off *Pseudomonas*
#' * `virus_attract` (nu):
#'     - `virus_attract` > 1 always increases chance of outbreak
#'     - this effect is strongest when virus starts on *Pseudomonas*
#' * `sd_N`:
#'     - effect is negligible
#' * `wt_vp` (Pseudomonas location in relation to virus):
#'     - outbreaks far less likely when virus starts on *Pseudomonas*
#'     - effect is strongest when `virus_attract` = 1 and `pseudo_repel` = 5
#' * `wt_pp` (Pseudomonas clustering):
#'     - effect is negligible
#'
#'
#'
#' # ------------
#' When `wasp_resp = "strong"`
#' # ------------
#'
#' * `pseudo_repel` (rho):
#'     - `pseudo_repel` = 5 increases chance of outbreak when virus starts off
#'       *Pseudomonas*, especially at higher *Pseudomonas* densities (>= 50%)
#'     - `pseudo_repel` = 5 decreases chance of outbreak when virus starts
#'       on *Pseudomonas*, but this effect levels off with greater
#'       *Pseudomonas* such that at 90% *Pseudomonas* plants, it is gone
#'       entirely
#' * `virus_attract` (nu):
#'     - greater value always increases chance of outbreak
#'     - effect is strongest when virus starts on *Pseudomonas*
#' * `sd_N`:
#'     - `sd_N` > 1 can slightly reduce chances of outbreaks, especially at
#'       lower *Pseudomonas* densities, so it slightly reduces the
#'       outbreak-inhibiting effect of *Pseudomonas*.
#'     - effect is weak
#' * `wt_vp` (Pseudomonas location in relation to virus):
#'     - chance of outbreak reduced when *Pseudomonas* starts on virus
#'     - this effect is most pronounced when `virus_attract` = 1 and
#'       when `pseudo_repel` = 5; when both are true, it's strongest
#' * `wt_pp` (Pseudomonas clustering):
#'     - effect is negligible
#'
#'
#'
#' # =============================================
#' From this, I think it makes sense to...
#' # =============================================
#'
#'
#' use these in main text:
#' - `pseudo_repel`
#' - `virus_attract`
#' - `wt_vp`
#'
#' and relegate `sd_N` and `wt_pp` to the supplement
#'
#'



















# =============================================================================*
# =============================================================================*
# Save plots ----
# =============================================================================*
# =============================================================================*

# crossing(.t = levels(wasp_resp_fct),
#          .v = c(1, 5)) |>
#     pmap(\(.t, .v) {
#         .vt <- serify("", sprintf("%s = %.0f", pretty_params("virus_attract", TRUE), .v), "")
#         .tt <- paste0(scenario_title(.t, TRUE, TRUE), "<br>", .vt)
#         if (.v == 5) .tt <- .vt
#     big_land_plotter(yvar = "p_emerge", wasp_resp = .t,
#                      w1 = FALSE,
#                      col_fct = "pseudo_repel",
#                      shp_lty_fct = "wt_vp",
#                      x_facet_fct = NULL,
#                      y_facet_fct = NULL,
#                      .title = .tt,
#                      fixed = list(wt_pp = "uniform", sd_N = 0,
#                                   virus_attract = .v),
#                      color_vals = scico(2, end = 0.8, palette = "hawaii"))
# }) |>
#     do.call(what = wrap_plots) +
#     plot_layout(nrow = 2, axis_titles = "collect", guides = "collect", byrow = FALSE)





if (.overwrite) {

    for (.t in levels(wasp_resp_fct)) {

        .p <- tibble(yvar = c("p_emerge", "outbreak_size")) |>
            mutate(virus_attract = map(1:n(), \(i) sort(unique(sim_df$virus_attract)))) |>
            unnest(virus_attract) |>
            pmap(\(yvar, virus_attract) {
                .ym <- if (yvar == "outbreak_size") NA else NULL
                .yb <- if (yvar == "outbreak_size") waiver() else c(0, 0.5, 1)
                big_land_plotter(yvar = yvar, wasp_resp = .t, outbreaks = "small",
                                 col_fct = "pseudo_repel",
                                 shp_lty_fct = "wt_vp",
                                 x_facet_fct = NULL,
                                 y_facet_fct = NULL,
                                 fixed = list(wt_pp = "uniform", sd_N = 0,
                                              virus_attract = virus_attract),
                                 color_vals = scico(2, end = 0.8, palette = "hawaii"),
                                 y_max = .ym, y_breaks = .yb) +
                    theme(plot.margin = margin(0, 0, 0, 0))
            }) |>
            do.call(what = wrap_plots) +
            plot_layout(design = "A#B\n###\nC#D", heights = c(1, 0.15, 1),
                        widths = c(1, 0.1, 1),
                        axis_titles = "collect", guides = "collect", byrow = TRUE) &
            illustrator_theme
        .f <- sprintf("_plots/large-plantscapes-%s.pdf", .t)
        save_plot(.f, .p, width = 3.25, height = 3)
    }; rm(.t, .p, .f)
}






# =============================================================================*
# =============================================================================*
# Big outbreaks  ----
# =============================================================================*
# =============================================================================*


big_outbreaks_p <- crossing(.wasp_resp = levels(wasp_resp_fct),
                            .outbreaks = c("big_zh", "big_zl")) |>
    arrange(.outbreaks, .wasp_resp) |>
    pmap(\(.wasp_resp, .outbreaks) {

        .Y0 <- ifelse(.outbreaks == "big_zh", 200, 400)
        .yt <- serify("",sprintf("%s = %.0f", pretty_params("Y0",TRUE), .Y0),"")

        .tt <- scenario_title(.wasp_resp, TRUE, TRUE)
        if (.outbreaks == "big_zl") .tt <- waiver()

        big_land_plotter(yvar = "outbreak_size", wasp_resp = .wasp_resp,
                         outbreaks = .outbreaks,
                         col_fct = "pseudo_repel",
                         shp_lty_fct = "wt_vp",
                         x_facet_fct = "virus_attract",
                         y_facet_fct = NULL,
                         fixed = list(wt_pp = "uniform", sd_N = 0),
                         color_vals = scico(2, end = 0.8, palette = "hawaii"),
                         .title = .tt,
                         .subtitle = .yt,
                         y_max = NA) +
            theme(plot.margin = margin(0, 0, 0, 0))
    }) |>
    do.call(what = wrap_plots) +
    plot_layout(design = "A#B\nC#D", widths = c(1, 0.1, 1),
                axis_titles = "collect", guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"),
          plot.subtitle = element_markdown(size = 14))

# big_outbreaks_p

if (.overwrite) {
    save_plot("_plots/large-plantscapes-big.pdf", big_outbreaks_p,
              width = 7.5, height = 7)
}




# =============================================================================*
# =============================================================================*
# Effect of virus_attract  ----
# =============================================================================*
# =============================================================================*

virus_attract_plotter <- function(wasp_resp, yvar, outbreaks = "small", w1 = FALSE,
                                  v100 = FALSE,
                                  .title = NULL, .subtitle = waiver(),
                                  xvar = "pseudo_percent",
                                  y_lims = NULL, x_lims = NULL,
                                  y_lines = NULL, log_rel = TRUE) {
    # wasp_resp = "weak"; yvar = "outbreak_size"; outbreaks = "big"; w1 = TRUE; v100 = FALSE
    # .title = waiver(); xvar = "log_alates"
    # rm(wasp_resp, yvar, outbreaks, w1, v100, .title, .subtitle)
    # rm(xvar, y_lims, x_lims, y_lines, log_rel)
    # rm(y_line, y, y_lab, shp_lty_ttl, dd)

    stopifnot(xvar %in% c("pseudo_percent", "mean", colnames(sim_df$sim[[1]])[-1]))

    if (v100 && !w1) stop("v100 implies w1")
    if ((v100 || w1) && outbreaks == "small")
        stop("outbreaks are large when v100 or w1")

    wasp_resp <- paste(wasp_resp)
    outbreaks <- paste(outbreaks)
    yvar <- paste(yvar)
    if (is.null(.title)) {
        .title <- scenario_title(wasp_resp, TRUE, TRUE)
    }

    y <- "rel"
    y_lab <- paste(first_cap(yvar_desc[[yvar]]),
                   "with<br>virus attract / without")
    shp_lty_ttl <- pretty_params("wt_vp", cap1 = TRUE) |>
        str_replace(" ", "<br>")
    if (log_rel) {
        y <- "log_rel"
        y_lab <- paste0("log<sub>2</sub>(", yvar_desc[[yvar]],
                        " with<br>virus attract / without)")
    }
    if (w1) {
        dd <- sim_w1_df |>
            filter(wasp_resp == .env$wasp_resp, virus_attract < 100,
                   wt_pp %in% c(NA, "uniform"))
    } else if (v100) {
        dd <- sim_w1_df |>
            filter(wasp_resp == .env$wasp_resp, virus_attract %in% c(1, 100),
                   wt_pp %in% c(NA, "uniform"))
    } else {
        dd <- sim_df |>
            filter(wasp_resp == .env$wasp_resp,
                   outbreaks == .env$outbreaks,
                   wt_pp %in% c(NA, "uniform"),
                   sd_N == 0)
    }

    if (xvar %in% colnames(sim_df$sim[[1]])[-1]) {
        dd[[xvar]] <- map_dbl(dd$sim, \(x) mean(x[[xvar]]))
    }
    x_breaks <- if (xvar == "pseudo_percent") { c(0, 10, 30, 50, 70, 90)
        } else {
        waiver()
    }
    x_lab <- if (xvar == "pseudo_percent") {"Percent *Pseudomonas* plants"
    } else if (xvar == "mean") paste("Mean", yvar_desc[[yvar]]) else {
        first_cap(yvar_desc[[xvar]])
    }

    dd |>
        add_no_pseudo_points() |>
        mutate(virus_attract = case_when(virus_attract == max(virus_attract) ~ "big",
                                         virus_attract == min(virus_attract) ~ "small",
                                         .default = NA)) |>
        select(wt_vp, n_pseudo, pseudo_repel, virus_attract, any_of(c(yvar, xvar))) |>
        pivot_wider(names_from = "virus_attract", values_from = yvar,
                    id_cols = wt_vp:pseudo_repel, unused_fn = mean) |>
        mutate(rel = big / small,
               log_rel = log2(big / small),
               mean = (big + small) / 2,
               pseudo_percent = n_pseudo / 10e3 * 100) |>
        mutate(pseudo_repel = factor(pseudo_repel)) |>
        ggplot(aes(.data[[xvar]], .data[[y]], color = pseudo_repel)) +
        geom_hline(yintercept = y_lines, color = "gray70") +
        geom_point(aes(shape = wt_vp)) +
        geom_line(aes(linetype = wt_vp)) +
        scale_x_continuous(breaks = x_breaks) +
        labs(x = x_lab, y = y_lab, title = .title, subtitle = .subtitle) +
        coord_cartesian(ylim = y_lims, xlim = x_lims) +
        scale_linetype_manual(shp_lty_ttl, values = c("solid", "22")) +
        scale_shape_manual(shp_lty_ttl, values = c(19, 17)) +
        scale_color_manual(pretty_params("pseudo_repel") |>
                               str_replace(" ", "<br>"),
                           values = scico(2, end = 0.8, palette = "hawaii")) +
        theme(plot.subtitle = element_markdown())
}

virus_attract_p <- crossing(yvar = c("p_emerge", "outbreak_size") |>
                                (\(x) factor(x, levels = x))(),
                            wasp_resp = wasp_resp_fct) |>
    mutate(.title = map(yvar, \(y) {
        if (y == "p_emerge") return(NULL) else return(waiver())})) |>
    mutate(y_lines = 0) |>
    pmap(virus_attract_plotter) |>
    wrap_plots(nrow = 2, guides = "collect", axis_titles = "collect", byrow = TRUE) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))

big_virus_attract_p <- crossing(wasp_resp = wasp_resp_fct,
                                outbreaks = sim_df$outbreaks |> unique() |>
                                    sort() |> tail(-1) |> fct_drop()) |>
    arrange(outbreaks, wasp_resp) |>
    mutate(.title = map(outbreaks, \(ob) {
        if (ob == levels(sim_df$outbreaks)[2]) return(NULL) else return(waiver())}),
        .subtitle = map_chr(outbreaks, \(ob) {
            .Y0 <- ifelse(ob == "big_zh", 200, 400)
            serify("",sprintf("%s = %.0f", pretty_params("Y0",TRUE), .Y0),"")})) |>
    mutate(yvar = "outbreak_size", y_lines = 0) |>
    pmap(virus_attract_plotter) |>
    wrap_plots(nrow = 2, guides = "collect", axis_titles = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))

# virus_attract_p
# big_virus_attract_p

if (.overwrite) {
    save_plot("_plots/virus_attract-small.pdf", virus_attract_p,
              width = 7.5, height = 6.5)
    save_plot("_plots/virus_attract-big.pdf", big_virus_attract_p,
              width = 7.5, height = 6.5)
}




#
# Effect of virus_attract = 5 when w = 1? ----
#

va_w1_obs_plots <- map(levels(wasp_resp_fct), \(.t) {
    big_land_plotter(yvar = "outbreak_size", wasp_resp = .t, outbreaks = "big",
                     col_fct = "pseudo_repel",
                     shp_lty_fct = "wt_vp",
                     x_facet_fct = "virus_attract",
                     y_facet_fct = NULL,
                     fixed = list(wt_pp = "uniform", sd_N = 0),
                     color_vals = scico(2, end = 0.8, palette = "hawaii"),
                     w1 = TRUE,
                     y_max = NA)
})

va_w1_vaf_plots <- crossing(x = c("pseudo_percent", "mean") |>
                                (\(x) factor(x, levels = x))(),
                            ty = wasp_resp_fct) |>
    pmap(\(x, ty) {
        virus_attract_plotter(yvar = "outbreak_size", wasp_resp = paste(ty),
                              w1 = TRUE, outbreaks = "big", .title = waiver(),
                              xvar = paste(x), x_lims = NULL, y_lines = 0)
    })


va_w1_p <- c(va_w1_obs_plots, va_w1_vaf_plots) |>
    do.call(what = wrap_plots) +
    plot_layout(ncol = 2, axis_titles = "collect", guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))


if (.overwrite) {
    save_plot("_plots/virus_attract-w1.pdf", va_w1_p, width = 7.5, height = 6.5)
}

#
# Effect of virus_attract = 100 when w = 1? ----
#

va_v100_obs_plots <- map(levels(wasp_resp_fct), \(.t) {
    big_land_plotter(yvar = "outbreak_size", wasp_resp = .t, outbreaks = "big",
                     col_fct = "pseudo_repel",
                     shp_lty_fct = "wt_vp",
                     x_facet_fct = "virus_attract",
                     y_facet_fct = NULL,
                     fixed = list(wt_pp = "uniform", sd_N = 0),
                     color_vals = scico(2, end = 0.8, palette = "hawaii"),
                     v100 = TRUE, w1 = TRUE,
                     y_max = NA)
})

va_v100_vaf_plots <- crossing(x = c("pseudo_percent", "mean") |>
                                  (\(x) factor(x, levels = x))(),
                            ty = wasp_resp_fct) |>
    pmap(\(x, ty) {
        virus_attract_plotter(yvar = "outbreak_size", wasp_resp = paste(ty),
                              outbreaks = "big",
                              v100 = TRUE, w1 = TRUE, .title = waiver(),
                              xvar = paste(x), x_lims = NULL, y_lines = 0)
    })


va_v100_p <- c(va_v100_obs_plots, va_v100_vaf_plots) |>
    do.call(what = wrap_plots) +
    plot_layout(ncol = 2, axis_titles = "collect", guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))

# This figure doesn't say anything that's already said in the w=1 plot
# if (.overwrite) {
#     save_plot("_plots/virus_attract-v100.pdf", va_v100_p, width = 6.5, height = 6.5)
# }



# =============================================================================*
# =============================================================================*
# Why pseudo_repel = 5 increase outbreaks when virus starts off Pseudomonas? ----
# =============================================================================*
# =============================================================================*

.wasp_resp <- "weak"
.sd_N <- 0
.virus_attract <- 1

.landscape <- read_rds("_scripts/interm-data/large-landscapes.rds") |>
    filter(n_pseudo == 7000,
           wt_pp == 1,
           ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") == "off *Pseudo.*") |>
    getElement("landscape") |> getElement(1) |>
    base::`[`(,,1:10)

# Takes ~30 sec for both
set.seed(1083022697)
sims_p1 <- large_simmer(.landscape, .wasp_resp, .sd_N, .virus_attract,
                        pseudo_repel = 1, n_sims = dim(.landscape)[3], summ = "all",
                        max_t = 50,
                        out_dispersals = TRUE)
sims_p5 <- large_simmer(.landscape, .wasp_resp, .sd_N, .virus_attract,
                        pseudo_repel = 5, n_sims = dim(.landscape)[3], summ = "all",
                        max_t = 50,
                        out_dispersals = TRUE)


# In the `disps` column the column indicates the plant the alate came from,
# and the row indicates the plant the alate dispersed to.
#
# Below, this is the total number of incoming alates to the initially
# infected plant:

tibble(p1 = sims_p1$disps |> map_int(\(x) sum(x[1,])),
       p5 = sims_p5$disps |> map_int(\(x) sum(x[1,])))
# # A tibble: 10 × 2
#       p1    p5
#    <int> <int>
#  1    46   138
#  2    47   133
#  3    74   153
#  4    55   140
#  5    53   131
#  6    52   123
#  7    46   130
#  8    61   108
#  9    62   135
# 10    57   149




# =============================================================================*
# =============================================================================*
# Why virus_attract has greatest effect when Pseudomonas at mid densities (at wasp_resp = weak)? ----
# =============================================================================*
# =============================================================================*

.wasp_resp <- "weak"
.sd_N <- 0
.pseudo_repel <- 1

.landscapes <- read_rds("_scripts/interm-data/large-landscapes.rds") |>
    filter(n_pseudo %in% c(1000, 7000, 9000),
           wt_pp == 1,
           ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") == "on *Pseudo.*") |>
    arrange(n_pseudo) |>
    getElement("landscape") |>
    set_names(c("low", "med", "high")) |>
    map(\(x) x[,,1:10,drop=FALSE])

# Takes ~2 min
set.seed(500860318)
sims_va <- imap(.landscapes, \(l, n) {

    sims_va1 <- large_simmer(l, .wasp_resp, .sd_N, virus_attract = 1,
                             pseudo_repel = .pseudo_repel, n_sims = dim(l)[3],
                             summ = "all", # max_t = 50,
                             out_dispersals = TRUE) |>
        #
        # In the `disps` column the column indicates the plant the alate came
        # from, and the row indicates the plant the alate dispersed to.
        #
        # The first line below is the total number of incoming alates to the
        # initially infected plant:
        #
        mutate(disps_in = map_int(disps, \(x) sum(x[1,])),
               disps_out = map_int(disps, \(x) sum(x[,1]))) |>
        select(-disps) |>
        mutate(virus_attract = 1)
    sims_va5 <- large_simmer(l, .wasp_resp, .sd_N, virus_attract = 5,
                             pseudo_repel = .pseudo_repel, n_sims = dim(l)[3],
                             summ = "all", # max_t = 50,
                             out_dispersals = TRUE) |>
        mutate(disps_in = map_int(disps, \(x) sum(x[1,])),
               disps_out = map_int(disps, \(x) sum(x[,1]))) |>
        select(-disps) |>
        mutate(virus_attract = 5)
    bind_rows(sims_va1, sims_va5) |>
        mutate(pseudo = factor(n, levels = names(.landscapes))) |>
        select(pseudo, virus_attract, rep, everything())
}, .progress = .prog_args) |>
    list_rbind()


sims_va |>
    group_by(pseudo, virus_attract) |>
    summarize(across(starts_with("disps"), mean), .groups = "drop")
# # A tibble: 6 × 4
#   pseudo virus_attract disps_in disps_out
#   <fct>          <dbl>    <dbl>     <dbl>
# 1 low                1     50.3      41.6
# 2 low                5    232.      190.
# 3 med                1     52.8      41.6
# 4 med                5    250       202.
# 5 high               1     20.7      17.8
# 6 high               5     99.8      81.4


sims_va |>
    group_by(pseudo, virus_attract) |>
    summarize(across(starts_with("disps"), mean), .groups = "drop") |>
    group_by(pseudo) |>
    summarize(across(starts_with("disps"), \(x) x[virus_attract == 5] / x[virus_attract == 1]))
# # A tibble: 3 × 3
#   pseudo disps_in disps_out
#   <fct>     <dbl>     <dbl>
# 1 low        4.61      4.56
# 2 med        4.73      4.85
# 3 high       4.82      4.57








# =============================================================================*
# =============================================================================*
# Why does clustering increase outbreak sizes? ----
# =============================================================================*
# =============================================================================*


land_plotter <- function(lands, .title = waiver(), .expand_axes = TRUE) {

    stopifnot(is.array(lands) || is.list(lands))

    if (is.array(lands)) {
        stopifnot(length(dim(lands)) == 3L)
        n_lands <- dim(lands)[3]
        lands <- map(1:n_lands, \(i) lands[,,i])
    } else n_lands <- length(lands)

    # combining pseudo and virus colors for both:
    both_col <- col2rgb(color_pal[c("pseudo", "virus")]) |>
        apply(1, mean) |>
        (\(x) x / 255)() |>
        as.list() |>
        do.call(what = rgb)
    type_pal <- c("gray80", color_pal[c("virus", "pseudo")], both_col) |>
        set_names(c("none", "virus", "pseudo", "both"))

    if (.expand_axes) {
        .plot_margin <- margin(0,0,0,r=6)
        .axis_expand <- waiver()
    } else {
        .plot_margin <- margin(0,0,0,0)
        .axis_expand <- c(0,0)
    }

    p <- lands |>
        imap(\(x, i) {
            colnames(x) <- paste(1:ncol(x))
            x |>
                as_tibble() |>
                mutate(x = 1:n()) |>
                select(x, everything()) |>
                pivot_longer(-x, names_to = "y", values_to = "type") |>
                mutate(y = as.integer(y),
                       type = factor(type, levels = 0:3,
                                     labels = c("none", "virus", "pseudo", "both")),
                       rep = i)
        }) |>
        list_rbind() |>
        mutate(rep = factor(rep)) |>
        ggplot(aes(x, y, fill = type)) +
        geom_raster(show.legend = TRUE) +
        scale_fill_manual(values = type_pal, drop = FALSE) +
        scale_x_continuous(expand = .axis_expand) +
        scale_y_reverse(expand = .axis_expand) +
        theme_void() +
        theme(plot.margin = .plot_margin,
              plot.title = element_markdown(size = 16),
              strip.text = element_markdown(size = 10, lineheight = 0.7)) +
        labs(title = .title) +
        coord_equal()
    if (n_lands > 1) p <- p + facet_wrap(~ rep, nrow = round(sqrt(n_lands)))

    return(p)

}


circle_maker <- function(center = c(0,0), radius = 0.5, npoints = 100,
                      rad_min = 0, rad_max = 2 * pi){
    tt <- seq(rad_min, rad_max, length.out = npoints)
    xx <- center[1] + radius * cos(tt)
    yy <- center[2] + radius * sin(tt)
    return(tibble(x = xx, y = yy, z = 1:length(xx)))
}



big_land_plotter("outbreak_size", "strong", col_fct = "wt_pp",
                 shp_lty_fct = "wt_vp",
                 y_facet_fct = "sd_N",
                 x_facet_fct = c("virus_attract", "pseudo_repel"),
                 facet_ncol = 4, facet_scales = "free_y")

# When virus starts off Pseudomonas:
# --------------------*

# clustering has much longer tail with a few outbreaks being really large:
big_off_df <- sim_df |>
    filter(wasp_resp == "weak", sd_N == 0, virus_attract == 1, pseudo_repel == 5,
           n_pseudo == 9000, wt_vp == "off *Pseudo.*", outbreaks == "small")

big_off_df |>
    select(wt_pp, sim) |>
    mutate(outbreak_size = map(sim, \(x) x$n_infected[x$n_infected > 1])) |>
    select(-sim) |>
    unnest(outbreak_size) |>
    ggplot(aes(outbreak_size, after_stat(density),
               color = wt_pp, linewidth = wt_pp, linetype = wt_pp)) +
    geom_freqpoly(bins = 25) +
    scale_linetype_manual(values = c(clustered = "22", uniform = "solid")) +
    scale_linewidth_manual(values = c(clustered = 1.25, uniform = 1)) +
    scale_color_manual(values = c(clustered = "dodgerblue", uniform = "gray60"))

# Biggest 9 outbreaks for clustered and uniform landscapes:
big_off_lands <- big_off_df |>
    (\(x) {z <- x$sim; names(z) <- x$wt_pp; return(z)})() |>
    map(\(x) {
        x |>
            arrange(desc(n_infected)) |>
            slice_head(n = 9) |>
            getElement("rep")
    }) |>
    (\(x) {
        dd <- list.files("_scripts/interm-data", "large-plantscapes-.?.?.rds",
                         full.names = TRUE) |>
            map(read_rds) |>
            list_rbind() |>
            filter(wasp_resp == big_off_df$wasp_resp[[1]],
                   sd_N == big_off_df$sd_N[[1]],
                   virus_attract == big_off_df$virus_attract[[1]],
                   pseudo_repel == big_off_df$pseudo_repel[[1]],
                   n_pseudo == big_off_df$n_pseudo[[1]],
                   ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") == big_off_df$wt_vp[[1]],
                   outbreaks == big_off_df$outbreaks[[1]]) |>
            mutate(wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
                       factor(levels = c("uniform", "clustered")))
        imap(x, \(xx, i) {
            dd$landscape[[which(dd$wt_pp == i)]][,,xx]
        })
    })()

# land_plotter(big_off_lands$uniform) +
#     geom_path(data = circle_maker(c(0.5, 0.5), radius = an_environ$radius,
#                                   rad_min = 0, rad_max = 0.5 * pi),
#               aes(x, y), inherit.aes = FALSE, color = "white")
# land_plotter(big_off_lands$clustered) +
#     geom_path(data = circle_maker(c(0.5, 0.5), radius = an_environ$radius,
#                                   rad_min = 0, rad_max = 0.5 * pi),
#               aes(x, y), inherit.aes = FALSE, color = "white")
#
# land_plotter(big_off_lands$uniform[1:8, 1:8,])
# land_plotter(big_off_lands$clustered[1:8, 1:8,])



# When virus starts on Pseudomonas:
# --------------------*

# clustering has much longer tail with a few outbreaks being really large:
big_on_df <- sim_df |>
    filter(wasp_resp == "weak", sd_N == 0, virus_attract == 1, pseudo_repel == 5,
           n_pseudo == 9000, wt_vp == "on *Pseudo.*", outbreaks == "small")

big_on_df |>
    select(wt_pp, sim) |>
    mutate(outbreak_size = map(sim, \(x) x$n_infected[x$n_infected > 1])) |>
    select(-sim) |>
    unnest(outbreak_size) |>
    ggplot(aes(outbreak_size, after_stat(density),
               color = wt_pp, linewidth = wt_pp, linetype = wt_pp)) +
    geom_freqpoly(bins = 25) +
    scale_linetype_manual(values = c(clustered = "22", uniform = "solid")) +
    scale_linewidth_manual(values = c(clustered = 1.25, uniform = 1)) +
    scale_color_manual(values = c(clustered = "dodgerblue", uniform = "gray60")) +
    scale_color_manual(values = c(clustered = "dodgerblue", uniform = "gray60"))

# Landscapes for biggest 2 outbreaks for clustered and uniform landscapes:
big_on_lands <- big_on_df |>
    (\(x) {z <- x$sim; names(z) <- x$wt_pp; return(z)})() |>
    map(\(x) {
        x |>
            arrange(desc(n_infected)) |>
            slice_head(n = 2) |>
            getElement("rep")
    }) |>
    (\(x) {
        dd <- list.files("_scripts/interm-data", "large-plantscapes-.?.?.rds",
                         full.names = TRUE) |>
            map(read_rds) |>
            list_rbind() |>
            filter(wasp_resp == big_on_df$wasp_resp[[1]],
                   sd_N == big_on_df$sd_N[[1]],
                   virus_attract == big_on_df$virus_attract[[1]],
                   pseudo_repel == big_on_df$pseudo_repel[[1]],
                   n_pseudo == big_on_df$n_pseudo[[1]],
                   ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") == big_on_df$wt_vp[[1]],
                   outbreaks == big_on_df$outbreaks[[1]]) |>
            mutate(wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
                       factor(levels = c("uniform", "clustered")))
        lands <- imap(x, \(xx, i) {
            dd$landscape[[which(dd$wt_pp == i)]][,,xx]
        })
    })()



# land_plotter(big_on_lands$uniform) +
#     geom_path(data = circle_maker(c(0.5, 0.5), radius = an_environ$radius,
#                                   rad_min = 0, rad_max = 0.5 * pi),
#               aes(x, y), inherit.aes = FALSE, color = "white")
# land_plotter(big_on_lands$clustered) +
#     geom_path(data = circle_maker(c(0.5, 0.5), radius = an_environ$radius,
#                                   rad_min = 0, rad_max = 0.5 * pi),
#               aes(x, y), inherit.aes = FALSE, color = "white")
#
# land_plotter(big_on_lands$uniform[1:8, 1:8,])
# land_plotter(big_on_lands$clustered[1:8, 1:8,])



