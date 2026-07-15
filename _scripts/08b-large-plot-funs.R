
suppressPackageStartupMessages({
    library(grid)
    library(gridtext)
})

#
# Add np = 0 for each landscape type.
# Doing this to dataframes below allows filtering by wt_vp and wt_pp and
# also retrieving n_pseudo == 0
#
add_no_pseudo_points <- function(data_df) {
    stopifnot(all(c("wt_vp", "wt_pp") %in% colnames(data_df)))
    split_cols <- c("wasp_resp", "p_load", "sd_N", "virus_attract",
                    "pseudo_repel", "fly_p") |>
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






#
# Create color scale and geom objects for `baseline_plotter`.
# Used inside `make_baseline_objs`.
#
make_baseline_col_geom <- function(def_vals,
                                   col_fct,
                                   lty_fct,
                                   color_vals,
                                   filt_df,
                                   inters,
                                   incl_vals,
                                   multiline_col_title,
                                   env) {

    if (!is.null(lty_fct)) stopifnot(length(lty_fct) == 1 && lty_fct %in% colnames(filt_df))
    if (!is.null(col_fct)) stopifnot(all(col_fct %in% colnames(filt_df)))

    title_maker <- function(varbl) {
        varbl |>
            pretty_params(cap1 = TRUE) |>
            str_replace_all("Pseudomonas", "Pseudo.") |>
            (\(x) {
                if (!multiline_col_title) return(x)
                xs <- str_split(x, "\\(")[[1]]
                xs[1] <- str_replace_all(xs[1], " ", "<br>")
                return(str_c(xs, collapse = "("))
            })()
    }
    add_lty <- function(line_geom, lty_fct) {
        if (!is.null(lty_fct)) {
            line_geom[[1]][["mapping"]][["linetype"]] <- quo(.data[[lty_fct]])
            if (!is.null(col_fct)) {
                line_geom[[2]][["mapping"]][["group"]] <- quo(interaction(.data[[lty_fct]],
                                                                          .data[[col_fct]]))
            } else {
                line_geom[[2]][["mapping"]][["group"]] <- quo(.data[[lty_fct]])
            }
            line_geom[[3]] <- scale_linetype(title_maker(lty_fct))
        }
        return(line_geom)
    }


    col_scale <- NULL
    line_geom <- list(geom_line(linewidth = 0.5),
                      geom_ribbon(aes(ymin = lower, ymax = upper),
                                  fill = "black", alpha = 0.25, color = NA))
    if (is.null(col_fct)) {
        line_geom <- add_lty(line_geom, lty_fct)
        assign("col_scale", col_scale, envir = env)
        assign("line_geom", line_geom, envir = env)
        return(filt_df)
    }


    if (length(col_fct) == 1L) {

        if (!is.factor(filt_df[[col_fct]]))
            filt_df[[col_fct]] <- factor(filt_df[[col_fct]])
        col_title <- title_maker(col_fct)
        if (!is.null(color_vals)) {
            stopifnot(length(color_vals) == length(levels(filt_df[[col_fct]])))
            col_scale <- scale_color_manual(col_title, values = color_vals,
                                            aesthetics = c("colour","fill"))
        } else {
            col_scale <- scale_color_scico_d(col_title, end = 0.8,
                                             palette = "hawaii",
                                             aesthetics = c("colour","fill"))
        }
        line_geom <- list(geom_line(aes(color = .data[[col_fct]]), linewidth = 0.5),
                          geom_ribbon(aes(ymin = lower, ymax = upper,
                                          fill = .data[[col_fct]]),
                                      alpha = 0.25, color = NA))

    } else {

        if (!all(col_fct %in% names(def_vals))) {
            stop("col_fct must be one of ", paste(names(def_vals),
                                                  collapse = ", "))
        }

        if (inters) {

            col_title <- "Param(s).<br>that<br>differ:"
            if (!multiline_col_title) col_title <- "Param(s). that differ:"
            for (x in col_fct) {
                if (length(unique(filt_df[[x]])) > 2L) stop(x, " has > 2 levels")
                if (!is.factor(filt_df[[x]])) filt_df[[x]] <- factor(filt_df[[x]])
                if (! def_vals[[x]] %in% levels(filt_df[[x]])) {
                    stop("levels for ", x, " should contain '",
                         def_vals[[x]], "'")
                }
            }
            filt_df[["col_fct"]] <- map_chr(1:nrow(filt_df), \(i) {
                z <- keep(col_fct, \(x) filt_df[[x]][[i]] != def_vals[[x]])
                if (length(z) == 0) return("none")
                paste(pretty_params(z, TRUE, serif = TRUE), collapse = "&amp;")
            })

        } else {

            col_title <- "Differences:"

            n_diff_lvls <- rep(0L, nrow(filt_df))
            for (x in col_fct) {
                if (length(unique(filt_df[[x]])) > 2L) stop(x, " has > 2 levels")
                if (!is.factor(filt_df[[x]])) filt_df[[x]] <- factor(filt_df[[x]])
                if (! def_vals[[x]] %in% levels(filt_df[[x]])) {
                    stop("levels for ", x, " should contain '",
                         def_vals[[x]], "'")
                }
                n_diff_lvls <- n_diff_lvls + as.integer(filt_df[[x]] != def_vals[[x]])
            }
            filt_df <- filt_df[n_diff_lvls <= 1L,]
            filt_df[["col_fct"]] <- map_chr(1:nrow(filt_df), \(i) {
                z <- keep(col_fct, \(x) filt_df[[x]][[i]] != def_vals[[x]])
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
                              values = c(0.75, rep(0.5, length(levels(filt_df[["col_fct"]]))-1)),
                              guide = "none"))

        col_fct <- "col_fct"

    }

    line_geom <- add_lty(line_geom, lty_fct)

    assign("col_scale", col_scale, envir = env)
    assign("line_geom", line_geom, envir = env)
    assign("col_fct", col_fct, envir = env)

    return(filt_df)

}


#
# Do background creating objects for baseline plotter
#
make_baseline_objs <- function(data_df,
                               obs_breaks,
                               outcomes,
                               col_fct,
                               lty_fct,
                               color_vals,
                               inters,
                               incl_vals,
                               obs_max,
                               p_title,
                               p_subtitle,
                               p_tag,
                               non_defaults,
                               multiline_col_title,
                               env) {

    if (is.null(obs_breaks)) obs_breaks <- scales::breaks_extended(n = 3)
    stopifnot(length(outcomes) == 1L && outcomes %in%
                  c("all", "both", "p_emerge", "outbreak_size", "n_infected"))
    if (outcomes == "all") outcomes <- c("p_emerge", "outbreak_size", "n_infected")
    else if (outcomes == "both") outcomes <- c("p_emerge", "outbreak_size")
    stopifnot(all(c("wasp_resp", "n_pseudo", outcomes, "boots") %in% colnames(data_df)))

    def_vals <- list(p_load = 0.5,
                     fly_p = 0.1,
                     sd_N = 0,
                     virus_attract = 1,
                     pseudo_repel = 1,
                     wt_vp = "off *Pseudo.*",
                     wt_pp = "uniform")
    # Replace with non-defaults for filtering:
    if (length(non_defaults) > 0) {
        stopifnot(is.list(non_defaults) && !is.null(names(non_defaults)))
        stopifnot(!any(names(non_defaults) %in% col_fct))
        stopifnot(!any(names(non_defaults) %in% lty_fct))
        stopifnot(all(names(non_defaults) %in% names(def_vals)))
        def_vals[names(non_defaults)] <- non_defaults
    }

    filt_df <- data_df
    for (x in names(def_vals)) {
        if (x %in% col_fct || x %in% lty_fct) next
        if (x %in% colnames(filt_df)) {
            stopifnot(def_vals[[x]] %in% filt_df[[x]])
            filt_df <- filt_df[filt_df[[x]] == def_vals[[x]], ]
        }
    }

    # Plot panel rows, columns, respectively:
    panel_dims <- c(length(outcomes), length(unique(filt_df$wasp_resp)))
    if (!panel_dims[1] %in% 1:3){
        stop("\n# of unique outcomes must be 1, 2, or 3. You have ",
             panel_dims[1], ".")
    }
    if (!panel_dims[2] %in% 1:2){
        stop("\n# of unique wasp_resp levels must be 1 or 2. You have ",
             panel_dims[2], ".")
    }
    dsn <- map(1:panel_dims[1], \(i)
               paste(LETTERS[((i-1L)*panel_dims[2]+1L):(i*panel_dims[2])],
                     collapse = "#")) |>
        c(sep = "\n") |> do.call(what = paste)
    dsn_widths <- head(rep(c(1, 0.05), panel_dims[2]), -1)
    patchwork_layout <- plot_layout(design = dsn, widths = dsn_widths,
                                    guides = "collect", axis_titles = "collect")

    if (length(p_title) == 1) p_title <- rep(p_title, prod(panel_dims))
    stopifnot(length(p_title) == prod(panel_dims))
    if (length(p_subtitle) == 1) p_subtitle <- rep(p_subtitle, prod(panel_dims))
    stopifnot(length(p_subtitle) == prod(panel_dims))
    if (length(p_tag) == 1) p_tag <- rep(p_tag, prod(panel_dims))
    if (length(p_tag) != prod(panel_dims)) {
        cat("length(p_tag) = ", length(p_tag), "\n")
        cat("prod(panel_dims) = ", prod(panel_dims), "\n")
    }
    stopifnot(length(p_tag) == prod(panel_dims))


    assign("obs_breaks", obs_breaks, envir = env)
    assign("outcomes", outcomes, envir = env)
    assign("panel_dims", panel_dims, envir = env)
    assign("patchwork_layout", patchwork_layout, envir = env)

    filt_df <- make_baseline_col_geom(def_vals, col_fct, lty_fct, color_vals, filt_df,
                                      inters, incl_vals, multiline_col_title,
                                      env)

    filt_df <- filt_df |>
        select(wasp_resp, n_pseudo, all_of(env$col_fct), all_of(env$lty_fct),
               all_of(outcomes), boots) |>
        pivot_longer(all_of(outcomes), names_to = "outcome") |>
        mutate(outcome = factor(outcome, levels = c("p_emerge", "outbreak_size",
                                                    "n_infected") |>
                                    keep(\(x) x %in% outcomes))) |>
        mutate(lower = map2_dbl(boots, outcome, \(b, o) b[[paste(o)]][["Lower"]]),
               upper = map2_dbl(boots, outcome, \(b, o) b[[paste(o)]][["Upper"]]))

    if (any(outcomes != "p_emerge")) {
        max_data_obs <- max(filt_df$value, filt_df$upper)
        if (!is.null(obs_max)) {
            if (max_data_obs > obs_max) {
                stop("`obs_max` too low. Max value in data is ", max_data_obs)
            }
        } else {
            obs_max <- max_data_obs
        }
    } else {
        obs_max <- NA_real_
    }

    filt_df_split <- filt_df |>
        split(~ wasp_resp + outcome, drop = TRUE)

    names(p_title) <- names(filt_df_split)
    names(p_subtitle) <- names(filt_df_split)
    names(p_tag) <- names(filt_df_split)

    assign("obs_max", obs_max, envir = env)
    assign("filt_df_split", filt_df_split, envir = env)
    assign("p_title", p_title, envir = env)
    assign("p_subtitle", p_subtitle, envir = env)
    assign("p_tag", p_tag, envir = env)


    invisible(NULL)
}





baseline_plotter <- function(col_fct = NULL,
                             lty_fct = NULL,
                             outcomes = "all",
                             color_vals = NULL,
                             inters = FALSE,
                             fix_y = TRUE,
                             obs_breaks = NULL,
                             obs_max = NULL,
                             incl_vals = FALSE,
                             return_list = FALSE,
                             p_title = list(waiver()),
                             p_subtitle = list(waiver()),
                             p_tag = list(waiver()),
                             non_defaults = list(),
                             multiline_ylab = FALSE,
                             multiline_col_title = TRUE,
                             data_df = sim_df) {

    # col_fct = NULL
    # lty_fct = NULL
    # # col_fct = c("sd_N", "wt_pp")
    # outcomes = "n_infected"
    # color_vals = NULL
    # # color_vals = c("black", "dodgerblue", "firebrick2")
    # inters = FALSE
    # fix_y = TRUE; obs_breaks = NULL; obs_max = NULL; incl_vals = TRUE; return_list = FALSE
    # p_title = list(waiver()); p_subtitle = list(waiver()); p_tag = list(waiver())
    # non_defaults = list(); multiline_ylab = TRUE
    # data_df = sim_df
    #
    # rm(col_fct, outcomes, color_vals, inters, fix_y, obs_breaks, obs_max)
    # rm(incl_vals, return_list, p_title, p_subtitle, p_tag, non_defaults, multiline_ylab, data_df)
    # rm(filt_df_split, panel_dims, patchwork_layout, col_scale, line_geom)
    # rm(p_list)


    if (!is.null(obs_max) && !fix_y && outcomes %in% c("p_emerge", "both", "all")) {
        stop("Cannot specify `obs_max` if `fix_y = FALSE`")
    }
    # Creates objects `filt_df_split`, `panel_dims`, `patchwork_layout`, `col_scale`,
    # and `line_geom`
    #
    # Edits objects `obs_breaks`, `outcomes`, `p_title`, `p_subtitle`, `p_tag`,
    # `obs_max`, and (optionally) `col_fct`
    make_baseline_objs(data_df, obs_breaks, outcomes, col_fct, lty_fct, color_vals,
                       inters, incl_vals, obs_max, p_title, p_subtitle, p_tag,
                       non_defaults, multiline_col_title, environment())

    p_list <- filt_df_split |>
        imap(\(dd, n) {
            yvar <- paste(dd$outcome[[1]])
            y_lines <- switch(yvar, p_emerge = c(0, 1),
                              outbreak_size = 2, n_infected = 1)
            y_lims <- switch(yvar, p_emerge = c(0, 1.1),
                             outbreak_size = c(2, obs_max),
                             n_infected = c(1, obs_max))
            free_y <- !fix_y && yvar != "p_emerge"
            if (free_y) y_lims[2] <- NA
            yb <- switch(yvar, p_emerge = 0:2/2,
                         outbreak_size = obs_breaks,
                         n_infected = obs_breaks)
            y_lab <- yvar_desc[[yvar]] |>
                (\(x) ifelse(yvar != "p_emerge", paste("mean", x), x))() |>
                first_cap()
            if (multiline_ylab) y_lab <- str_replace_all(y_lab, " ", "<br>")
            p <- dd |>
                ggplot(aes(n_pseudo / 10e3 * 100, value)) +
                geom_hline(yintercept = y_lines, color = "gray70") +
                line_geom +
                labs(x = "Percent *Pseudomonas* plants", y = y_lab,
                     title = p_title[[n]], subtitle = p_subtitle[[n]],
                     tag = p_tag[[n]]) +
                scale_y_continuous(breaks = yb) +
                scale_x_continuous(breaks = n_pseudo_lvls / 10e3 * 100) +
                col_scale +
                coord_cartesian(ylim = y_lims) +
                theme(plot.margin = margin(0,0,0,0))
            if (fix_y && panel_dims[2] > 1L &&
                dd$wasp_resp[[1]] == levels(wasp_resp_fct)[[2]]) {
                p <- p + theme(axis.text.y = element_blank())
            }
            if (panel_dims[1] > 1L && yvar != tail(outcomes, 1)) {
                p <- p + theme(axis.text.x = element_blank())
            }
            return(p)
        })

    if (return_list) return(p_list)

    if (prod(panel_dims) == 1) return(p_list[[1]])

    p_list |>
        wrap_plots() +
        patchwork_layout

}





add_top_labels <- function(plot_list, add_bot_labs = TRUE) {
    top_labs <- levels(wasp_resp_fct) |>
        map(\(x) {
            grob <- richtext_grob(scenario_title(x),
                                  gp = gpar(fontsize = 13, lineheight = 0.8))
            wrap_elements(panel = grob)
        })
    bot_labs <- map(rep(c("uninhabited", "*Pseudomonas*"), 2),
                    \(x) {
                        txt <- paste0("Virus starts on<br>", x, "<br>plant")
                        grob <- richtext_grob(txt, gp = gpar(fontsize = 10,
                                                             lineheight = 0.8))
                        wrap_elements(panel = grob)
                    })
    if (add_bot_labs) out <- c(plot_list, top_labs, bot_labs)
    else out <- c(plot_list, top_labs)
    return(out)
}


