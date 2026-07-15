#'
#' Heatmaps and line graphs of outbreak sizes for simulations of two
#' levels of wasp responsiveness to aphid densities:
#' "weak" (so Pseudomonas increases outbreak size)
#' "strong" (so Pseudomonas decreases outbreak size)
#' ... where we manipulate one or two parameters at a time.
#'


source("_scripts/00-preamble.R")

.overwrite <- FALSE




# ============================================================================*
# ============================================================================*

# >> 1-par manips ====

# ============================================================================*
# ============================================================================*


# --------------------------------------*
# _ sims ----
# --------------------------------------*

manip_pars <- list(Y0 = seq(1, 10, 0.5),
                   N0 = 1:20 * 10,
                   sd_N = 0:20 * 5,
                   K = 5:25 * 1000,
                   virus_attract = 0:19 / 2 + 1,
                   pseudo_repel = 0:19 / 2 + 1,
                   pseudo_surv = seq(0.83, 1, 0.01),
                   zeta = seq(0, 1, 0.05)) |>
    map(\(x) round(x, 2))



one_manip_sim <- function(wasp_resp, x_name, x_val) {

    # wasp_resp = "strong"
    # x_name = "Y0"
    # x_val = manip_pars[[x_name]][[1]]
    # rm(wasp_resp, x_name, x_val, args, sims_p, sims_np, out)

    args <- list(wasp_resp = wasp_resp, n_pseudo = 3L, large_sims = TRUE)
    args[[x_name]] <- x_val

    sims_p <- do.call(run_sim_combos, args) |>
        select(-rep) |>
        mutate(n_pseudo = 3L)
    sims_np <- do.call(run_sim_combos, list_assign(args, n_pseudo = 0L)) |>
        select(-rep) |>
        mutate(n_pseudo = 0L)

    out <- bind_rows(sims_p, sims_np) |>
        mutate(wasp_resp = .env$wasp_resp,
               par_name = factor(x_name, levels = names(manip_pars)),
               par_val = x_val) |>
        select(wasp_resp, n_pseudo, par_name, par_val, everything())

    return(out)

}


if (.overwrite || !file.exists(interm_files$extreme_manip)) {

    # Takes ~6 min
    set.seed(1531497906)
    manip_sims <- levels(wasp_resp_fct) |>
        map(\(wasp_resp) {
            manip_pars |>
                (\(x) {
                    out <- tibble(vals = do.call(c, manip_pars) |> unname())
                    out$par_name <- manip_pars |>
                        imap(\(x, n) rep(n, length(x))) |>
                        list_c()
                    return(out)
                })() |>
                pmap(\(vals, par_name) {
                    map(vals, \(v) one_manip_sim(wasp_resp, par_name, v)) |>
                        list_rbind()
                }, .progress = .prog_args) |>
                list_rbind()
        }) |>
        list_rbind() |>
        mutate(wasp_resp = factor(wasp_resp, levels = levels(wasp_resp_fct)))

    write_rds(manip_sims, interm_files$extreme_manip, compress = "gz")

} else {

    manip_sims <- read_rds(interm_files$extreme_manip)

}






# --------------------------------------*
# _ plotter ----
# --------------------------------------*


one_manip_plotter <- function(wasp_resp, par_name,
                              outcomes = c("p_emerge", "outbreak_size", "n_infected"),
                              .og_vals = TRUE,
                              .add_pts = FALSE,
                              .return_list = FALSE,
                              .x_pos = "bottom",
                              .seed = NULL) {
    # wasp_resp = "weak"; par_name = "zeta"; .og_vals = TRUE; .add_pts = FALSE
    # outcomes = c("p_emerge", "outbreak_size", "n_infected")
    # .return_list = FALSE; .x_pos = "bottom"
    # rm(wasp_resp, par_name, outcomes, out_lvls, .og_vals, .add_pts, .return_list)
    # rm(.x_pos, dd, dd_og, .points_lines, plot_list, x_lab)

    if (!is.null(.seed)) set.seed(.seed)

    out_lvls <- c("p_emerge", "outbreak_size", "n_infected") |>
        keep(\(x) x %in% outcomes)

    dd <- manip_sims |>
        filter(wasp_resp %in% .env$wasp_resp, par_name == .env$par_name) |>
        mutate(n_pseudo = factor(n_pseudo)) |>
        group_by(n_pseudo, par_val) |>
        summarize(boots = list(ci_booter(n_infected, outcomes)),
                  p_emerge = mean(n_infected > 1),
                  outbreak_size = mean(n_infected[n_infected > 1]),
                  n_infected = mean(n_infected),
                  .groups = "drop") |>
        pivot_longer(p_emerge:n_infected, names_to = "outcome") |>
        filter(outcome %in% outcomes) |>
        mutate(outcome = factor(outcome, levels = out_lvls),
               lower = map2_dbl(boots, outcome, \(.b, .o) .b[[paste(.o)]][["Lower"]]),
               upper = map2_dbl(boots, outcome, \(.b, .o) .b[[paste(.o)]][["Upper"]])) |>
        select(-boots)
    if (par_name == "K") dd$par_val <- carrying_capacity(K = dd$par_val, pred_surv = 1)

    # Get original value(s)
    dd_og <- map_dbl(wasp_resp, \(wr) {
        og <- run_sim_combos(wasp_resp = wr, n_pseudo = 0, return_args = TRUE) |>
            getElement(par_name)
        if (is.null(og)) og <- formals(lil_plantscape)[[par_name]]
        if (is.null(og)) og <- eval(formals(sim_plantscape)[[par_name]])
        if (is.null(og)) og <- eval(formals(make_insects_ptr)[[par_name]])
        if (is.null(og)) og <- eval(formals(make_disease_ptr)[[par_name]])
        if (par_name == "K") og <- carrying_capacity(K = og, pred_surv = 1)
        return(og)
    })

    .points_lines <- list(geom_point(), geom_line(linewidth = 0.5))
    if (!.add_pts) {
        .points_lines <- .points_lines[2]
        .points_lines[[1]][["aes_params"]][["linewidth"]] <- 1
    }

    if (par_name == "K") {
        x_lab <- paste("Aphid carrying capacity",
                       "(<span style=\"font-family: serif\"><i>K</i>",
                       "(<i>&lambda;</i> &minus; 1)</span>)")
    } else x_lab <- pretty_params(par_name) |> first_cap()


    plot_list <- dd |>
        filter(!is.na(value)) |>
        split(~ outcome) |>
        map(\(ddd) {
            yvar <- paste(ddd$outcome[[1]])
            yl <- switch(yvar, p_emerge = c(0, 1), outbreak_size = c(2, 9),
                         n_infected = c(1,9))
            yb <- switch(yvar, p_emerge = 0:2/2, outbreak_size = c(3,6,9),
                         n_infected = c(1,5,9))
            # yp <- switch(yvar, p_emerge = "left", outbreak_size = "right")
            y_lab <- yvar_desc[[yvar]] |>
                (\(x) ifelse(yvar == "p_emerge", x, paste("mean", x)))() |>
                first_cap()
            p <- ddd |>
                ggplot(aes(par_val, value, color = n_pseudo)) +
                geom_hline(yintercept = yl, color = "gray70")
            if (.og_vals) p <- p + geom_vline(xintercept = dd_og, color = "gray70",
                                              linetype = "33")
            p +
                .points_lines +
                geom_ribbon(aes(ymin = lower, ymax = upper, fill = n_pseudo),
                            alpha = 0.25, color = NA) +
                labs(x = x_lab, y = y_lab) +
                scale_y_continuous(breaks = yb) +
                # scale_y_continuous(breaks = yb, position = yp) +
                scale_x_continuous(position = .x_pos) +
                scale_color_manual("*Pseudo.*<br>plants", values = np_pal,
                                   aesthetics = c("color","fill")) +
                coord_cartesian(xlim = range(dd$par_val)) +
                theme(plot.margin = margin(0,0,0,0),
                      axis.title.x.top = element_markdown(color = "black", size = 11),
                      axis.text.x.top = element_markdown(color = "black", size = 9))
        })

    if (.return_list) return(plot_list)

    if (length(plot_list) == 1) return(plot_list[[1]])

    wrap_plots(plot_list) +
        plot_layout(ncol = 1, guides = "collect", axes = "collect")
}




# one_manip_plotter("weak", "zeta", c("n_infected"))
# one_manip_plotter(c("weak", "strong"), "zeta")

# --------------------------------------------*
# _ supp. plots ----
# --------------------------------------------*


# Line plot for a parameter, with fill split into prob. emergence and outbreak size
# If parameter is zeta, it doesn't split by wasp_resp bc that only differs by zeta.
supp_line_plotter <- function(par_name, seeds, .incl_ninf = FALSE) {
    stopifnot(length(par_name) == 1 && is.character(par_name))
    stopifnot(length(.incl_ninf) == 1 && is.logical(.incl_ninf))
    stopifnot(is.numeric(seeds))
    if (par_name != "zeta") {
        stopifnot(length(seeds) == 2)
        .outcomes <- c("p_emerge", "outbreak_size")
        if (.incl_ninf) .outcomes <- c(.outcomes, "n_infected")
        no <- length(.outcomes)

        pl <- map2(levels(wasp_resp_fct), seeds,
                   \(wr, seed) {
                       pl <- one_manip_plotter(wr, par_name,
                                               outcomes = .outcomes,
                                               .return_list = TRUE, .seed = seed)
                       pl[[1]] <- pl[[1]] + labs(title = scenario_title(wr))
                       pl[[1]] <- pl[[1]] + theme(axis.text.x = element_blank())
                       if (wr == tail(levels(wasp_resp_fct),1)) {
                           pl <- map(pl, \(x) x + theme(axis.text.y=element_blank()))
                       }
                       return(pl)
                   }) |>
            do.call(what = c) |>
            (\(x) {
                if (.incl_ninf) {
                    z <- list(x[[1]], x[[4]], plot_spacer(), plot_spacer(),
                              x[[2]], x[[5]], plot_spacer(), plot_spacer(),
                              x[[3]], x[[6]])
                } else {
                    z <- list(x[[1]], x[[3]], plot_spacer(), plot_spacer(),
                              x[[2]], x[[4]])
                }
                return(z)
            })()
        hts <- c(1, 0.05, 1)
    } else {
        stopifnot(length(seeds) == 1)
        pl <- one_manip_plotter(levels(wasp_resp_fct), par_name,
                                outcomes = .outcomes,
                                .return_list = TRUE, .seed = seeds) |>
            (\(x) {
                x[[1]] <- x[[1]] + theme(axis.text.x = element_blank())
                #
                if (.incl_ninf) {
                    z <- list(plot_spacer(), x[[1]], plot_spacer(), x[[2]],
                              plot_spacer(), x[[3]])
                } else {
                    z <- list(plot_spacer(), x[[1]], plot_spacer(), x[[2]])
                }
                return(z)
            })()
        hts <- c(0.05, 1, 0.1, 1)
    }
    if (.incl_ninf) hts <- c(hts, tail(hts, 2))
    pl |>
        wrap_plots(nrow = length(hts), guides = "collect", axis_titles = "collect",
                   heights = hts) +
        plot_annotation(tag_level = "A") &
        theme(plot.tag = element_markdown(face = "bold"),
              plot.tag.location = "panel",
              plot.tag.position = c(0.05, 1.05))
}


# supp_line_plotter("K", c(1609768752, 1975481712))
# supp_line_plotter("pseudo_surv", c(1546463762, 915602074), .incl_ninf = TRUE)



if (.overwrite) {
    save_plot("_plots/extremes-manip-lines-pseudo_surv-all-outcomes.pdf",
              supp_line_plotter("pseudo_surv", c(1546463762, 915602074)),
              width = 6.5, height = 5)
    # Because K isn't shown in main text, we also want to plot n_infected
    save_plot("_plots/extremes-manip-lines-K-all-outcomes.pdf",
              supp_line_plotter("K", c(1609768752, 1975481712), .incl_ninf = TRUE),
              width = 6.5, height = 7)
    # Because the `wasp_resp` only differs by zeta, we plot that differently:
    save_plot("_plots/extremes-manip-lines-zeta-all-outcomes.pdf",
              supp_line_plotter("zeta", 65130342),
              width = 4, height = 5)
}





# ============================================================================*
# ============================================================================*

# >> 2-par manips ====

# ============================================================================*
# ============================================================================*


# --------------------------------------*
# _ sims ----
# --------------------------------------*




one_manip2_sim <- function(wasp_resp, par_name_a, par_val_a, par_name_b, par_val_b) {

    # wasp_resp = "strong"
    # par_name_a = "Y0"; par_val_a = manip_pars[[par_name_a]][[1]]
    # par_name_b = "N0"; par_val_b = manip_pars[[par_name_b]][[10]]
    # rm(wasp_resp, par_name_a, par_val_a, par_name_b, par_val_b, args, out)
    #

    # In case these are factors:
    par_name_a <- paste(par_name_a)
    par_name_b <- paste(par_name_b)

    args <- list(wasp_resp = wasp_resp, large_sims = TRUE)
    args[[par_name_a]] <- par_val_a
    args[[par_name_b]] <- par_val_b

    out <- c(3L, 0L) |>
        map(\(np) {
            .np <- as.integer(np)
            do.call(run_sim_combos, list_assign(args, n_pseudo = .np)) |>
                select(-rep) |>
                mutate(p_emerge = as.integer(n_infected > 1),
                       outbreak_size = ifelse(n_infected > 1, n_infected, NA)) |>
                summarize(boots = list(ci_booter(n_infected, "all")),
                          across(-boots, \(x) mean(x, na.rm = TRUE))) |>
                mutate(n_pseudo = .np) |>
                relocate(boots, .after = last_col())
        }) |>
        list_rbind() |>
        mutate(wasp_resp = .env$wasp_resp,
               par_name_a = factor(.env$par_name_a, levels = names(manip_pars)),
               par_name_b = factor(.env$par_name_b, levels = names(manip_pars)),
               par_val_a = .env$par_val_a,
               par_val_b = .env$par_val_b) |>
        select(wasp_resp, n_pseudo, starts_with("par_"), everything())

    return(out)

}




if (!file.exists(interm_files$extreme_manip2) || .overwrite) {

    # Takes ~15 min
    set.seed(2025929231)
    manip2_sims <- tibble(par_name_a = "Y0",
                          par_name_b = "N0") |>
        mutate(wasp_resp = list(levels(wasp_resp_fct))) |>
        unnest(wasp_resp) |>
        mutate(vals = map2(par_name_a, par_name_b,
                           \(par_name_a, par_name_b) {
                               crossing(par_val_a = manip_pars[[par_name_a]],
                                        par_val_b = manip_pars[[par_name_b]])
                           })) |>
        unnest(vals) |>
        arrange(wasp_resp, par_name_a, par_name_b, par_val_a, par_val_b) |>
        select(wasp_resp, par_name_a, par_val_a, par_name_b, par_val_b) |>
        pmap(one_manip2_sim, .progress = .prog_args) |>
        list_rbind() |>
        mutate(wasp_resp = factor(wasp_resp, levels = levels(wasp_resp_fct))) |>
        rename(Y0 = par_val_a, N0 = par_val_b) |>
        select(-starts_with("par_name_"))

    write_rds(manip2_sims, interm_files$extreme_manip2, compress = "gz")

} else {

    manip2_sims <- read_rds(interm_files$extreme_manip2)

}






#
# _ N0 / Y0 thresholds ----
#

# ... for when outbreaks never occur
manip2_sims |>
    group_by(wasp_resp, N0) |>
    filter(Y0 == min(Y0[p_emerge == 0])) |>
    group_by(wasp_resp) |>
    summarize(rel = median(Y0 / N0))
manip2_sims |>
    group_by(wasp_resp, N0) |>
    filter(Y0 == min(Y0[p_emerge == 0])) |>
    ungroup() |>
    ggplot(aes(Y0, N0)) +
    geom_point() +
    facet_wrap(~ wasp_resp, nrow = 1)

# ... for when outbreaks always occur
manip2_sims |>
    filter(p_emerge == 1) |>
    group_by(wasp_resp, N0) |>
    filter(Y0 == max(Y0)) |>
    group_by(wasp_resp) |>
    summarize(rel = median(Y0 / N0))
manip2_sims |>
    filter(p_emerge == 1) |>
    group_by(wasp_resp, N0) |>
    filter(Y0 == max(Y0)) |>
    ungroup() |>
    ggplot(aes(Y0, N0)) +
    geom_point() +
    facet_wrap(~ wasp_resp, nrow = 1)



#
# _ above / below thresholds ----
#
set.seed(778424076)
thresh_sims <- crossing(wasp_resp = sort(unique(manip_sims$wasp_resp)),
                        n_pseudo = c(0L, 3L),
                        ny = factor(0:1, labels = c("below", "above"))) |>
    pmap(\(wasp_resp, n_pseudo, large_sims, ny) {
        Y0  <- ifelse(ny == "below", 3, 1)
        N0 <- ifelse(ny == "below", 50, 150)
        out <- run_sim_combos(wasp_resp = wasp_resp, n_pseudo = n_pseudo,
                              Y0 = Y0, N0 = N0,
                              large_sims = TRUE, summ = "time",
                              out_stages = "two") |>
            mutate(aphids = aphids_juv + aphids_adu + alates_juv + alates_adu + parasitized,
                   alates = alates_adu) |>
            group_by(time) |>
            summarize(across(c(aphids, alates, wasps, virus), mean),
                      .groups = "drop") |>
            mutate(wasp_resp = .env$wasp_resp,
                   n_pseudo = .env$n_pseudo,
                   thresh = ny) |>
            select(thresh, wasp_resp, n_pseudo, everything())
    }) |>
    list_rbind()


thresh_empty_data <- thresh_sims |>
    select(aphids:virus) |>
    pivot_longer(aphids:virus, names_to = "species", values_to = "density") |>
    mutate(species = factor(species, levels =
                                c("aphids", "alates", "wasps", "virus"))) |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 1)

thresh_p <- thresh_sims |>
    distinct(thresh, wasp_resp) |>
    arrange(thresh, wasp_resp) |>
    pmap(\(wasp_resp, thresh) {
        # wasp_resp = "strong"; thresh = "above"; .title = NULL
        # rm(wasp_resp, thresh, .title, .subtitle, .strip)
        wasp_resp <- paste(wasp_resp)
        thresh <- paste(thresh)
        if (thresh == "below") {
            .title <- scenario_title(wasp_resp)
        } else .title <- waiver()
        .subtitle <- serify("", "*Y*<sub>0</sub> / *N*<sub>0</sub>",
                           list(below = " too high", above = " too low")[[thresh]])
        # if (wasp_resp == levels(wasp_resp_fct)[1]) {
        if (wasp_resp == levels(wasp_resp_fct)[2]) {
            .strip <- element_blank()
        } else .strip <- element_markdown(size = 9, angle = 0, hjust = 1)
        thresh_sims |>
            filter(wasp_resp == .env$wasp_resp, thresh == .env$thresh) |>
            select(-wasp_resp, -thresh) |>
            pivot_longer(aphids:virus, names_to = "species",
                         values_to = "density") |>
            mutate(n_pseudo = factor(n_pseudo),
                   species = factor(species, levels =
                                        c("aphids", "alates", "wasps", "virus"))) |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_hline(data = filter(thresh_empty_data, species == "virus"),
                       aes(yintercept = density), color = "gray70") +
            geom_blank(data = thresh_empty_data) +
            geom_line(aes(color = n_pseudo), linewidth = 1) +
            scale_color_manual("*Pseudo.*<br>plants", values = np_pal) +
            # scale_linetype_manual("*Pseudo.*<br>plants", values = np_linetypes) +
            # scale_linewidth_manual("*Pseudo.*<br>plants", values = np_linewidths) +
            labs(x = "Time (days)", y = "Density", title = .title,
                 subtitle = .subtitle) +
            scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
            coord_cartesian(clip = "off") +
            facet_wrap(~ species, ncol = 1, scales = "free_y",
                       strip.position = "left") +
            theme(strip.text.y = .strip,
                  strip.text.y.left = .strip,
                  strip.placement = "outside",
                  axis.text.x = element_markdown(size = 7),
                  axis.text.y = element_markdown(size = 7),
                  plot.subtitle = element_markdown())
    }) |>
    wrap_plots(ncol = 2, guides = "collect", axes = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"),
          plot.tag.location = c("panel", "plot", "margin")[1],
          plot.tag.position = c(0.05, 0.95))

# thresh_p

if (.overwrite) {
    save_plot("_plots/extreme-manips-thresholds.pdf", thresh_p, width = 6.5, height = 6)
}








# --------------------------------------*
# _ heatmap fxns ----
# --------------------------------------*


source("_scripts/02b-heatmap-funs.R")



# #
# # Minimum and maximums for the effect of Pseudomonas on each measure.
# # Used for in `pseudo_eff_heatmap` above.
# #
# manip2_sims |>
#     select(wasp_resp, Y0, N0, n_pseudo, n_infected:outbreak_size) |>
#     group_by(wasp_resp, Y0, N0) |>
#     summarize(across(n_infected:outbreak_size,
#                      \(x) x[n_pseudo != "0"] - x[n_pseudo == "0"]),
#               .groups = "drop") |>
#     summarize(across(n_infected:outbreak_size,
#                      \(x) list(range(x, na.rm = TRUE)))) |>
#     unnest(n_infected:outbreak_size)
# #   n_infected p_emerge outbreak_size
# #        <dbl>    <dbl>         <dbl>
# # 1      -2.70   -0.438         -2.29
# # 2       5.18    0.96           3.97



# outbreak_heatmap("p_emerge", "weak", "Y0", "N0") /
#     outbreak_heatmap("p_emerge", "strong", "Y0", "N0") +
#     plot_layout(guides = "collect")
#
# outbreak_heatmap("outbreak_size", "weak", "Y0", "N0") /
#     outbreak_heatmap("outbreak_size", "strong", "Y0", "N0") +
#     plot_layout(guides = "collect")
#
#
# pseudo_eff_heatmap("p_emerge", "weak", "Y0", "N0") +
#     pseudo_eff_heatmap("p_emerge", "strong", "Y0", "N0") +
#     plot_layout(nrow = 1, guides = "collect")
#
# pseudo_eff_heatmap("outbreak_size", "weak", "Y0", "N0") +
#     pseudo_eff_heatmap("outbreak_size", "strong", "Y0", "N0") +
#     plot_layout(nrow = 1, guides = "collect")
#
# pseudo_eff_heatmap("n_infected", "weak", "Y0", "N0") +
#     pseudo_eff_heatmap("n_infected", "strong", "Y0", "N0") +
#     plot_layout(nrow = 1, guides = "collect")
#
# pseudo_eff_heatmap("outbreak_size", "weak", "Y0", "N0") +
#     pseudo_eff_heatmap("outbreak_size", "strong", "Y0", "N0") +
#     plot_layout(nrow = 1, guides = "collect")
#
# pseudo_eff_heatmap("n_infected", "weak", "Y0", "N0") +
#     pseudo_eff_heatmap("n_infected", "strong", "Y0", "N0") +
#     plot_layout(nrow = 1, guides = "collect")





#
# _ extremes ----
#


manip2_sims |>
    select(wasp_resp, Y0, N0, n_pseudo, n_infected:outbreak_size) |>
    pivot_longer(n_infected:outbreak_size, names_to = "outcome") |>
    # filter(!is.na(value)) |> ## << can't do this bc it mucks up first summarize
    mutate(outcome = factor(outcome, levels = c("p_emerge", "outbreak_size",
                                          "n_infected"))) |>
    group_by(wasp_resp, Y0, N0, outcome) |>
    summarize(value = value[n_pseudo == 3] - value[n_pseudo == 0],
              .groups = "drop") |>
    group_by(wasp_resp, outcome) |>
    summarize(min_val = min(value, na.rm = TRUE),
              max_val = max(value, na.rm = TRUE), .groups = "drop")




#
# _ supp. plots ----
#








# Heatmaps of Y0 and N0, with fill split into prob. emergence and outbreak size
yn_hm_p <- {pseudo_eff_heatmap("p_emerge", "weak", "Y0", "N0", .z_pal = "broc",
                               .add_title = scenario_title("weak")) +
        geom_segment(data = tibble(Y0 = 1, N0 = 75, Y0_s = Y0 + 0.5, N0_s = 150),
                     aes(Y0_s, N0_s, xend = Y0, yend = N0),
                     inherit.aes = FALSE,
                     color = "black", linewidth = 2, linejoin = "mitre",
                     arrow = arrow(type = "closed", length = unit(0.1, "inches"))) +
        geom_text(data = tibble(Y0 = 1, N0 = 200), aes(Y0, N0),
                  label = "original\nvalues", inherit.aes = FALSE,
                  hjust = 0, vjust = 1, color = "black",
                  size.unit = "pt", size = 14, lineheight = 0.7)} +
    pseudo_eff_heatmap("p_emerge", "strong", "Y0", "N0", .z_pal = "broc",
                       .add_title = scenario_title("strong")) +
    {pseudo_eff_heatmap("outbreak_size", "weak", "Y0", "N0", .z_pal = "vik") +
            geom_text(data = tibble(Y0 = 5.5, N0 = 60), aes(Y0, N0),
                      label = "no outbreaks", inherit.aes = FALSE,
                      hjust = 0.5, vjust = 0.5, color = "white",
                      size.unit = "pt", size = 14)} +
    pseudo_eff_heatmap("outbreak_size", "strong", "Y0", "N0", .z_pal = "vik") +
    plot_layout(nrow = 2, guides = "collect", axes = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"),
          plot.tag.location = "panel",
          plot.tag.position = c(0.05, 1.0))

if (.overwrite) {
    save_plot("_plots/extremes-manip-Y0-N0-heatmaps.pdf", yn_hm_p, 6.5, 5.5)
}



# ============================================================================*
# ============================================================================*
# Histograms of 2D diffs ----
# ============================================================================*
# ============================================================================*

manip2_sims |>
    select(wasp_resp, Y0, N0, n_pseudo, p_emerge, outbreak_size) |>
    group_by(wasp_resp, Y0, N0) |>
    summarize(p_emerge = p_emerge[n_pseudo > 0] - p_emerge[n_pseudo == 0],
              outbreak_size = outbreak_size[n_pseudo > 0] -
                  outbreak_size[n_pseudo == 0],
              .groups = "drop") |>
    pivot_longer(p_emerge:outbreak_size, names_to = "outcome") |>
    mutate(outcome = factor(outcome, levels = c("p_emerge", "outbreak_size"),
                            labels = c("Probabililty of emergence",
                                       "Mean outbreak size"))) |>
    filter(!is.na(value)) |>
    (\(x) {
        x |>
            group_by(wasp_resp, outcome) |>
            summarize(min = min(value),
                      mean = mean(value),
                      med = median(value),
                      max = max(value), .groups = "drop") |>
            print()
        return(x)
    })() |>
    ggplot(aes(value, after_stat(density))) +
    geom_vline(xintercept = 0, color = "gray70", linewidth = 1) +
    geom_hline(yintercept = 0, color = "gray70", linewidth = 1) +
    geom_freqpoly(aes(color = wasp_resp, linetype = wasp_resp, linewidth = wasp_resp),
                  bins = 25) +
    scale_linetype_manual("Parasitoid<br>responsiveness<br>to aphids",
                          values = c(strong = "22", weak = "solid")) +
    scale_linewidth_manual("Parasitoid<br>responsiveness<br>to aphids",
                           values = c(strong = 1.25, weak = 1)) +
    scale_color_manual("Parasitoid<br>responsiveness<br>to aphids",
                       values = c(strong = "firebrick4", weak = "dodgerblue")) +
    labs(x = "Effect of <i>Pseudomonas</i>") +
    facet_wrap(~ outcome, ncol = 1, scales = "free")






# ============================================================================*
# ============================================================================*
# zeta breakpoint ----
# ============================================================================*
# ============================================================================*

manip_sims |>
    filter(par_name == "zeta") |>
    group_by(par_val, n_pseudo) |>
    summarize(p_emerge = mean(n_infected > 1),
              outbreak_size = mean(n_infected[n_infected > 1]),
              .groups = "drop") |>
    pivot_longer(p_emerge:outbreak_size, names_to = "outcome") |>
    mutate(outcome = factor(outcome, levels = c("p_emerge", "outbreak_size"))) |>
    group_by(par_val, outcome) |>
    summarize(diff = value[n_pseudo > 0] - value[n_pseudo == 0],
              .groups = "drop") |>
    filter(par_val > 0.6, par_val < 0.8) |>
    ggplot(aes(par_val, diff)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_vline(xintercept = c(0.7, 0.688), color = "gray70", linetype = "22") +
    geom_line() +
    geom_point() +
    facet_wrap(~ outcome, ncol = 1)




# ============================================================================*
# ============================================================================*

# >> Main text plots ====

# ============================================================================*
# ============================================================================*

#' Variables:
#'
#' Y0
#' N0
#' K
#' pseudo_surv
#' zeta
#'
#' Of these, heatmaps:
#' Y0 + N0
#'

if (!dir.exists("_plots/extremes-manip-subs")) dir.create("_plots/extremes-manip-subs")





if (.overwrite) {
    # Line plot for zeta:
    save_plot("_plots/extremes-manip-subs/lines-zeta.pdf",
              one_manip_plotter(levels(manip2_sims$wasp_resp), "zeta",
                                outcomes = "n_infected", .x_pos = "top",
                                .seed = 65130342) &
                  illustrator_theme,
              width = 4, height = 2)
    # Heatmaps for N0 vs Y0:
    .p <- levels(wasp_resp_fct) |>
        map(\(.wr) pseudo_eff_heatmap("n_infected", .wr, "Y0", "N0"))  |>
        wrap_plots(design = "A#B", widths = c(1, 0.05, 1),
                   guides = "collect", axes = "collect")
    save_plot("_plots/extremes-manip-subs/heatmap-Y0-N0.pdf",
              .p & illustrator_theme, width = 4.15, height = 2)
    # save_plot("_plots/extremes-manip-subs/heatmap-Y0-N0-legend.pdf",
    #           cowplot::get_legend(.p & theme(legend.title = element_blank())),
    #           width = 1, height = 2)
    rm(.p)
}





ps_manip_p <- map2(levels(wasp_resp_fct), c(188673947, 899304975),
                   \(wr, seed) {
                       one_manip_plotter(wr, "pseudo_surv",
                                         outcomes = "n_infected",
                                         .seed = seed)
                   }) |>
    wrap_plots(design = "A#B", widths = c(1, 0.05, 1),
               guides = "collect", axes = "collect")
# ps_manip_p

if (.overwrite) {
    save_plot("_plots/extremes-manip-lines-pseudo_surv.pdf",
              ps_manip_p & illustrator_theme, width = 6, height = 2)
}
