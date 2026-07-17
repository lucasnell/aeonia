


# function to generate objects par_name_a, par_name_b, dd, pars_og, and .title,
# used for all heatmaps:
heatmaps_make_objs <- function(yvar,
                               wasp_resp,
                               par_name_a,
                               par_name_b,
                               .add_title,
                               env) {

    # In case these are factors:
    par_name_a <- paste(par_name_a)
    par_name_b <- paste(par_name_b)

    stopifnot(!"par_name_a" %in% colnames(manip2_sims))
    dd <- manip2_sims |>
        filter(wasp_resp == .env$wasp_resp)
    if (nrow(dd) == 0) {
        stop("\nCombination of wasp_resp, par_name_a, and par_name_b not found!")
    }
    dd <- dd |>
        select(n_pseudo, all_of(c(par_name_a, par_name_b, yvar, "p_emerge"))) |>
        mutate(n_pseudo = factor(n_pseudo))

    pars_og <- run_sim_combos(wasp_resp = wasp_resp, n_pseudo = 0,
                              return_args = TRUE) |>
        base::`[`(c(par_name_a, par_name_b)) |>
        as_tibble()

    .title <- waiver()
    if (is.character(.add_title)) {
        .title <- .add_title
    } else if (is.logical(.add_title) && .add_title) {
        .title <- scenario_title(wasp_resp)
    }

    assign("par_name_a", par_name_a, envir = env)
    assign("par_name_b", par_name_b, envir = env)
    assign("dd", dd, envir = env)
    assign("pars_og", pars_og, envir = env)
    assign(".title", .title, envir = env)

    invisible(NULL)

}


# function to add shared plot parts for heatmaps:
heatmaps_shared <- function(d, yvar, pars_og, par_name_a, par_name_b, .contour_args,
                            ..tag = waiver(), ..title = waiver(),
                            .shorten_K = FALSE) {

    x_lab <- pretty_params(par_name_a) |> first_cap()
    y_lab <- pretty_params(par_name_b) |> first_cap()
    if (.shorten_K && (par_name_a == "K" || par_name_b == "K")) {
        if (par_name_a == "K") {
            d$par_val_a <- d$par_val_a / 1000
            x_lab <- "K &divide; 1000"
            pars_og[["par_val_a"]] <- pars_og[["par_val_a"]] / 1000
        }
        if (par_name_b == "K") {
            d$par_val_b <- d$par_val_b / 1000
            y_lab <- "K &divide; 1000"
            pars_og[["par_val_b"]] <- pars_og[["par_val_b"]] / 1000
        }
    }


    pp <- d |>
        ggplot(aes(.data[[par_name_a]], .data[[par_name_b]])) +
        geom_raster(aes(fill = .data[[yvar]]), na.rm = TRUE)
    if (!isTRUE(is.na(.contour_args))) {
        .contour_args$mapping <- aes(z = .data[[yvar]])
        if (!"color" %in% names(.contour_args)) .contour_args$color <- "white"
        if (!"na.rm" %in% names(.contour_args)) .contour_args$na.rm <- TRUE
        # check for whether values in `d[[yvar]]` overlap >= 1 break
        do_contours <- TRUE
        if ("breaks" %in% names(.contour_args)) {
            z_breaks <- .contour_args$breaks
            min_d <- min(d[[yvar]], na.rm = TRUE)
            max_d <- max(d[[yvar]], na.rm = TRUE)
            do_contours <- sum(z_breaks <= max_d & z_breaks >= min_d) >= 1
        }
        # If this is the case (or if no breaks provided), add contours:
        if (do_contours) pp <- pp + do.call(geom_contour, .contour_args)
    }
    pp <- pp +
        # geom_vline(xintercept = pars_og[[par_name_a]], linetype = "22",
        #            color = "white", linewidth = 1) +
        # geom_hline(yintercept = pars_og[[par_name_b]], linetype = "22",
        #            color = "white", linewidth = 1) +
        geom_point(data = pars_og, size = 2, shape = 23, color = "black",
                   fill = "white", stroke = 1) +
        labs(x = x_lab, y = y_lab, tag = ..tag, title = ..title)

    return(pp)
}



#' `yvar` is Pr(emerge), outbreak size, or # infected
#' facets are by n_pseudo:
outbreak_heatmap <- function(yvar, wasp_resp, par_name_a, par_name_b,
                             .contour_args = list(),
                             .tag = waiver(),
                             .add_title = FALSE,
                             .n_pseudo = NA,
                             .shorten_K = FALSE) {

    # yvar = "p_emerge"; wasp_resp = "strong"; par_name_a = "Y0"; par_name_b = "N0"
    # .contour_args = list(); .tag = waiver(); .add_title = FALSE
    # .n_pseudo = NA; .shorten_K = FALSE
    # rm(yvar, wasp_resp, par_name_a, par_name_b, .contour_args, .tag, .add_title)
    # rm(.n_pseudo, .shorten_K, dd, pars_og, .title, z_lab, z_breaks, z_lims, p)

    # generate objects par_name_a, par_name_b, dd, pars_og, and .title
    heatmaps_make_objs(yvar, wasp_resp, par_name_a, par_name_b, .add_title,
                       environment())

    if (yvar == "outbreak_size") {
        z_lab <- "Outbreak<br>size"
        z_breaks <- 1:4 * 2 + 1
        z_lims <- c(2, 9)
    } else if (yvar == "p_emerge") {
        z_lab <- "Emergence<br>prob."
        z_breaks <- 0:4 * 0.25
        z_lims <- c(0, 1)
    } else if (yvar == "n_infected") {
        z_lab <- "Infected<br>plants"
        z_breaks <- c(1, 3, 5, 7, 9)
        z_lims <- c(1, 9)
    } else {
        stop("\nonly yvar == \"p_emerge\", \"outbreak_size\", and ",
             "\"n_infected\" are programmed")
    }
    if (is.list(.contour_args) && !"breaks" %in% names(.contour_args)) {
        .contour_args$breaks <- z_breaks
    }

    if (is.na(.n_pseudo)) {
        p <- dd |>
            mutate(n_pseudo = factor(paste(n_pseudo), levels = levels(n_pseudo),
                                     labels = sprintf("n<sub>P</sub> = %s",
                                                      levels(n_pseudo)))) |>
            heatmaps_shared(yvar, pars_og, par_name_a, par_name_b,
                            .contour_args = .contour_args,
                            ..tag = .tag, ..title = .title, .shorten_K = .shorten_K) +
            facet_wrap(~ n_pseudo, nrow = 1)
    } else {
        p <- dd |>
            filter(n_pseudo == .n_pseudo) |>
            heatmaps_shared(yvar, pars_og, par_name_a, par_name_b,
                            .contour_args = .contour_args,
                            ..tag = .tag, ..title = .title, .shorten_K = .shorten_K)
    }
    p <- p +
        scale_fill_scico(z_lab, limits = z_lims, breaks = z_breaks,
                         palette = "tokyo", direction = -1)

    return(p)
}




#' `yvar` is Pr(emerge), outbreak size, or # infected
pseudo_eff_heatmap <- function(yvar, wasp_resp, par_name_a, par_name_b,
                               .contour_args = NA, .tag = waiver(),
                               .add_title = FALSE, .shorten_K = FALSE,
                               .z_pal = "bam") {

    # yvar = "n_infected"; wasp_resp = "strong"; par_name_a = "Y0"; par_name_b = "N0"
    # .contour_args = NA; .tag = waiver(); .add_title = FALSE
    # .shorten_K = FALSE
    # rm(yvar, wasp_resp, par_name_a, par_name_b, .contour_args, .tag, .add_title)
    # rm(.shorten_K, dd, pars_og, .title, z_lab, z_breaks, z_lims, y_summ)

    # generate objects par_name_a, par_name_b, dd, pars_og, and .title
    heatmaps_make_objs(yvar, wasp_resp, par_name_a, par_name_b, .add_title,
                       environment())

    z_dir = -1
    if (yvar == "outbreak_size") {
        z_lab <- "Effect of<br>*Pseudo.* on<br>outbreak size"
        z_breaks <- -2:2 * 2
        z_lim <- c(-1, 1) * 4
    } else if (yvar == "p_emerge") {
        z_lab <- "Effect of<br>*Pseudo.* on<br>emergence prob."
        z_breaks <- -2:2 * 0.5
        z_lim <- c(-1, 1) * 1
    } else if (yvar == "n_infected") {
        z_lab <- "Effect of<br>*Pseudo.* on<br>peak infected<br>plants"
        z_breaks <- -1:2 * 2
        z_lim <- c(-2.75, 5.2)
    } else {
        stop("\nonly yvar == \"p_emerge\", \"outbreak_size\", and ",
             "\"n_infected\" are programmed")
    }
    if (is.list(.contour_args) && !"breaks" %in% names(.contour_args)) {
        .contour_args$breaks <- z_breaks
    }

    y_summ <- \(y, np) round(y[np != "0"] - y[np == "0"], 3)
    dd |>
        group_by(across(all_of(c(par_name_a, par_name_b)))) |>
        summarize(across(all_of(yvar), \(y) y_summ(y, np = n_pseudo)),
                  .groups = "drop") |>
        heatmaps_shared(yvar, pars_og, par_name_a, par_name_b,
                        .contour_args = .contour_args,
                        ..tag = .tag, ..title = .title, .shorten_K = .shorten_K) +
        scale_fill_scico(z_lab, palette = .z_pal, midpoint = 0,
                         direction = z_dir, breaks = z_breaks, limits = z_lim,
                         guide = guide_colourbar(theme = theme(
                             # legend.key.size = unit(12, "in"),
                             legend.key.height = unit(1.2, "in"),
                             legend.key.width  = unit(0.24, "in")
                         ))) +
        # This fixes a very annoying bug where using
        # `legend.text = ggtext::element_markdown(...)` makes the legend
        # really large and doesn't allow me to change its size
        replace_theme(theme_get(),
                      legend.text = element_text(color = "black", size = 9))
}


