
#'
#' Heatmaps and line graphs of outbreak sizes for simulations of two scenarios:
#' one where Pseudomonas decreases outbreaks ("low"),
#' and another where it increases outbreaks ("high")
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
                   pseudo_surv = seq(0.85, 1, 0.01),
                   zeta = seq(0, 1, 0.05)) |>
    map(\(x) round(x, 2))



one_manip_sim <- function(type, x_name, x_val) {

    # type = "low"
    # x_name = "Y0"
    # x_val = manip_pars[[x_name]][[1]]
    # rm(type, x_name, x_val, args, sims_p, sims_np, out)

    args <- list(type = type, n_pseudo = 3L, large_sims = TRUE)
    args[[x_name]] <- x_val

    sims_p <- do.call(run_sim_combos, args) |>
        select(-rep) |>
        mutate(n_pseudo = 3L)
    sims_np <- do.call(run_sim_combos, list_assign(args, n_pseudo = 0L)) |>
        select(-rep) |>
        mutate(n_pseudo = 0L)

    out <- bind_rows(sims_p, sims_np) |>
        mutate(type = .env$type,
               par_name = factor(x_name, levels = names(manip_pars)),
               par_val = x_val) |>
        select(type, n_pseudo, par_name, par_val, everything())

    return(out)

}


if (.overwrite || !file.exists(rds_files$extreme_manip)) {

    # Takes ~6 min
    set.seed(1531497906)
    manip_sims <- c("low", "high") |>
        map(\(type) {
            manip_pars |>
                (\(x) {
                    out <- tibble(vals = do.call(c, manip_pars) |> unname())
                    out$par_name <- manip_pars |>
                        imap(\(x, n) rep(n, length(x))) |>
                        list_c()
                    return(out)
                })() |>
                pmap(\(vals, par_name) {
                    map(vals, \(v) one_manip_sim(type, par_name, v)) |>
                        list_rbind()
                }, .progress = .prog_args) |>
                list_rbind()
        }) |>
        list_rbind() |>
        mutate(type = factor(type, levels = c("low", "high")))

    write_rds(manip_sims, rds_files$extreme_manip, compress = "gz")

} else {

    manip_sims <- read_rds(rds_files$extreme_manip)

}






# --------------------------------------*
# _ supp. plots ----
# --------------------------------------*


one_manip_plotter <- function(type, par_name, .og_vals = TRUE, .md_vals = TRUE) {
    # type = "low"; par_name = "K"; .og_vals = TRUE; .md_vals = TRUE
    # rm(type, par_name, .og_vals, .md_vals, dd, pv_og, dd_og, yvar_pal, ..LINES)
    # rm(x_lvl_labs, x_scale, dd_md, p)
    dd <- manip_sims |>
        filter(type == .env$type, par_name == .env$par_name) |>
        mutate(n_pseudo = factor(n_pseudo)) |>
        group_by(n_pseudo, par_val) |>
        # (n_infected - 2) / 7  below is to have them on the same scale for
        # plotting. Axes labels will show differences properly
        summarize(p_emerge = mean(n_infected > 1),
                  outbreak_size = (mean(n_infected[n_infected > 1]) - 2) / 7,
                  .groups = "drop") |>
        pivot_longer(p_emerge:outbreak_size, names_to = "outcome")
    if (par_name == "K") dd$par_val <- carrying_capacity(K = dd$par_val, pred_surv = 1)
    pv_og <- run_sim_combos(type = type, n_pseudo = 0, return_args = TRUE) |>
        getElement(par_name)
    if (is.null(pv_og)) pv_og <- formals(lil_plantscape)[[par_name]]
    if (is.null(pv_og)) pv_og <- eval(formals(sim_plantscape)[[par_name]])
    if (is.null(pv_og)) pv_og <- eval(formals(make_insects_ptr)[[par_name]])
    if (is.null(pv_og)) pv_og <- eval(formals(make_disease_ptr)[[par_name]])
    if (par_name == "K") pv_og <- carrying_capacity(K = pv_og, pred_surv = 1)
    dd_og <- tibble(par_val = pv_og)

    yvar_pal <- brewer.pal(8, "Dark2")[c(1,8)] |>  # viridis(100)[c(10, 70)] |>
        set_names("outbreak_size", "p_emerge")

    ..LINES <- list(geom_line(aes(linetype = n_pseudo), linewidth = 0.5),
                    scale_linetype_manual("Number of<br>*Pseudomonas*<br>plants",
                                          values = np_linetypes))
    x_lvl_labs <- NULL
    x_scale <- scale_x_continuous()
    dd_md <- dd |>
        group_by(outcome, par_val) |>
        summarize(value0 = value[n_pseudo == 0L],
                  value = value[n_pseudo != 0L],
                  .groups = "drop_last") |>
        mutate(max_diff = abs(value0 - value) == max(abs(value0 - value),
                                                     na.rm = TRUE)) |>
        ungroup() |>
        filter(max_diff)


    if (par_name == "K") {
        x_lab <- paste("Aphid carrying capacity",
                       "(<span style=\"font-family: serif\"><i>K</i>",
                       "(<i>&lambda;</i> &minus; 1)</span>)")
    } else x_lab <- pretty_params(par_name) |> first_cap()


    p <- dd |>
        filter(!is.na(value)) |>
        ggplot(aes(par_val, value, color = outcome)) +
        geom_hline(yintercept = c(0, 1), color = "gray70")
    if (.og_vals) p <- p + geom_vline(data = dd_og[1,], aes(xintercept = par_val),
                                    color = "gray70", linewidth = 0.75)
    if (.md_vals) p <- p +
        geom_linerange(data = dd_md,
                      aes(x = par_val, ymin = value0, ymax = value),
                      linetype = "21", linewidth = 1)
    p +
        geom_point(aes(shape = n_pseudo)) +
        ..LINES +
        labs(x = x_lab, y = yvar_desc[["p_emerge"]] |> first_cap()) +
        x_scale +
        scale_y_continuous(breaks = 0:2/2,
                           sec.axis = sec_axis(\(x) 7 * x + 2,
                                               paste("Mean",
                                                     yvar_desc[["outbreak_size"]]),
                                               breaks = c(3,6,9))) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_color_manual(values = yvar_pal, guide = "none") +
        scale_shape_manual("Number of<br>*Pseudomonas*<br>plants",
                           values = np_shapes) +
        guides(shape = guide_legend(override.aes = list(color = "black"))) +
        theme(axis.title.y.left = element_text(color = yvar_pal[["p_emerge"]]),
              axis.text.y.left = element_text(color = yvar_pal[["p_emerge"]]),
              axis.ticks.y.left = element_line(color = yvar_pal[["p_emerge"]]),
              axis.title.y.right = element_text(color = yvar_pal[["outbreak_size"]]),
              axis.text.y.right = element_text(color = yvar_pal[["outbreak_size"]]),
              axis.ticks.y.right = element_line(color = yvar_pal[["outbreak_size"]]))
}




# one_manip_plotter("high", "zeta", TRUE, FALSE)


# manip_plots <- c("low", "high") |>
#     set_names() |>
#     map(\(type) {
#
#         # type = "low"
#         # rm(type, .title, plot_list)
#
#         .title <- scenario_title(type)
#
#         plot_list <- levels(manip_sims$par_name) |>
#             map(\(x) one_manip_plotter(type, x))
#
#         plot_list |>
#             do.call(what = wrap_plots) +
#             plot_layout(nrow = 4, guides = "collect", axes = "collect") +
#             plot_annotation(title = .title)
#     })
#
#
# # manip_plots$low
# # manip_plots$high
#
# if (.overwrite) {
#     for (n in names(manip_plots)) {
#         save_plot(sprintf("_plots/extreme-manips-all-%s.pdf", n),
#                   manip_plots[[n]],
#                   width = 6.5, height = 6)
#         # save_plot(sprintf("_plots/extremes-manip-all-illustrator-%s.pdf", n),
#         #           manip_plots[[n]] & illustrator_theme &
#         #               theme(axis.title.x = element_markdown()),
#         #           width = 6.5, height = 4)
#     }; rm(n)
# }









# ============================================================================*
# ============================================================================*

# >> 2-par manips ====

# ============================================================================*
# ============================================================================*


# --------------------------------------*
# _ sims ----
# --------------------------------------*




one_manip2_sim <- function(type, par_name_a, par_val_a, par_name_b, par_val_b) {

    # type = "low"
    # par_name_a = "Y0"; par_val_a = manip_pars[[par_name_a]][[1]]
    # par_name_b = "N0"; par_val_b = manip_pars[[par_name_b]][[1]]
    # rm(type, par_name_a, par_val_a, par_name_b, par_val_b, args, sims_p, sims_np, out)
    #

    # In case these are factors:
    par_name_a <- paste(par_name_a)
    par_name_b <- paste(par_name_b)

    args <- list(type = type, large_sims = TRUE)
    args[[par_name_a]] <- par_val_a
    args[[par_name_b]] <- par_val_b

    out <- c(3L, 0L) |>
        map(\(np) {
            .np <- as.integer(np)
            do.call(run_sim_combos, list_assign(args, n_pseudo = .np)) |>
                select(-rep) |>
                mutate(p_emerge = as.integer(n_infected > 1),
                       outbreak_size = ifelse(n_infected > 1, n_infected, NA)) |>
                summarize(across(everything(), \(x) mean(x, na.rm = TRUE))) |>
                mutate(n_pseudo = .np)
        }) |>
        list_rbind() |>
        mutate(type = .env$type,
               par_name_a = factor(.env$par_name_a, levels = names(manip_pars)),
               par_name_b = factor(.env$par_name_b, levels = names(manip_pars)),
               par_val_a = .env$par_val_a,
               par_val_b = .env$par_val_b) |>
        select(type, n_pseudo, starts_with("par_"), everything())

    return(out)

}




if (!file.exists(rds_files$extreme_manip2) || .overwrite) {

    # Takes ~12 min
    set.seed(2025929231)
    manip2_sims <- tibble(par_name_a = "Y0",
                          par_name_b = "N0") |>
        mutate(type = list(c("low", "high"))) |>
        unnest(type) |>
        mutate(vals = map2(par_name_a, par_name_b,
                           \(par_name_a, par_name_b) {
                               crossing(par_val_a = manip_pars[[par_name_a]],
                                        par_val_b = manip_pars[[par_name_b]])
                           })) |>
        unnest(vals) |>
        arrange(type, par_name_a, par_name_b, par_val_a, par_val_b) |>
        select(type, par_name_a, par_val_a, par_name_b, par_val_b) |>
        pmap(one_manip2_sim, .progress = .prog_args) |>
        list_rbind() |>
        mutate(type = factor(type, levels = c("low", "high"))) |>
        rename(Y0 = par_val_a, N0 = par_val_b) |>
        select(-starts_with("par_name_"))

    write_rds(manip2_sims, rds_files$extreme_manip2, compress = "gz")

} else {

    manip2_sims <- read_rds(rds_files$extreme_manip2)

}



#
# _ N0 / Y0 thresholds ----
#

# ... for when outbreaks never occur
manip2_sims |>
    group_by(type, N0) |>
    filter(Y0 == min(Y0[p_emerge == 0])) |>
    group_by(type) |>
    summarize(rel = median(Y0 / N0))
manip2_sims |>
    group_by(type, N0) |>
    filter(Y0 == min(Y0[p_emerge == 0])) |>
    ungroup() |>
    ggplot(aes(Y0, N0)) +
    geom_point() +
    facet_wrap(~ type, nrow = 1)

# ... for when outbreaks always occur
manip2_sims |>
    filter(p_emerge == 1) |>
    group_by(type, N0) |>
    filter(Y0 == max(Y0)) |>
    group_by(type) |>
    summarize(rel = median(Y0 / N0))
manip2_sims |>
    filter(p_emerge == 1) |>
    group_by(type, N0) |>
    filter(Y0 == max(Y0)) |>
    ungroup() |>
    ggplot(aes(Y0, N0)) +
    geom_point() +
    facet_wrap(~ type, nrow = 1)



#
# _ above / below thresholds ----
#
set.seed(778424076)
thresh_sims <- crossing(type = sort(unique(manip_sims$type)),
                        n_pseudo = c(0L, 3L),
                        large_sims = FALSE,
                        ny = factor(0:1, labels = c("below", "above"))) |>
    pmap(\(type, n_pseudo, large_sims, ny) {
        Y0  <- ifelse(ny == "below", 3, 1)
        N0 <- ifelse(ny == "below", 50, 150)
        out <- run_sim_combos(type = type,
                              n_pseudo = n_pseudo,
                              large_sims = large_sims,
                              Y0 = Y0,
                              N0 = N0) |>
            select(-rep) |>
            mutate(type = .env$type,
                   n_pseudo = .env$n_pseudo,
                   thresh = ny)
    }) |>
    list_rbind()


thresh_empty_data <- thresh_sims |>
    filter(is.na(x)) |>
    mutate(aphids = aphids + parasitized) |>
    select(aphids, alates, wasps, virus) |>
    pivot_longer(aphids:virus, names_to = "species", values_to = "density") |>
    mutate(species = factor(species, levels =
                                c("aphids", "alates", "wasps", "virus"))) |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 1)

thresh_p_list <- thresh_sims |>
    distinct(thresh, type) |>
    arrange(thresh, type) |>
    pmap(\(type, thresh) {
        # type = "low"; thresh = "above"; .title = NULL
        # rm(type, thresh, .title, .subtitle, .strip)
        type <- paste(type)
        thresh <- paste(thresh)
        if (thresh == "below") {
            .title <- scenario_title(type, TRUE, TRUE)
        } else .title <- waiver()
        .subtitle <- serify("", "*Y*<sub>0</sub> / *N*<sub>0</sub>",
                           list(below = " too high", above = " too low")[[thresh]])
        if (type == "low") {
            .strip <- element_blank()
        } else .strip <- element_markdown(size = 8)
        thresh_sims |>
            filter(type == .env$type, thresh == .env$thresh, is.na(x)) |>
            mutate(aphids = aphids + parasitized) |>
            select(n_pseudo, time, aphids, alates, wasps, virus) |>
            pivot_longer(aphids:virus, names_to = "species",
                         values_to = "density") |>
            mutate(n_pseudo = factor(n_pseudo),
                   species = factor(species, levels =
                                        c("aphids", "alates", "wasps", "virus"))) |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_blank(data = thresh_empty_data) +
            geom_line(aes(color = species, linetype = n_pseudo,
                          linewidth = n_pseudo)) +
            scale_color_manual(values = color_pal, guide = "none") +
            scale_linetype_manual("*Pseudo.*<br>plants", values = np_linetypes) +
            scale_linewidth_manual("*Pseudo.*<br>plants", values = np_linewidths) +
            labs(x = "Time (days)", y = "Density", title = .title,
                 subtitle = .subtitle) +
            scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
            coord_cartesian(clip = "off") +
            facet_wrap(~ species, ncol = 1, scales = "free_y",
                       strip.position = "right") +
            guides(linetype = guide_legend(
                override.aes = list(color = "black"))) +
            theme(strip.text.y = .strip,
                  axis.text.x = element_markdown(size = 7),
                  axis.text.y = element_markdown(size = 7),
                  plot.subtitle = element_markdown())
    })




thresh_p <- wrap_plots(thresh_p_list, ncol = 2, guides = "collect", axes = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))


if (.overwrite) {
    save_plot("_plots/extreme-manips-thresholds.pdf", thresh_p, width = 6.5, height = 6)
}








# --------------------------------------*
# _ supp. plots ----
# --------------------------------------*




# function to generate objects par_name_a, par_name_b, dd, pars_og, and .title,
# used for all heatmaps:
heatmaps_make_objs <- function(yvar,
                               type,
                               par_name_a,
                               par_name_b,
                               .add_title,
                               env) {

    # In case these are factors:
    par_name_a <- paste(par_name_a)
    par_name_b <- paste(par_name_b)

    stopifnot(!"par_name_a" %in% colnames(manip2_sims))
    dd <- manip2_sims |>
        filter(type == .env$type)
    if (nrow(dd) == 0) {
        stop("\nCombination of type, par_name_a, and par_name_b not found!")
    }
    dd <- dd |>
        select(n_pseudo, all_of(c(par_name_a, par_name_b, yvar))) |>
        mutate(n_pseudo = factor(n_pseudo))

    pars_og <- run_sim_combos(type = type, n_pseudo = 0,
                              return_args = TRUE) |>
        base::`[`(c(par_name_a, par_name_b)) |>
        as_tibble()

    .title <- waiver()
    if (.add_title) {
        .title <- scenario_title(type)
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



# z is outbreak size or Pr(emerge), facets by n_pseudo:
outbreak_heatmap <- function(yvar, type, par_name_a, par_name_b,
                             .contour_args = list(),
                             .tag = waiver(),
                             .add_title = FALSE,
                             .n_pseudo = NA,
                             .shorten_K = FALSE) {

    # yvar = "p_emerge"; type = "low"; par_name_a = "Y0"; par_name_b = "N0"
    # .contour_args = list(); .tag = waiver(); .add_title = FALSE
    # .n_pseudo = NA; .shorten_K = FALSE
    # rm(yvar, type, par_name_a, par_name_b, .contour_args, .tag, .add_title)
    # rm(.n_pseudo, .shorten_K, dd, pars_og, .title, z_lab, z_breaks, z_lims, p)

    # generate objects par_name_a, par_name_b, dd, pars_og, and .title
    heatmaps_make_objs(yvar, type, par_name_a, par_name_b, .add_title,
                       environment())

    if (yvar == "outbreak_size") {
        z_lab <- "Outbreak<br>size"
        z_breaks <- 1:4 * 2 + 1
        z_lims <- c(2, 9)
    } else if (yvar == "p_emerge") {
        z_lab <- "Emergence<br>prob."
        z_breaks <- 0:4 * 0.25
        z_lims <- c(0, 1)
    } else stop("only yvar == \"p_emerge\" or \"outbreak_size\" is programmed")
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


# z is effect of n_pseudo on outbreak size or Pr(emerge):
pseudo_eff_heatmap <- function(yvar, type, par_name_a, par_name_b,
                               .contour_args = list(), .tag = waiver(),
                               .add_title = FALSE, .shorten_K = FALSE) {

    # yvar = "p_emerge"; type = "low"; par_name_a = "Y0"; par_name_b = "N0"
    # .contour_args = list(); .tag = waiver(); .add_title = FALSE
    # .shorten_K = FALSE
    # rm(yvar, type, par_name_a, par_name_b, .contour_args, .tag, .add_title)
    # rm(.shorten_K, dd, pars_og, .title, z_lab, z_breaks, z_lims, y_summ)

    # generate objects par_name_a, par_name_b, dd, pars_og, and .title
    heatmaps_make_objs(yvar, type, par_name_a, par_name_b, .add_title,
                       environment())

    if (yvar == "outbreak_size") {
        z_lab <- "Effect of<br>*Pseudomonas* on<br>outbreak size"
        z_breaks <- -2:2 * 1.5
        z_lim <- c(-1, 1) * 3
        z_pal = "vik"
        z_dir = 1
    } else if (yvar == "p_emerge") {
        z_lab <- "Effect of<br>*Pseudomonas* on<br>emergence prob."
        z_breaks <- -2:2 * 0.4
        z_lim <- c(-1, 1) * 0.8
        z_pal = "bam"
        z_dir = -1
    } else stop("only yvar == \"p_emerge\" or \"outbreak_size\" is programmed")
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
        scale_fill_scico(z_lab, palette = z_pal, midpoint = 0,
                         direction = z_dir, breaks = z_breaks, limits = z_lim)
}




one_manip2_plotter <- function(type, par_name_a, par_name_b,
                               .contour_args = list(), .tag = waiver(),
                               .add_title = FALSE, .shorten_K = FALSE) {

    # generate objects par_name_a, par_name_b, dd, pars_og, and .title
    heatmaps_make_objs(type, par_name_a, par_name_b, .add_title, environment())

    # Separate by n_pseudo:
    p1 <- outbreak_heatmap(type, par_name_a, par_name_b, .contour_args, .tag) +
        theme(plot.title = element_blank())
    # Effect of n_pseudo:
    p2 <- pseudo_eff_heatmap(type, par_name_a, par_name_b, .contour_args) +
        theme(plot.title = element_blank())

    p1 + p2 +
        plot_layout(nrow = 1, widths = c(2, 1), axes = "collect") +
        plot_annotation(title = .title)
}


manip2_sims |>
    group_by(type, Y0, N0) |>
    summarize(outbreak_size = outbreak_size[n_pseudo != "0"] -
                  outbreak_size[n_pseudo == "0"],
              p_emerge = p_emerge[n_pseudo != "0"] -
                  p_emerge[n_pseudo == "0"],
              .groups = "drop") |>
    summarize(os_min = min(outbreak_size, na.rm = TRUE),
              os_max = max(outbreak_size, na.rm = TRUE),
              pe_min = min(p_emerge),
              pe_max = max(p_emerge))
#   os_min os_max pe_min pe_max
#    <dbl>  <dbl>  <dbl>  <dbl>
# 1  -2.33   1.96  -0.49  0.652



# outbreak_heatmap("p_emerge", "low", "Y0", "N0") /
#     outbreak_heatmap("p_emerge", "high", "Y0", "N0") +
#     plot_layout(guides = "collect")
#
# outbreak_heatmap("outbreak_size", "low", "Y0", "N0") /
#     outbreak_heatmap("outbreak_size", "high", "Y0", "N0") +
#     plot_layout(guides = "collect")
#
#
# (pseudo_eff_heatmap("p_emerge", "low", "Y0", "N0") |
#     pseudo_eff_heatmap("p_emerge", "high", "Y0", "N0")) +
#     plot_layout(guides = "collect")
#
# (pseudo_eff_heatmap("outbreak_size", "low", "Y0", "N0") |
#     pseudo_eff_heatmap("outbreak_size", "high", "Y0", "N0")) +
#     plot_layout(guides = "collect")



# ============================================================================*
# ============================================================================*
# Histograms of 2D diffs
# ============================================================================*
# ============================================================================*

yvar <- "p_emerge"

c("low", "high") |>
    set_names() |>
    map(\(type) {
        y_summ <- \(y, np) round(y[np != "0"] - y[np == "0"], 3)
        manip2_sims |>
            filter(type == .env$type) |>
            group_by(type, Y0, N0) |>
            summarize(diff = y_summ(.data[[yvar]], np = n_pseudo), .groups = "drop")
    }) |>
    list_rbind() |>
    filter(!is.na(diff)) |>
    (\(x) {
        x |>
            group_by(type) |>
            summarize(min = min(diff),
                      mean = mean(diff),
                      med = median(diff),
                      max = max(diff)) |>
            print()
        return(x)
    })() |>
    ggplot(aes(diff, after_stat(density))) +
    geom_vline(xintercept = 0, linetype = "22", color = "gray70", linewidth = 1) +
    geom_freqpoly(aes(color = type), bins = 25, linewidth = 1) +
    scale_color_manual(values = c(low = "firebrick4", high = "dodgerblue")) +
    labs(x = paste("Effect of <i>Pseudomonas</i> on", yvar_desc[[yvar]]))




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
    for (.t in c("low", "high")) {
        for (.v in c("outbreak_size", "p_emerge")) {
            .f <- sprintf("_plots/extremes-manip-subs/heatmap-%s-%s-Y0-N0.pdf",
                          list(outbreak_size = "obs", p_emerge = "pem")[[.v]], .t)
            .p <- pseudo_eff_heatmap(.v, .t, "Y0", "N0", .contour_args = NA)
            save_plot(.f, .p + illustrator_theme, 2, 2)
        }
        for (.v in c("K", "pseudo_surv", "zeta")) {
            .f <- sprintf("_plots/extremes-manip-subs/lines-%s-%s.pdf", .t, .v)
            .p <- one_manip_plotter(.t, .v, TRUE, FALSE)
            save_plot(.f, .p + illustrator_theme, 2.5, 1.25)
        }
    }; rm(.t, .f, .p, .v)
}



# for (.v in c("outbreak_size", "p_emerge")) {
#     cowplot::get_legend(pseudo_eff_heatmap(.v, "low", "Y0", "N0",
#                                            .contour_args = NA)) |>
#         save_plot(filename = sprintf("~/Desktop/legend-%s.pdf", .v), width = 1, height = 8)
# }
