
#'
#' Heatmaps and line graphs of outbreak sizes for simulations of two scenarios:
#' one where Pseudomonas decreases outbreaks ("low"),
#' and another where it increases outbreaks ("high")
#' ... where we manipulate one or two parameters at a time.
#'


source("_scripts/00-preamble.R")


# 1000 sims for each type of scenario:
extreme_sims <- read_rds(rds_files$extreme_large) |>
    mutate(outbreak_size = map_dbl(outbreak_size, mean))


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
                   zeta = seq(0, 1, 0.05),
                   spat_config = 1:length(spat_config_lvls)) |>
    map(\(x) if (inherits(x, "numeric")) round(x, 2) else x)



one_manip_sim <- function(type, x_name, x_val) {

    # type = "low"
    # x_name = "Y0"
    # x_val = manip_pars[[x_name]][[1]]
    # rm(type, x_name, x_val, args, sims_p, sims_np, out)

    args <- list(type = type, n_pseudo = 3L, large_sims = TRUE)
    # spat_config takes a character, but needs to be numeric to be compatible
    # with other parameter types:
    if (x_name == "spat_config") {
        args[[x_name]] <- spat_config_lvls[[x_val]]
    } else args[[x_name]] <- x_val

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


if (!file.exists(rds_files$extreme_manip)) {

    # Takes ~50 sec
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


one_manip_plotter <- function(type, par_name, .point = TRUE, .vline = TRUE) {
        # type = "low"; par_name = "spat_config"
        # rm(type, par_name, dd, dd_og, ..LINES, x_lvl_labs, dd_max_diff, p)
        dd <- manip_sims |>
            filter(type == .env$type, par_name == .env$par_name) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            group_by(n_pseudo, par_val) |>
            summarize(outbreak_size = mean(outbreak_size),
                      .groups = "drop")
        dd_og <- extreme_sims |>
            filter(type == .env$type) |>
            mutate(par_val = run_sim_combos(type = .env$type, n_pseudo = 0,
                                            return_args = TRUE) |>
                       getElement(par_name)) |>
            mutate(n_pseudo = factor(n_pseudo))
        ..LINES <- geom_line()
        if (par_name == "spat_config") {
            # x_lvl_labs <- spat_config_abbrevs
            x_lvl_labs <- tolower(as.roman(1:length(spat_config_lvls)))
            dd$par_val <- map_chr(dd$par_val, \(i) spat_config_lvls[[i]]) |>
                factor(levels = spat_config_lvls,
                       labels = x_lvl_labs)
            dd_og$par_val <- factor(dd_og$par_val, levels = spat_config_lvls,
                                    labels = x_lvl_labs)
            ..LINES <- theme(legend.position = "none")
        }
        dd_max_diff <- dd |>
            group_by(par_val) |>
            summarize(dos = outbreak_size[n_pseudo == 0L] -
                          outbreak_size[n_pseudo != 0L],
                      .groups = "drop") |>
            filter(abs(dos) == max(abs(dos))) |>
            getElement("par_val")
        p <- dd |>
            ggplot(aes(par_val, outbreak_size, color = n_pseudo)) +
            geom_hline(yintercept = c(1, 9), color = "gray70") +
            geom_point() +
            ..LINES
        if (.point) p <- p + geom_point(data = dd_og, size = 3, shape = 8)
        if (.vline) {
            p <- p +
                geom_vline(xintercept = dd_max_diff, color = "gray70",
                           linetype = "22")
        }
        p +
            labs(x = pretty_params(par_name) |> first_cap(),
                 y = "Outbreak size") +
            scale_y_continuous(breaks = c(1, 5, 9)) +
            coord_cartesian(ylim = c(1, 9)) +
            scale_color_manual("Number of<br>*Pseudomonas*<br>patches",
                               values = np_pal) +
            guides(color = guide_legend(override.aes = list(shape = 19)))
}




manip_plots <- c("low", "high") |>
    set_names() |>
    map(\(type) {

        # type = "high"
        # rm(type, .title, plot_list)

        .title <- sprintf("*Pseudomonas* %s outbreaks",
                          ifelse(type == "low", "inhibits", "promotes"))

        plot_list <- levels(manip_sims$par_name) |>
            map(\(x) one_manip_plotter(type, x))


        plot_list |>
            do.call(what = wrap_plots) +
            plot_layout(guides = "collect", axes = "collect") +
            plot_annotation(title = .title)
    })


# manip_plots$low
# manip_plots$high

# for (n in names(manip_plots)) {
#     save_plot(sprintf("_plots/extreme-manips-all-%s.pdf", n),
#               manip_plots[[n]],
#               width = 8, height = 5)
#     # save_plot(sprintf("_plots/extremes-manip-all-illustrator-%s.pdf", n),
#     #           manip_plots[[n]] & illustrator_theme &
#     #               theme(axis.title.x = element_markdown()),
#     #           width = 6.5, height = 4)
# }; rm(n)









# ============================================================================*
# ============================================================================*

# >> 2-par manips ====

# ============================================================================*
# ============================================================================*


# --------------------------------------*
# _ sims ----
# --------------------------------------*


manip2_pars <- combn(9, 2) |>
    t() |>
    (\(x) {colnames(x) <- paste0("par_name_", letters[1:2]); return(x)})() |>
    as_tibble() |>
    mutate(across(everything(), \(x) map_chr(x, \(i) names(manip_pars)[[i]]))) |>
    mutate(across(everything(), \(x) factor(x, levels = names(manip_pars))))




one_manip2_sim <- function(type, par_name_a, par_val_a, par_name_b, par_val_b) {

    # type = "low"
    # par_name_a = "Y0"; par_val_a = manip_pars[[par_name_a]][[1]]
    # par_name_b = "N0"; par_val_b = manip_pars[[par_name_b]][[1]]
    # rm(type, par_name_a, par_val_a, par_name_b, par_val_b, args, sims_p, sims_np, out)

    # In case these are factors:
    par_name_a <- paste(par_name_a)
    par_name_b <- paste(par_name_b)

    args <- list(type = type, n_pseudo = 3L, large_sims = TRUE)
    args[[par_name_a]] <- par_val_a
    args[[par_name_b]] <- par_val_b
    # spat_config takes a character, but needs to be numeric to be compatible
    # with other parameter types:
    if (par_name_a == "spat_config") args[[par_name_a]] <- spat_config_lvls[[par_val_a]]
    if (par_name_b == "spat_config") args[[par_name_b]] <- spat_config_lvls[[par_val_b]]

    sims_p <- do.call(run_sim_combos, args) |>
        select(-rep) |>
        summarize(across(everything(), mean)) |>
        mutate(n_pseudo = 3L)
    sims_np <- do.call(run_sim_combos, list_assign(args, n_pseudo = 0L)) |>
        select(-rep) |>
        summarize(across(everything(), mean)) |>
        mutate(n_pseudo = 0L)

    out <- bind_rows(sims_p, sims_np) |>
        mutate(type = .env$type,
               par_name_a = factor(.env$par_name_a, levels = names(manip_pars)),
               par_name_b = factor(.env$par_name_b, levels = names(manip_pars)),
               par_val_a = .env$par_val_a,
               par_val_b = .env$par_val_b) |>
        select(type, n_pseudo, starts_with("par_"), everything())

    return(out)

}




if (!file.exists(rds_files$extreme_manip2)) {

    # Takes ~1 hr
    set.seed(2025929231)
    manip2_sims <- manip2_pars |>
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
        mutate(type = factor(type, levels = c("low", "high")))

    write_rds(manip2_sims, rds_files$extreme_manip2, compress = "gz")

} else {

    manip2_sims <- read_rds(rds_files$extreme_manip2)

}




# --------------------------------------*
# _ supp. plots ----
# --------------------------------------*




# function to generate objects par_name_a, par_name_b, dd, pars_og, and .title,
# used for all heatmaps:
heatmaps_make_objs <- function(type, par_name_a, par_name_b, .add_title, env) {

    # In case these are factors:
    par_name_a <- paste(par_name_a)
    par_name_b <- paste(par_name_b)

    dd <- manip2_sims |>
        filter(type == .env$type,
               par_name_a == .env$par_name_a,
               par_name_b == .env$par_name_b)
    if (nrow(dd) == 0) {
        stop("\nCombination of type, par_name_a, and par_name_b not found!")
    }
    dd <- dd |>
        select(n_pseudo, starts_with("par_val_"), outbreak_size) |>
        mutate(n_pseudo = factor(n_pseudo))

    pars_og <- run_sim_combos(type = type, n_pseudo = 0,
                              return_args = TRUE) |>
        base::`[`(c(par_name_a, par_name_b)) |>
        set_names(c("par_val_a", "par_val_b")) |>
        as_tibble()
    if ("spat_config" %in% names(pars_og)) {
        pars_og[["spat_config"]] <- which(spat_config_lvls == pars_og[["spat_config"]])
    }

    .title <- waiver()
    if (.add_title) {
        .title <- sprintf("*Pseudomonas* %s outbreaks",
                          ifelse(type == "low", "inhibits", "promotes"))
    }

    assign("par_name_a", par_name_a, envir = env)
    assign("par_name_b", par_name_b, envir = env)
    assign("dd", dd, envir = env)
    assign("pars_og", pars_og, envir = env)
    assign(".title", .title, envir = env)

    invisible(NULL)

}


# function to add shared plot parts for heatmaps:
heatmaps_shared <- function(d, pars_og, par_name_a, par_name_b, .contour_args,
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
        ggplot(aes(par_val_a, par_val_b)) +
        geom_raster(aes(fill = outbreak_size))
    if (!isTRUE(is.na(.contour_args))) {
        .contour_args$mapping <- aes(z = outbreak_size)
        if (!"color" %in% names(.contour_args)) .contour_args$color <- "white"
        pp <- pp + do.call(geom_contour, .contour_args)
    }
    pp <- pp +
        # geom_vline(xintercept = pars_og[[par_name_a]], linetype = "22",
        #            color = "white", linewidth = 1) +
        # geom_hline(yintercept = pars_og[[par_name_b]], linetype = "22",
        #            color = "white", linewidth = 1) +
        geom_point(data = pars_og, size = 3, shape = 8, color = "white") +
        labs(x = x_lab, y = y_lab, tag = ..tag, title = ..title)

    spat_scale_fun <- NULL
    if (par_name_a == "spat_config") spat_scale_fun <- scale_x_continuous
    if (par_name_b == "spat_config") spat_scale_fun <- scale_y_continuous
    if (!is.null(spat_scale_fun)) {
        pp <- pp +
            spat_scale_fun(breaks = 1:length(spat_config_lvls),
                           labels = 1:length(spat_config_lvls) |>
                               as.roman() |> tolower())
    }

    return(pp)
}



# z is outbreak size, facets by n_pseudo:
outbreak_heatmap <- function(type, par_name_a, par_name_b,
                             .contour_args = list(), .tag = waiver(), .add_title = FALSE,
                             .n_pseudo = NA, .shorten_K = FALSE) {

    # generate objects par_name_a, par_name_b, dd, pars_og, and .title
    heatmaps_make_objs(type, par_name_a, par_name_b, .add_title, environment())

    z_breaks <- 0:4 * 2 + 1
    if (is.list(.contour_args) && !"breaks" %in% names(.contour_args)) {
        .contour_args$breaks <- z_breaks
    }

    if (is.na(.n_pseudo)) {
        p <- dd |>
            mutate(n_pseudo = factor(paste(n_pseudo), levels = levels(n_pseudo),
                                     labels = sprintf("n<sub>P</sub> = %s",
                                                      levels(n_pseudo)))) |>
            heatmaps_shared(pars_og, par_name_a, par_name_b,
                            .contour_args = .contour_args,
                            ..tag = .tag, ..title = .title, .shorten_K = .shorten_K) +
            facet_wrap(~ n_pseudo, nrow = 1)
    } else {
        p <- dd |>
            filter(n_pseudo == .n_pseudo) |>
            heatmaps_shared(pars_og, par_name_a, par_name_b,
                            .contour_args = .contour_args,
                            ..tag = .tag, ..title = .title, .shorten_K = .shorten_K)
    }
    p <- p +
        scale_fill_scico("Outbreak<br>size", limits = c(1,9), breaks = z_breaks,
                         palette = "tokyo", direction = -1)

    return(p)
}



# z is effect of n_pseudo on outbreak size:
pseudo_eff_heatmap <- function(type, par_name_a, par_name_b,
                               .contour_args = list(), .tag = waiver(),
                               .add_title = FALSE, .shorten_K = FALSE) {

    # generate objects par_name_a, par_name_b, dd, pars_og, and .title
    heatmaps_make_objs(type, par_name_a, par_name_b, .add_title, environment())

    z_breaks <- seq(-5, 5, 2.5)
    if (is.list(.contour_args) && !"breaks" %in% names(.contour_args)) {
        .contour_args$breaks <- z_breaks
    }

    dd |>
        group_by(across(starts_with("par_val")))  |>
        summarize(outbreak_size = outbreak_size[n_pseudo != "0"] -
                      outbreak_size[n_pseudo == "0"],
                  .groups = "drop") |>
        mutate(outbreak_size = round(outbreak_size, 3)) |>
        heatmaps_shared(pars_og, par_name_a, par_name_b, .contour_args = .contour_args,
                        ..tag = .tag, ..title = .title, .shorten_K = .shorten_K) +
        scale_fill_scico("Effect of<br>*Pseudomonas* on<br>outbreak size",
                         palette = "vik", midpoint = 0,
                         breaks = z_breaks, limits = c(-5, 5))
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






# if (!dir.exists("_plots/extremes-manip2")) dir.create("_plots/extremes-manip2", recursive = TRUE)
# for (x in distinct(manip2_sims, type, par_name_a, par_name_b) |>
#      filter(par_name_a != "spat_config", par_name_b != "spat_config") |>
#      (\(x) split(x, 1:nrow(x)))()) {
#     p <- one_manip2_plotter(x$type, x$par_name_a, x$par_name_b)
#     fn <- sprintf("_plots/extremes-manip2/%s-%s-%s.pdf", x$type, x$par_name_a, x$par_name_b)
#     save_plot(fn, p, 8, 3)
# }; rm(x, p, fn)

manip2_full_plots <- c("low", "high") |>
    set_names() |>
    map(\(type) {
        # type = "high"
        # rm(type, .title, plot_list)
        .title <- sprintf("*Pseudomonas* %s outbreaks",
                          ifelse(type == "low", "inhibits", "promotes"))
        plot_list <- manip2_pars |>
            filter(par_name_a != "spat_config", par_name_b != "spat_config") |>
            pmap(\(par_name_a, par_name_b) {
                pseudo_eff_heatmap(type, par_name_a, par_name_b,
                                   .shorten_K = TRUE) +
                    labs(x = ifelse(par_name_a == "K",
                                    serify("", "K &divide; 1000", ""),
                                    pretty_params(par_name_a, TRUE,
                                                  serif = TRUE)),
                         y = ifelse(par_name_b == "K",
                                    serify("", "K &divide; 1000", ""),
                                    pretty_params(par_name_b, TRUE,
                                                  serif = TRUE))) +
                    theme(legend.direction = "horizontal",
                          legend.title.position = "top")
            })
        plot_list |>
            c(list(guide_area())) |>
            do.call(what = wrap_plots) +
            plot_layout(guides = "collect", ncol = 6) +
            plot_annotation(title = .title) &
            theme(axis.title.x = element_markdown(size = 6),
                  axis.title.y = element_markdown(size = 6),
                  axis.text.x = element_markdown(size = 5),
                  axis.text.y = element_markdown(size = 5))
    })


# manip2_full_plots$low
# manip2_full_plots$high


# outbreak_heatmap("high", "Y0", "N0")
# pseudo_eff_heatmap("low", "Y0", "N0")



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

for (.t in c("low", "high")) {
    .f <- sprintf("_plots/extremes-manip-subs/heatmap-%s-Y0-N0.pdf", .t)
    .p <- pseudo_eff_heatmap(.t, "Y0", "N0", .contour_args = NA)
    save_plot(.f, .p + illustrator_theme, 2.5, 2.5)
    for (.v in c("K", "pseudo_surv", "zeta")) {
        .f <- sprintf("_plots/extremes-manip-subs/lines-%s-%s.pdf", .t, .v)
        .p <- one_manip_plotter(.t, .v, TRUE, FALSE)
        save_plot(.f, .p + illustrator_theme, 2.5, 1.25)
    }
}; rm(.t, .f, .p, .v)



