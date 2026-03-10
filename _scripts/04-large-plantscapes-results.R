
#'
#' Plots for larger landscape simulations
#'

source("_scripts/00-preamble.R")




#' Because below returns an empty tibble, it shows that when we set
#' `wt_vp = 1e-6` it does indeed result in no overlap between
#' *Pseudomonas* and virus.
#' Unsurprisingly, when we set `wt_vp = 100` and manually have them overlap,
#' there is always overlap.
#'
# read_rds("_scripts/interm-data/large-plantscapes.rds") |>
#     select(n_pseudo:wt_pp, landscape) |>
#     filter(!is.na(wt_pp)) |>
#     # Predicted number of instances of both Pseudomonas and virus
#     # across 100 sims:
#     mutate(n_both_pred = ifelse(wt_vp < 1, 0L, 100L)) |>
#     # Observed:
#     mutate(n_both_obs = map_int(landscape, \(x) sum(x == 3L))) |>
#     select(n_pseudo:wt_pp, starts_with("n_both")) |>
#     filter(n_both_pred != n_both_obs)





# from 03-large-plantscapes.sh:
sim_df <- read_rds("_scripts/interm-data/large-plantscapes.rds") |>
    # Remove bc this vector is >1 GB
    select(-landscape) |>
    mutate(outbreak_size = map_dbl(sim, \(x) mean(x$outbreak_size[x$outbreak_size > 1])),
           log_outbreak_size = map_dbl(sim, \(x) mean(log10(x$outbreak_size[x$outbreak_size > 1]))),
           p_emerge = map_dbl(sim, \(x) mean(x$outbreak_size > 1))) |>
    mutate(wt_vp = ifelse(wt_vp < 1, "off virus", "on virus") |>
               factor(levels = c("off virus", "on virus")),
           wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
               factor(levels = c("uniform", "clustered")))


# Used for a couple of plot fxns below:
pretty_greek_fct <- function(x, letter, suffix = "") {
    lvls <- sort(unique(x))
    factor(x, levels = lvls,
           labels = serify("", paste0("&", letter, ";", suffix),
                           sprintf(" = %.1f", lvls)))
}
make_land_type <- function(wt_vp, wt_pp, sep = "<br>") {
    lvls <- paste(rep(c("off virus", "on virus"), each = 2),
                  rep(c("uniform", "clustered"), 2), sep = "<br>")
    factor(paste0(wt_vp, sep, wt_pp), levels = lvls)
}


#' Basic check to make sure that the two scenarios still do what we expect
#' them to: Pseudomonas decreases outbreak size in "low" and increases in "high"

sim_df |>
    filter(sd_N == 0 & virus_attract == 1 & pseudo_repel == 1 &
               ((wt_vp == "off virus" & wt_pp == "uniform" & n_pseudo == 7e3) |
                    n_pseudo == 0)) |>
    select(type, n_pseudo, outbreak_size) |>
    arrange(type, n_pseudo)



big_land_plotter <- function(yvar,
                             type,
                             col_fct = "sd_N",
                             x_facet_fct = "land_type",
                             y_facet_fct = c("virus_attract", "pseudo_repel"),
                             shp_lty_fct = NULL,
                             fixed = NULL,
                             color_vals = NULL,
                             facet_scales = "fixed",
                             facet_ncol = NULL,
                             y_breaks = waiver(),
                             y_max = NULL) {

    # yvar = "outbreak_size"; type = "low"; col_fct = "pseudo_repel"
    # x_facet_fct = "virus_attract"; y_facet_fct = NULL; shp_lty_fct = "wt_vp"
    # fixed = list(wt_pp = "uniform", sd_N = 0)
    # color_vals = NULL; facet_scales = "fixed"; facet_ncol = NULL
    # y_breaks = waiver(); y_max = NULL
    # rm(yvar, type, col_fct, x_facet_fct, y_facet_fct, shp_lty_fct, fixed)
    # rm(color_vals, facet_scales, facet_ncol, dd, facet_fct, facet_form)
    # rm(points_lines, shp_lty_ttl, facets_coord, y_breaks, y_max)

    facet_scales <- match.arg(facet_scales, c("fixed", "free", "free_x", "free_y"))

    dd <- sim_df |>
        filter(type == .env$type) |>
        select(-type, -sim) |>
        split(~ sd_N + virus_attract + pseudo_repel, drop = TRUE) |>
        # Add np = 0 for each landscape type:
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
        list_rbind() |>
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
        shp_lty_ttl <- pretty_params(shp_lty_fct, TRUE) |>
            str_replace_all("\\s+", "<br>")
    } else {
        points_lines <- list(geom_point(), geom_line(aes(group = grp)))
        shp_lty_ttl <- NULL
    }

    if (length(facet_fct) > 0) {
        facet_form <- paste(paste(y_facet_fct, collapse = "+"), "~",
                            paste(x_facet_fct, collapse = "+")) |>
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
        col_scale <- scale_color_manual(pretty_params(col_fct, TRUE, serif = TRUE),
                                        values = color_vals)
    } else {
        col_scale <- scale_color_viridis_d(pretty_params(col_fct, TRUE, serif = TRUE),
                                           begin = 0.2, end = 0.8, option = "plasma")
    }

    dd |>
        ggplot(aes(n_pseudo / 10e3 * 100, .data[[yvar]], color = .data[[col_fct]])) +
        geom_hline(yintercept = 0, color = "gray70") +
        points_lines +
        col_scale +
        scale_linetype_manual(shp_lty_ttl, values = c("solid", "22")) +
        # scale_linewidth_manual(shp_lty_ttl, values = c(1, 0.5)) +
        scale_shape_manual(shp_lty_ttl, values = c(19, 17)) +
        scale_x_continuous(breaks = c(0, 10, 30, 50, 70, 90)) +
        scale_y_continuous(breaks = y_breaks) +
        facets_coord +
        labs(x = "Percent *Pseudomonas* patches",
             y = yvar_desc[[yvar]] |> first_cap(),
             title = sprintf("*Pseudomonas* %s outbreaks",
                             ifelse(type == "low", "inhibits", "promotes"))) +
        theme(legend.text = element_markdown())
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
#' When `type = "high"`
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
#'       effect of *Pseudomonas* on outbreak sizes at 90% *Pseudomonas* patches
#'       compared to 70%
#'         - this is especially pronounced virus starts on *Pseudomonas*  and
#'           when landscape is uniform
#' * `sd_N`:
#'     - `sd_N` > 0 can cause slight increase in outbreak size but the effect
#'       is pretty neglible
#' * `wt_vp` (Pseudomonas location in relation to virus):
#'     - outbreak sizes always lower when virus starts on *Pseudomonas*
#'     - this effect is strongest when `pseudo_repel` = 5.0
#'     - this effect is stronger when in a clustered landscape only
#'       when `virus_attract` = 1
#'     - can cause *Pseudomonas* to have little effect on outbreak size,
#'       or to even slightly decrease them
#' * `wt_pp` (Pseudomonas clustering):
#'     - clustered causes larger outbreaks at high *Pseudomonas* densities
#'     - only occurs when `virus_attract` or `pseudo_repel` equals 1,
#'       with `virus_attract` = 1 having the stronger effect
#'
#'
#'
#' # ------------
#' When `type = "low"`
#' # ------------
#'
#' * `pseudo_repel` (rho):
#'     - `pseudo_repel` = 5 increases outbreak sizes when virus starts off
#'       *Pseudomonas*
#'     - when virus starts on *Pseudomonas*, `pseudo_repel` = 5 decreases
#'       outbreak sizes at low *Pseudomonas*  densities and has little effect
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
#'     - does nothing
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
#' When `type = "high"`
#' # ------------
#' * `pseudo_repel` (rho):
#'     - `pseudo_repel` = 5 decreases chance of outbreak when virus starts on
#'       *Pseudomonas*
#'     - `pseudo_repel` = 5 has little effect when virus starts off *Pseudomonas*
#' * `virus_attract` (nu):
#'     - `virus_attract` > 1 always increases chance of outbreak
#'     - this effect is strongest when virus starts on *Pseudomonas*
#' * `sd_N`:
#'     - `sd_N` > 0 can make outbreaks less likely, mostly at low *Pseudomonas*
#'       densities
#'     - overall effect is weak
#' * `wt_vp` (Pseudomonas location in relation to virus):
#'     - outbreaks far less likely when virus starts on *Pseudomonas*
#'     - effect is strongest when `virus_attract` = 1 and `pseudo_repel` = 5
#' * `wt_pp` (Pseudomonas clustering):
#'     - effect is negligible
#'
#'
#'
#' # ------------
#' When `type = "low"`
#' # ------------
#'
#' * `pseudo_repel` (rho):
#'     - `pseudo_repel` = 5 increases chance of outbreak when virus starts off
#'       *Pseudomonas*, especially at higher *Pseudomonas* densities (>= 50%)
#'     - `pseudo_repel` = 5 decreases chance of outbreak when virus starts
#'       on *Pseudomonas*, but this effect levels off with greater
#'       *Pseudomonas* such that at 90% *Pseudomonas* patches, it is gone
#'       entirely
#' * `virus_attract` (nu):
#'     - greater value always increases chance of outbreak
#'     - effect is strongest when virus starts on *Pseudomonas*
#' * `sd_N`:
#'     - effect is negligible
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






for (.t in c("low", "high")) {
    for (.v in c("outbreak_size", "p_emerge")) {
        .ym <- NULL
        .yb <- c(0, 0.5, 1)
        if (.v == "outbreak_size") {
            .ym <- NA
            .yb <- waiver()
        }
        .pl <- map(sort(unique(sim_df$virus_attract)), \(.virus_attract) {
            big_land_plotter(yvar = .v, type = .t,
                             col_fct = "pseudo_repel",
                             shp_lty_fct = "wt_vp",
                             x_facet_fct = NULL,
                             y_facet_fct = NULL,
                             fixed = list(wt_pp = "uniform", sd_N = 0,
                                          virus_attract = .virus_attract),
                             color_vals = scico(2, end = 0.8, palette = "hawaii"),
                             y_max = .ym, y_breaks = .yb) +
                theme(plot.margin = margin(0, 0, 0, 0))
        })
        .p <- wrap_plots(.pl, nrow = 1, guides = "collect", axis_titles = "collect") &
            illustrator_theme
        .f <- sprintf("_plots/large-plantscapes-%s-%s.pdf", .t, .v)
        save_plot(.f, .p, width = 3, height = 2)
    }
}; rm(.t, .v, .ym, .yb, .pl, .p, .f)





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


big_land_plotter("outbreak_size", "high", col_fct = "wt_pp",
                 shp_lty_fct = "wt_vp",
                 y_facet_fct = "sd_N",
                 x_facet_fct = c("virus_attract", "pseudo_repel"),
                 facet_ncol = 4, facet_scales = "free_y")

# When virus starts off Pseudomonas:
# --------------------*

# clustering has much longer tail with a few outbreaks being really large:
big_off_df <- sim_df |>
    filter(type == "high", sd_N == 0, virus_attract == 1, pseudo_repel == 5,
           n_pseudo == 9000, wt_vp == "off virus")

big_off_df |>
    select(wt_pp, sim) |>
    mutate(outbreak_size = map(sim, \(x) x$outbreak_size[x$outbreak_size > 1])) |>
    select(-sim) |>
    unnest(outbreak_size) |>
    ggplot(aes(outbreak_size, after_stat(density), color = wt_pp)) +
    geom_freqpoly(bins = 25, linewidth = 1) +
    scale_color_manual(values = c(clustered = "dodgerblue", uniform = "gray60"))

# Biggest 9 outbreaks for clustered and uniform landscapes:
big_off_lands <- big_off_df |>
    (\(x) {z <- x$sim; names(z) <- x$wt_pp; return(z)})() |>
    map(\(x) {
        x |>
            arrange(desc(outbreak_size)) |>
            slice_head(n = 9) |>
            getElement("rep")
    }) |>
    (\(x) {
        dd <- read_rds("_scripts/interm-data/large-plantscapes.rds") |>
            filter(type == big_off_df$type[[1]],
                   sd_N == big_off_df$sd_N[[1]],
                   virus_attract == big_off_df$virus_attract[[1]],
                   pseudo_repel == big_off_df$pseudo_repel[[1]],
                   n_pseudo == big_off_df$n_pseudo[[1]],
                   ifelse(wt_vp < 1, "off virus", "on virus") == big_off_df$wt_vp[[1]]) |>
            mutate(wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
                       factor(levels = c("uniform", "clustered")))
        lands <- imap(x, \(xx, i) {
            dd$landscape[[which(dd$wt_pp == i)]][,,xx]
        })
    })()

land_plotter(big_off_lands$uniform)
land_plotter(big_off_lands$clustered)

land_plotter(big_off_lands$uniform[1:8, 1:8,])
land_plotter(big_off_lands$clustered[1:8, 1:8,])



# When virus starts on Pseudomonas:
# --------------------*

# clustering has much longer tail with a few outbreaks being really large:
big_on_df <- sim_df |>
    filter(type == "high", sd_N == 0, virus_attract == 1, pseudo_repel == 5,
           n_pseudo == 9000, wt_vp == "on virus")

big_on_df |>
    select(wt_pp, sim) |>
    mutate(outbreak_size = map(sim, \(x) x$outbreak_size[x$outbreak_size > 1])) |>
    select(-sim) |>
    unnest(outbreak_size) |>
    ggplot(aes(outbreak_size, after_stat(density), color = wt_pp)) +
    geom_freqpoly(bins = 25, linewidth = 1) +
    scale_color_manual(values = c(clustered = "dodgerblue", uniform = "gray60"))

# Landscapes for biggest 2 outbreaks for clustered and uniform landscapes:
big_on_lands <- big_on_df |>
    (\(x) {z <- x$sim; names(z) <- x$wt_pp; return(z)})() |>
    map(\(x) {
        x |>
            arrange(desc(outbreak_size)) |>
            slice_head(n = 2) |>
            getElement("rep")
    }) |>
    (\(x) {
        dd <- read_rds("_scripts/interm-data/large-plantscapes.rds") |>
            filter(type == big_on_df$type[[1]],
                   sd_N == big_on_df$sd_N[[1]],
                   virus_attract == big_on_df$virus_attract[[1]],
                   pseudo_repel == big_on_df$pseudo_repel[[1]],
                   n_pseudo == big_on_df$n_pseudo[[1]],
                   ifelse(wt_vp < 1, "off virus", "on virus") == big_on_df$wt_vp[[1]]) |>
            mutate(wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
                       factor(levels = c("uniform", "clustered")))
        lands <- imap(x, \(xx, i) {
            dd$landscape[[which(dd$wt_pp == i)]][,,xx]
        })
    })()

circleFun <- function(center = c(0,0), radius = 0.5, npoints = 100,
                      rad_min = 0, rad_max = 2 * pi){
    tt <- seq(rad_min, rad_max, length.out = npoints)
    xx <- center[1] + radius * cos(tt)
    yy <- center[2] + radius * sin(tt)
    return(tibble(x = xx, y = yy, z = 1:length(xx)))
}



land_plotter(big_on_lands$uniform) +
    # geom_rect(xmin = 0.5, xmax = 8.5, ymin = -0.5, ymax = -8.5, fill = NA, color = "black")
    geom_path(data = circleFun(c(0.5, 0.5), radius = pop_info$radius,
                               rad_min = 0, rad_max = 0.5 * pi),
              aes(x, y), inherit.aes = FALSE)
land_plotter(big_on_lands$clustered)

land_plotter(big_on_lands$uniform[1:8, 1:8,])
land_plotter(big_on_lands$clustered[1:8, 1:8,])



