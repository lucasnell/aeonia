
#'
#' Plots for larger landscape simulations
#'

source("_scripts/00-preamble.R")



# from 03-large-plantscapes.sh:

set.seed(47127361) # for bootstrapping
sim_df <- read_rds("_scripts/interm-data/large-plantscapes.rds") |>
    mutate(outbreak_size = map_dbl(sim, \(x) mean(x$outbreak_size)),
           boots = map(sim, \(x) boot_ci(x$outbreak_size))) |>
    unnest(boots) |>
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
               ((wt_vp == "off virus" & wt_pp == "uniform" & n_pseudo == 3e3) |
                    n_pseudo == 0)) |>
    select(type, n_pseudo, outbreak_size) |>
    arrange(type, n_pseudo)

sim_df |>
    filter(sd_N == 0 & virus_attract == 1 & pseudo_repel == 1 &
               ((wt_vp == "on virus" & wt_pp == "uniform" & n_pseudo == 3e3) |
                    n_pseudo == 0)) |>
    select(type, n_pseudo, outbreak_size) |>
    arrange(type, n_pseudo)

big_land_plotter <- function(type,
                             col_fct = "sd_N",
                             x_facet_fct = "land_type",
                             y_facet_fct = c("virus_attract", "pseudo_repel"),
                             shp_lty_fct = NULL,
                             fixed = NULL,
                             color_vals = NULL) {

    # type = "low"; col_fct = "sd_N"; x_facet_fct = c("virus_attract", "pseudo_repel")
    # y_facet_fct = NULL; shp_lty_fct = NULL; fixed = list(wt_pp = "uniform")
    # rm(type, col_fct, x_facet_fct, y_facet_fct, shp_lty_fct, fixed)
    # rm(dd, facet_fct, facet_form)

    dd <- sim_df |>
        filter(type == .env$type) |>
        select(-type, -landscape, -sim) |>
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
        stopifnot(length(fixed) == 1L)
        dd <- dd[dd[[names(fixed)]] == fixed[[1]],]
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
        shp_lty_ttl <- pretty_params(shp_lty_fct, TRUE) |>
            str_replace_all("\\s+", "<br>")
    } else {
        points_lines <- list(geom_point(), geom_line(aes(group = grp)))
        shp_lty_ttl <- NULL
    }

    facet_form <- paste(paste(y_facet_fct, collapse = "+"), "~",
                        paste(x_facet_fct, collapse = "+")) |>
        as.formula()

    if (!is.null(color_vals)) {
        col_scale <- scale_color_manual(pretty_params(col_fct, TRUE, serif = TRUE),
                                        values = color_vals)
    } else {
        col_scale <- scale_color_viridis_d(pretty_params(col_fct, TRUE, serif = TRUE),
                                           begin = 0.2, end = 0.8, option = "plasma")
    }

    dd |>
        ggplot(aes(n_pseudo / 10e3 * 100, outbreak_size, color = .data[[col_fct]])) +
        geom_hline(yintercept = 0, color = "gray70") +
        points_lines +
        col_scale +
        scale_linetype_manual(shp_lty_ttl, values = c("solid", "22")) +
        scale_shape_manual(shp_lty_ttl, values = c(19, 17)) +
        scale_x_continuous(breaks = c(0, 10, 30, 50, 70, 90)) +
        facet_grid(facet_form) +
        coord_cartesian(ylim = c(0, max(sim_df$outbreak_size))) +
        labs(x = "Percent *Pseudomonas* patches",
             y = "Outbreak size",
             title = sprintf("*Pseudomonas* %s outbreaks",
                             ifelse(type == "low", "inhibits", "promotes"))) +
        theme(legend.text = element_markdown())
}






# big_land_plotter("high", col_fct = "virus_attract",
#                  # shp_lty_fct = "pseudo_repel",
#                  x_facet_fct = "land_type",
#                  y_facet_fct = "pseudo_repel", fixed = list(sd_N = 50))


big_land_plotter(type = "high",
                 col_fct = "sd_N",
                 shp_lty_fct = "wt_vp",
                 x_facet_fct = c("virus_attract", "pseudo_repel"),
                 y_facet_fct = NULL, fixed = list(wt_pp = "uniform"))

big_land_plotter(type = "high",
                 col_fct = "virus_attract",
                 shp_lty_fct = "pseudo_repel",
                 x_facet_fct = c("wt_vp", "sd_N"),
                 y_facet_fct = "wt_pp", # , fixed = list(wt_pp = "uniform"),
                 color_vals = c("black", "red"))

# for (.t in c("low", "high")) {
#     .p <- big_land_plotter(type = .t,
#                            col_fct = "virus_attract",
#                            shp_lty_fct = "pseudo_repel",
#                            x_facet_fct = c("sd_N", "wt_vp"),
#                            y_facet_fct = NULL, fixed = list(wt_pp = "uniform"),
#                            color_vals = scico(2, end = 0.8, palette = "hawaii")) +
#         illustrator_theme +
#         theme(strip.text.x = element_blank())
#     .f <- sprintf("_plots/large-plantscapes-%s.pdf", .t)
#     save_plot(.f, .p, width = 6, height = 2)
# }; rm(.t, .p, .f)


# LEFT OFF ----
# Use above (for both high and low) for figure?



#' Notes:
#'
#' When `type = "high"`
#' - `pseudo_repel`:
#'     - `pseudo_repel` = 5 usually decreases outbreak sizes and the effect
#'       of *Pseudomonas* on outbreak sizes
#'     - increases outbreak sizes and makes *Pseudomonas* more likely to promote
#'       outbreaks only when *Pseudomonas* starts off the virus and is
#'       uniformly distributed; this effect is minor
#'     - effects are minimal when *Pseudomonas* starts off the virus
#'     - when *Pseudomonas* starts on the virus, the reduction on outbreak
#'       sizes is strongest at low *Pseudomonas* densities, especially
#'       when `sd_N` = 50 (see description under `wt_vp`, too)
#' - `virus_attract`:
#'     - `virus_attract` > 1 always increases outbreak sizes
#'     - has little/no effect on the shape of the *Pseudomonas* ~ outbreak
#'       size relationship when *Pseudomonas* starts on the virus and
#'       when *Pseudomonas* density is low
#'     - when `virus_attract` = 1, there is an especially large drop in the
#'       effect of *Pseudomonas* on outbreak sizes at 90% *Pseudomonas* patches
#'       compared to 70%, potentially resulting in lower outbreak sizes
#'       compared to no *Pseudomonas* at all
#'     - this is especially pronounced when `sd_N` = 50 and when
#'       `pseudo_repel` = 1
#' - `sd_N`:
#'     - `sd_N` > 0 always decreases outbreak size overall, EXCEPT when
#'       `virus_attract` = 5 and `pseudo_repel` = 1
#'     - can result in negative effect of *Pseudomonas* on outbreak size,
#'     and usually causes
#'     - has especially strong effect at low *Pseudomonas* density
#'       (`n_pseudo` = 1000), when `virus_attract` = 1, and especially
#'       when `pseudo_repel` = 5
#' - `wt_vp` (Pseudomonas location in relation to virus):
#'     - outbreak sizes always lower when *Pseudomonas* not on starting virus
#'     - this effect only occurred when `pseudo_repel` = 5.0, and
#'       was most pronounced when `sd_N` = 0 was also true
#'        (see description under `pseudo_repel`, too)
#'     - because this only has effects when *Pseudomonas* is present, it
#'       caused a reduction in outbreak sizes at low *Pseudomonas* densities,
#'       often followed by a rebound at higher densities;
#'       this rebound often did not reach the outbreak size for no *Pseudomonas*
#' - `wt_pp` (Pseudomonas clustering):
#'     - no consistent effects
#'
#'
#' When `type = "low"`
#' - `pseudo_repel`:
#'     - has little to no effect when *Pseudomonas* does not start on virus
#'     - when *Pseudomonas* starts on virus, greater `pseudo_repel` causes
#'       *Pseudomonas* to have greater negative effect on outbreak size
#'       at lower densities (i.e., went from linear to convex curve),
#'       although it doesn't change the effect of *Pseudomonas* on outbreak
#'       size when measuring between `n_pseudo` = 0 and 9000
#' - `virus_attract`:
#'     - greater value increases outbreak size significantly
#'     - greater value also increases the magnitude of the effect of *Pseudomonas*
#' - `sd_N`:
#'     - does nothing
#' - `wt_vp` (Pseudomonas location in relation to virus):
#'     - outbreak size reduced when *Pseudomonas* starts on virus (this
#'       only obviously happens when there is *Pseudomonas* on landscape),
#'       resulting in a more convex curve instead of linear
#'     - this effect is most pronounced when `virus_attract` = 1 and
#'       especially when `pseudo_repel` = 5
#'     - when `pseudo_repel` = 1 and `virus_attract` = 5,
#'       this effect is gone entirely
#' - `wt_pp` (Pseudomonas clustering):
#'     - little to no effect
#'
#'
#' From this, I think it makes sense to...
#'
#' use these in main text:
#' - `pseudo_repel`
#' - `virus_attract`
#' - `sd_N`
#' - `wt_vp`
#'
#' and relegate `wt_pp` to the supplement
#'
#'









# ============================================================================*
# Landscape plots ----
# ============================================================================*


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



land_df <- sim_df |>
    filter(n_pseudo > 0) |>
    distinct(n_pseudo, wt_vp, wt_pp, landscape) |>
    mutate(wt_vp2 = str_replace_all(wt_vp, " ", "_")) |>
    arrange(n_pseudo, wt_vp, wt_pp)


land_df |>
    mutate(land_type = sprintf("n<sub>P</sub> = %i<br>%s<br>%s",
                               n_pseudo, wt_vp, wt_pp)) |>
    (\(x) set_names(x[["landscape"]], x[["land_type"]]))() |>
    land_plotter()

# if (!dir.exists("_plots/big-landscapes-xy")) dir.create("_plots/big-landscapes-xy")
# for (i in 1:nrow(land_df)) {
#     fn <- sprintf("_plots/big-landscapes-xy/%i_np-%s-%s.pdf",
#                   land_df$n_pseudo[[i]], land_df$wt_vp2[[i]], land_df$wt_pp[[i]])
#     p <- land_plotter(land_df$landscape[[i]], .expand_axes = FALSE) +
#         illustrator_theme
#     save_plot(fn, p, width = 4, height = 4)
# }; rm(i, fn, p)


