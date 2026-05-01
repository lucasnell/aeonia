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




large_simmer <- function(landscape, wasp_resp, sd_N, virus_attract, pseudo_repel,
                         outbreaks = "small", ...) {

    N0 <- 55
    p_load <- ifelse(outbreaks == "small", 0.05, 1)
    zeta <- ifelse(wasp_resp == "weak", 0.1, 0.9)
    Y0 <- ifelse(outbreaks == "big_zh", 200, 400)

    args <- list(landscape = landscape,
                 sd_N = sd_N,
                 virus_attract = virus_attract,
                 pseudo_repel = pseudo_repel,
                 Y0 = Y0,
                 N0 = N0,
                 zeta = zeta,
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = 100L,
                 summ = "all")

    args <- list_assign(args, ...)

    return(do.call(big_plantscape, args))

}




# from 06-large-plantscapes.sh:
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
           p_emerge = map_dbl(sim, \(x) mean(x$n_infected > 1))) |>
    mutate(wt_vp = ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") |>
               factor(levels = c("off *Pseudo.*", "on *Pseudo.*")),
           wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
               factor(levels = c("uniform", "clustered")),
           outbreaks = factor(outbreaks, levels = c("small", "big_zh", "big_zl")),
           wasp_resp = factor(wasp_resp, levels = c("strong", "weak")))


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
           wasp_resp = factor(wasp_resp, levels = c("strong", "weak")))


# sim_v100_df <- list.files("_scripts/interm-data", "large-plantscapes-v100-.?.?.rds",
#                      full.names = TRUE) |>
#     map(\(x) {
#         read_rds(x) |>
#             # Remove this vector immediately bc it's >1 GB
#             select(-landscape)
#     }) |>
#     list_rbind() |>
#     mutate(outbreak_size = map_dbl(sim, \(x) mean(x$n_infected[x$n_infected > 1])),
#            log_outbreak_size = map_dbl(sim, \(x) mean(log10(x$n_infected[x$n_infected > 1]))),
#            p_emerge = map_dbl(sim, \(x) mean(x$n_infected > 1))) |>
#     mutate(wt_vp = ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") |>
#                factor(levels = c("off *Pseudo.*", "on *Pseudo.*")),
#            wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
#                factor(levels = c("uniform", "clustered")),
#            wasp_resp = factor(wasp_resp, levels = c("strong", "weak")))




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
#  1 strong    small            0          2.67
#  2 strong    small         7000          2.31
#  3 strong    big_zh           0       6355.
#  4 strong    big_zh        7000       2412.
#  5 strong    big_zl           0       2220.
#  6 strong    big_zl        7000        606.
#  7 weak      small            0          3.03
#  8 weak      small         7000          4.40
#  9 weak      big_zh           0       6338.
# 10 weak      big_zh        7000       9337.
# 11 weak      big_zl           0       2370.
# 12 weak      big_zl        7000       7363.


#
# Add np = 0 for each landscape type. Used inside `big_land_plotter`
# and `virus_attract_plotter`.
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
# Density sims ----
# =============================================================================*
# =============================================================================*

if (! file.exists(rds_files$dens_sims) || .overwrite) {
    # Takes ~26 sec
    set.seed(314679353)
    dens_sims <- crossing(wr = sort(unique(sim_df$wasp_resp)),
                          np = sort(unique(sim_df$n_pseudo))) |>
        pmap(\(wr, np) {
            if (np > 0) {
                .landscape <- read_rds("_scripts/interm-data/large-landscapes.rds") |>
                    filter(n_pseudo == np,
                           wt_pp == 1,
                           ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") ==
                               "off *Pseudo.*") |>
                    getElement("landscape") |> getElement(1) |>
                    base::`[`(,,1,drop=FALSE)
            } else {
                .landscape <- array(c(1L, rep(0L, 100L*100L-1L)), c(100, 100, 1))
            }
            large_simmer(.landscape, wasp_resp = wr, sd_N = 0, virus_attract = 1,
                         pseudo_repel = 1, outbreaks = "small",
                         summ = "none", out_stages = TRUE,
                         n_sims = dim(.landscape)[3]) |>
                mutate(aphids = aphids_juv + aphids_adu + alates_juv + alates_adu + parasitized,
                       alates = alates_adu,
                       n_pseudo = np, wasp_resp = wr) |>
                select(wasp_resp, n_pseudo, everything())
        }) |>
        list_rbind()

    write_rds(dens_sims, rds_files$dens_sims, compress = "gz")

} else {

    dens_sims<- read_rds(rds_files$dens_sims)

}



total_dens_blank <- dens_sims |>
    filter(is.na(x)) |>
    select(wasp_resp, n_pseudo, time, aphids, alates, wasps) |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    mutate(species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 50)
max_plant_dens_blank <- dens_sims |>
    filter(!is.na(x)) |>
    select(time, aphids, alates, wasps) |>
    group_by(time) |>
    summarize(across(aphids:wasps, max), .groups = "drop") |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    mutate(species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 50)


total_dens_p_list <- levels(sim_df$wasp_resp) |>
    set_names() |>
    map(\(wr) {
        dens_sims |>
            filter(wasp_resp == wr) |>
            filter(is.na(x)) |>
            select(n_pseudo, time, aphids, alates, wasps) |>
            pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
            mutate(species = factor(species, levels = c("aphids", "alates", "wasps")),
                   n_pseudo = factor(n_pseudo)) |>
            ggplot(aes(time, density / 1e6)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_blank(data = total_dens_blank) +
            geom_line(aes(color = n_pseudo), linewidth = 0.75) +
            scale_color_viridis_d(begin = 0.1, end = 0.9, option = "plasma") +
            scale_y_continuous(n.breaks = 4L) +
            facet_wrap(~ species, ncol = 1, scales = "free_y") +
            # coord_cartesian(expand = FALSE, clip = FALSE) +
            labs(x = "Time (days)", y = "Total density (&times; 10<sup>6</sup>)",
                 title = scenario_title(wr, TRUE))

    })





max_plant_dens_p_list <- levels(sim_df$wasp_resp) |>
    set_names() |>
    map(\(wr) {
        dens_sims |>
            filter(wasp_resp == wr) |>
            filter(!is.na(x)) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            select(n_pseudo, time, aphids, alates, wasps) |>
            group_by(n_pseudo, time) |>
            summarize(across(aphids:wasps, max), .groups = "drop") |>
            pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
            mutate(species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_blank(data = max_plant_dens_blank) +
            geom_line(aes(color = n_pseudo), linewidth = 0.75) +
            scale_color_viridis_d(begin = 0.1, end = 0.9, option = "plasma") +
            facet_wrap(~ species, ncol = 1, scales = "free_y") +
            # coord_cartesian(expand = FALSE, clip = FALSE) +
            labs(x = "Time (days)", y = "Maximum per-plant density",
                 title = scenario_title(wr, TRUE))

    })



# wrap_plots(total_dens_p_list, nrow = 1, guides = "collect", axes = "collect")
# wrap_plots(max_plant_dens_p_list, nrow = 1, guides = "collect", axes = "collect")







# =============================================================================*
# =============================================================================*
# Save plots ----
# =============================================================================*
# =============================================================================*

# crossing(.t = c("strong", "weak"),
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

    .p <- total_dens_p_list[[1]] +
        (total_dens_p_list[[2]] + theme(axis.text.y = element_blank())) +
        plot_layout(design = "A#B", widths = c(1, 0.05, 1),
                    axis_titles = "collect", guides = "collect") &
        illustrator_theme &
        theme(panel.spacing.y = unit(0.5, "lines"),
              plot.margin = margin(0,0,0,0))
    save_plot("_plots/densities-large-plantscapes.pdf", .p,
              width = 6.3, height = 2)

    for (.t in c("strong", "weak")) {

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


big_outbreaks_p <- crossing(.wasp_resp = c("strong", "weak"),
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
                            wasp_resp = sim_df$wasp_resp |> unique() |> sort()) |>
    mutate(.title = map(yvar, \(y) {
        if (y == "p_emerge") return(NULL) else return(waiver())})) |>
    mutate(y_lines = 0) |>
    pmap(virus_attract_plotter) |>
    wrap_plots(nrow = 2, guides = "collect", axis_titles = "collect", byrow = TRUE) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))

big_virus_attract_p <- crossing(wasp_resp = sim_df$wasp_resp |> unique() |> sort(),
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

va_w1_obs_plots <- map(c("strong", "weak"), \(.t) {
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
                            ty = sim_df$wasp_resp |> unique() |> sort()) |>
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

va_v100_obs_plots <- map(c("strong", "weak"), \(.t) {
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
                            ty = sim_df$wasp_resp |> unique() |> sort()) |>
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



