
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



large_simmer <- function(landscape, type, sd_N, virus_attract, pseudo_repel,
                         outbreaks = "small", ...) {

    if (outbreaks == "big") {
        p_load <- 0.5
        Y0 <- 150
        N0 <- ifelse(type == "low", 150, 20)
    } else {
        p_load <- 0.05
        Y0 <- ifelse(type == "low", 220, 300)
        N0 <- ifelse(type == "low", 60, 35)
    }

    zeta <- ifelse(type == "low", 1, 0.1)

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
            # Remove this vector immediately bc it's >1 GB
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
           outbreaks = factor(outbreaks, levels = c("small", "big")),
           type = factor(type, levels = c("low", "high")))


sim_w1_df <- list.files("_scripts/interm-data", "large-plantscapes-w1-.?.?.rds",
                     full.names = TRUE) |>
    map(\(x) {
        read_rds(x) |>
            # Remove this vector immediately bc it's >1 GB
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
           type = factor(type, levels = c("low", "high")))




#' Basic check to make sure that the two scenarios still do what we expect
#' them to: Pseudomonas decreases outbreak size in "low" and increases in "high"

sim_df |>
    filter(sd_N == 0 & virus_attract == 1 & pseudo_repel == 1 &
               ((wt_vp == "off *Pseudo.*" & wt_pp == "uniform" & n_pseudo == 7e3) |
                    n_pseudo == 0)) |>
    select(type, outbreaks, n_pseudo, outbreak_size) |>
    arrange(type, outbreaks, n_pseudo)
#   type  outbreaks n_pseudo outbreak_size
#   <fct> <fct>        <dbl>         <dbl>
# 1 low   small            0          4.73
# 2 low   small         7000          2.6
# 3 low   big              0       7604.
# 4 low   big           7000       2736.
# 5 high  small            0          2.58
# 6 high  small         7000          5.18
# 7 high  big              0       2015.
# 8 high  big           7000       7745.


#
# Add np = 0 for each landscape type. Used for plotting.
#
add_no_pseudo_points <- function(data_df) {

    stopifnot(all(c("wt_vp", "wt_pp") %in% colnames(data_df)))
    split_cols <- c("type", "outbreaks", "sd_N", "virus_attract", "pseudo_repel") |>
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
                             type,
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
                             w1 = FALSE) {

    # yvar = "outbreak_size"; type = "low"; outbreaks = "small"; col_fct = "pseudo_repel"
    # x_facet_fct = "virus_attract"; y_facet_fct = NULL; shp_lty_fct = "wt_vp"
    # fixed = list(wt_pp = "uniform", sd_N = 0)
    # color_vals = NULL; facet_scales = "fixed"; facet_ncol = NULL
    # y_breaks = waiver(); y_max = NULL; w1 = FALSE
    # rm(yvar, type, outbreaks, col_fct, x_facet_fct, y_facet_fct, shp_lty_fct, fixed)
    # rm(color_vals, facet_scales, facet_ncol, dd, facet_fct, facet_form)
    # rm(points_lines, shp_lty_ttl, facets_coord, y_breaks, y_max, w1)
    # rm(pretty_greek_fct, make_land_type)

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

    if (is.null(.title)) {
        .title <- sprintf("*Pseudomonas*<br>%s outbreaks",
                          ifelse(type == "low", "inhibits", "promotes"))
    }

    facet_scales <- match.arg(facet_scales, c("fixed", "free", "free_x", "free_y"))

    if (w1) {
        dd <- sim_w1_df |>
            filter(type == .env$type) |>
            mutate(sd_N = 0)
    } else dd <- sim_df |>
        filter(type == .env$type, outbreaks == .env$outbreaks) |>
        select(-outbreaks)

    dd <- dd |>
        mutate(outbreak_size = ifelse(is.na(outbreak_size), -Inf, outbreak_size)) |>
        select(-type, -sim) |>
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
        labs(x = "Percent *Pseudomonas* patches",
             y = yvar_desc[[yvar]] |> first_cap(),
             title = .title) +
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


#'
#' Effects of `sd_N` and `wt_pp`:
#'
.y <- "outbreak_size"
.t <- "low"
.o <- "small"
big_land_plotter(yvar = .y, type = .t, outbreaks = .o,
                 col_fct = "wt_pp",
                 shp_lty_fct = "sd_N",
                 x_facet_fct = c("pseudo_repel", "wt_vp"),
                 y_facet_fct = "virus_attract",
                 facet_scales = "free", facet_ncol = 4L)


#'
#' Effects of others:
#'
.y <- "outbreak_size"
.t <- "low"
.o <- "big"
big_land_plotter(yvar = .y, type = .t, outbreaks = .o,
                 col_fct = "pseudo_repel",
                 shp_lty_fct = "wt_vp",
                 x_facet_fct = NULL,
                 y_facet_fct = "virus_attract",
                 fixed = list(wt_pp = "uniform", sd_N = 0),
                 color_vals = scico(2, end = 0.8, palette = "hawaii"),
                 facet_scales = "free")







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
    dens_sims <- crossing(ty = sort(unique(sim_df$type)),
                          np = sort(unique(sim_df$n_pseudo))) |>
        pmap(\(ty, np) {
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
            large_simmer(.landscape, type = ty, sd_N = 0, virus_attract = 1,
                         pseudo_repel = 1, outbreaks = "small",
                         summ = "none", n_sims = dim(.landscape)[3],
                         n_threads = 1L) |>
                mutate(n_pseudo = np, type = ty) |>
                select(type, n_pseudo, everything())
        }) |>
        list_rbind()

    write_rds(dens_sims, rds_files$dens_sims, compress = "gz")

} else {

    dens_sims<- read_rds(rds_files$dens_sims)

}



total_dens_blank <- dens_sims |>
    filter(is.na(x)) |>
    mutate(aphids = aphids + parasitized) |>
    select(type, n_pseudo, time, aphids, alates, wasps) |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    mutate(species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 50)
max_plant_dens_blank <- dens_sims |>
    filter(is.na(x)) |>
    mutate(aphids = aphids + parasitized) |>
    group_by(time) |>
    summarize(across(aphids:wasps, max), .groups = "drop") |>
    select(aphids, alates, wasps) |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    mutate(species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 50)


total_dens_p_list <- levels(sim_df$type) |>
    set_names() |>
    map(\(ty) {
        dens_sims |>
            filter(type == ty) |>
            filter(is.na(x)) |>
            mutate(aphids = aphids + parasitized) |>
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
                 title = sprintf("*Pseudomonas*<br>%s outbreaks",
                                 ifelse(ty == "low", "inhibits", "promotes")))
    })


max_plant_dens_p_list <- levels(sim_df$type) |>
    set_names() |>
    map(\(ty) {
        dens_sims |>
            filter(type == ty) |>
            filter(!is.na(x)) |>
            mutate(n_pseudo = factor(n_pseudo),
                   aphids = aphids + parasitized) |>
            group_by(n_pseudo, time) |>
            summarize(across(aphids:wasps, max), .groups = "drop") |>
            select(n_pseudo, time, aphids, alates, wasps) |>
            pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
            mutate(species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_blank(data = total_dens_blank) +
            geom_line(aes(color = n_pseudo), linewidth = 0.75) +
            scale_color_viridis_d(begin = 0.1, end = 0.9, option = "plasma") +
            facet_wrap(~ species, ncol = 1, scales = "free_y") +
            # coord_cartesian(expand = FALSE, clip = FALSE) +
            labs(x = "Time (days)", y = "Maximum per-plant density",
                 title = sprintf("*Pseudomonas*<br>%s outbreaks",
                                 ifelse(ty == "low", "inhibits", "promotes")))
    })



# wrap_plots(total_dens_p_list, nrow = 1, guides = "collect", axes = "collect")
# wrap_plots(max_plant_dens_p_list, nrow = 1, guides = "collect", axes = "collect")





# =============================================================================*
# =============================================================================*
# Save plots ----
# =============================================================================*
# =============================================================================*

# crossing(.t = c("low", "high"),
#          .v = c(1, 5)) |>
#     pmap(\(.t, .v) {
#     big_land_plotter(yvar = "outbreak_size", type = .t,
#                      w1 = TRUE,
#                      col_fct = "pseudo_repel",
#                      shp_lty_fct = "wt_vp",
#                      x_facet_fct = NULL,
#                      y_facet_fct = NULL,
#                      fixed = list(wt_pp = "uniform", sd_N = 0,
#                                   virus_attract = .v),
#                      color_vals = scico(2, end = 0.8, palette = "hawaii"))
# }) |>
#     do.call(what = wrap_plots) +
#     plot_layout(nrow = 2, axis_titles = "collect", guides = "collect")




if (.overwrite) {
    for (.t in c("low", "high")) {

        .p <- tibble(yvar = c("p_emerge", "outbreak_size")) |>
            mutate(virus_attract = map(1:n(), \(i) sort(unique(sim_df$virus_attract)))) |>
            unnest(virus_attract) |>
            pmap(\(yvar, virus_attract) {
                .ym <- if (yvar == "outbreak_size") NA else NULL
                .yb <- if (yvar == "outbreak_size") waiver() else c(0, 0.5, 1)
                big_land_plotter(yvar = yvar, type = .t, outbreaks = "small",
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
        save_plot(.f, .p, width = 3.25, height = 4.5)
    }; rm(.t, .p, .f)
}






# =============================================================================*
# =============================================================================*
# Big outbreaks  ----
# =============================================================================*
# =============================================================================*


big_outbreaks_p <- c("low", "high") |>
    map(\(type) {
        big_land_plotter(yvar = "outbreak_size", type = type, outbreaks = "big",
                         col_fct = "pseudo_repel",
                         shp_lty_fct = "wt_vp",
                         x_facet_fct = "virus_attract",
                         y_facet_fct = NULL,
                         fixed = list(wt_pp = "uniform", sd_N = 0),
                         color_vals = scico(2, end = 0.8, palette = "hawaii"),
                         y_max = NA) +
            theme(plot.margin = margin(0, 0, 0, 0))
    }) |>
    do.call(what = wrap_plots) +
    plot_layout(design = "A#B", widths = c(1, 0.1, 1),
                axis_titles = "collect", guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))

if (.overwrite) {
    save_plot("_plots/large-plantscapes-big.pdf", big_outbreaks_p,
              width = 7.5, height = 4)
}




# =============================================================================*
# =============================================================================*
# Effect of virus_attract  ----
# =============================================================================*
# =============================================================================*

virus_attract_plotter <- function(type, yvar, outbreaks = "small", w1 = FALSE,
                                  .title = NULL, xvar = "pseudo_percent",
                                  y_lims = NULL, x_lims = NULL,
                                  y_lines = NULL, log_rel = TRUE) {
    # type = "high"; yvar = "outbreak_size"; outbreaks = "big"; w1 = TRUE
    # .title = waiver(); xvar = "log_alates"
    # rm(type, yvar, outbreaks, w1, .title, xvar, y_lims, x_lims, y_lines)
    # rm(y_line, y, y_lab, shp_lty_ttl, dd)

    stopifnot(xvar %in% c("pseudo_percent", "mean", colnames(sim_df$sim[[1]])[-1]))

    type <- paste(type)
    yvar <- paste(yvar)
    if (is.null(.title)) {
        .title <- list(low = "*Pseudomonas*<br>inhibits outbreaks",
                       high = "*Pseudomonas*<br>promotes outbreaks")[[type]]
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
            filter(type == .env$type,
                   wt_pp %in% c(NA, "uniform"))
    } else {
        dd <- sim_df |>
            filter(type == .env$type,
                   outbreaks == .env$outbreaks,
                   wt_pp %in% c(NA, "uniform"),
                   sd_N == 0)
    }

    if (xvar %in% colnames(sim_df$sim[[1]])[-1]) {
        dd[[xvar]] <- map_dbl(dd$sim, \(x) mean(x[[xvar]]))
    }
    x_breaks <- if (xvar == "pseudo_percent") c(0, 10, 30, 50, 70, 90) else {
        waiver()
    }
    x_lab <- if (xvar == "pseudo_percent") {"Percent *Pseudomonas* patches"
    } else if (xvar == "mean") paste("Mean", yvar_desc[[yvar]]) else {
        first_cap(yvar_desc[[xvar]])
    }

    dd |>
        add_no_pseudo_points() |>
        select(wt_vp, n_pseudo, pseudo_repel, virus_attract, any_of(c(yvar, xvar))) |>
        pivot_wider(names_from = "virus_attract", values_from = yvar,
                    id_cols = wt_vp:pseudo_repel, unused_fn = mean) |>
        mutate(rel = `5` / `1`,
               log_rel = log2(`5` / `1`),
               mean = (`5` + `1`) / 2,
               pseudo_percent = n_pseudo / 10e3 * 100) |>
        mutate(pseudo_repel = factor(pseudo_repel)) |>
        ggplot(aes(.data[[xvar]], .data[[y]], color = pseudo_repel)) +
        geom_hline(yintercept = y_lines, color = "gray70") +
        geom_point(aes(shape = wt_vp)) +
        geom_line(aes(linetype = wt_vp)) +
        scale_x_continuous(breaks = x_breaks) +
        labs(x = x_lab,
             y = y_lab,
             title = .title) +
        coord_cartesian(ylim = y_lims, xlim = x_lims) +
        scale_linetype_manual(shp_lty_ttl, values = c("solid", "22")) +
        scale_shape_manual(shp_lty_ttl, values = c(19, 17)) +
        scale_color_manual(pretty_params("pseudo_repel") |>
                               str_replace(" ", "<br>"),
                           values = scico(2, end = 0.8, palette = "hawaii"))
}

virus_attract_p <- crossing(yvar = c("p_emerge", "outbreak_size") |>
                                (\(x) factor(x, levels = x))(),
                            type = sim_df$type |> unique() |> sort()) |>
    mutate(y_lines = 0) |>
    pmap(virus_attract_plotter) |>
    wrap_plots(nrow = 2, guides = "collect", axis_titles = "collect", byrow = TRUE) +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))

big_virus_attract_p <- tibble(type = sim_df$type |> unique() |> sort(),
                              yvar = "outbreak_size",
                              outbreaks = "big") |>
    mutate(y_lines = 0) |>
    pmap(virus_attract_plotter) |>
    wrap_plots(nrow = 1, guides = "collect", axis_titles = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))

# virus_attract_p
# big_virus_attract_p

if (.overwrite) {
    save_plot("_plots/virus_attract-small.pdf", virus_attract_p,
              width = 6.5, height = 6.5)
    save_plot("_plots/virus_attract-big.pdf", big_virus_attract_p,
              width = 6.5, height = 3.5)
}




#
# Effect of virus_attract = 5 when w = 1? ----
#
va_w1_p <- crossing(x = c("pseudo_percent", "mean") |> (\(x) factor(x, levels = x))(),
                    ty = sim_df$type |> unique() |> sort()) |>
    pmap(\(x, ty) {
        .t <- waiver()
        if (x == "pseudo_percent") .t <- NULL
        virus_attract_plotter(yvar = "outbreak_size", type = paste(ty),
                              w1 = TRUE, .title = .t,
                              xvar = paste(x), x_lims = NULL, y_lines = 0)
    }) |>
    do.call(what = wrap_plots) +
    plot_layout(nrow = 2, axis_titles = "collect", guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"))


if (.overwrite) {
    save_plot("_plots/virus_attract-w1.pdf", va_w1_p, width = 6.5, height = 5.5)
}






# =============================================================================*
# =============================================================================*
# Why pseudo_repel = 5 increase outbreaks when virus starts off Pseudomonas? ----
# =============================================================================*
# =============================================================================*

.type <- "high"
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
sims_p1 <- large_simmer(.landscape, .type, .sd_N, .virus_attract,
                        pseudo_repel = 1, n_sims = dim(.landscape)[3], summ = "all",
                        max_t = 50,
                        out_dispersals = TRUE)
sims_p5 <- large_simmer(.landscape, .type, .sd_N, .virus_attract,
                        pseudo_repel = 5, n_sims = dim(.landscape)[3], summ = "all",
                        max_t = 50,
                        out_dispersals = TRUE)


# In the `disps` column the column indicates the plant the alate came from,
# and the row indicates the plant the alate dispersed to.
#
# Below, this is the total number of incoming alates to the initially
# infected patch:

sims_p1$disps |> map_int(\(x) sum(x[1,]))
# [1] 55 52 53 58 54 51 51 64 56 48

sims_p5$disps |> map_int(\(x) sum(x[1,]))
# [1] 132 139 139 132 137 122 132 162 122 122




# =============================================================================*
# =============================================================================*
# Why virus_attract has greatest effect when Pseudomonas at mid densities (at type = high)? ----
# =============================================================================*
# =============================================================================*

.type <- "high"
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

# Takes ~2
set.seed(500860318)
sims_va <- imap(.landscapes, \(l, n) {

    sims_va1 <- large_simmer(l, .type, .sd_N, virus_attract = 1,
                             pseudo_repel = .pseudo_repel, n_sims = dim(l)[3],
                             summ = "all", # max_t = 50,
                             out_dispersals = TRUE) |>
        #
        # In the `disps` column the column indicates the plant the alate came
        # from, and the row indicates the plant the alate dispersed to.
        #
        # The first line below is the total number of incoming alates to the
        # initially infected patch:
        #
        mutate(disps_in = map_int(disps, \(x) sum(x[1,])),
               disps_out = map_int(disps, \(x) sum(x[,1]))) |>
        select(-disps) |>
        mutate(virus_attract = 1)
    sims_va5 <- large_simmer(l, .type, .sd_N, virus_attract = 5,
                             pseudo_repel = .pseudo_repel, n_sims = dim(l)[3],
                             summ = "all", # max_t = 50,
                             out_dispersals = TRUE) |>
        mutate(disps_in = map_int(disps, \(x) sum(x[1,])),
               disps_out = map_int(disps, \(x) sum(x[,1]))) |>
        select(-disps) |>
        mutate(virus_attract = 5)
    bind_rows(sims_va1, sims_va5) |>
        mutate(pseudo = factor(n, levels = c("low", "med", "high"))) |>
        select(pseudo, virus_attract, rep, everything())
}, .progress = .prog_args) |>
    list_rbind()

sims_va |>
    group_by(pseudo, virus_attract) |>
    summarize(across(starts_with("disps"), mean), .groups = "drop")
#   pseudo virus_attract disps_in disps_out
#   <fct>          <dbl>    <dbl>     <dbl>
# 1 low                1     47.5      38.4
# 2 low                5    233.      192.
# 3 med                1     63.2      52.2
# 4 med                5    302.      249.
# 5 high               1     27.3      22.1
# 6 high               5    133.      110.


sims_va |>
    group_by(pseudo, virus_attract) |>
    summarize(across(starts_with("disps"), mean), .groups = "drop") |>
    group_by(pseudo) |>
    summarize(across(starts_with("disps"), \(x) x[virus_attract == 5] / x[virus_attract == 1]))
#   pseudo disps_in disps_out
#   <fct>     <dbl>     <dbl>
# 1 low        4.91      5.01
# 2 med        4.77      4.76
# 3 high       4.89      4.96



sims_va |>
    group_by(pseudo, virus_attract) |>
    summarize(across(c("alates", "log_aphids", "disps_in", "disps_out"), mean), .groups = "drop") |>
    (\(x) {print(x); return(x)})() |>
    group_by(pseudo) |>
    summarize(across(everything(), \(x) x[virus_attract == 5] / x[virus_attract == 1]))



st <- c("low", "med", "high") |>
    map(\(n){
    l <- .landscapes[[n]]
    large_simmer(l, .type, .sd_N, virus_attract = 1,
                 pseudo_repel = .pseudo_repel, n_sims = dim(l)[3],
                 summ = "time") |>
        mutate(pseudo = factor(n, levels = c("low", "med", "high"))) |>
        select(pseudo, everything())
    }) |>
    list_rbind()

st |>
    filter(rep == 1) |>
    mutate(aphids = aphids + parasitized) |>
    select(pseudo, rep, time, aphids, alates, wasps) |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    mutate(id = interaction(pseudo, rep, species, drop = TRUE)) |>
    ggplot(aes(time, density, color = pseudo)) +
    geom_line(aes(group = id), linewidth = 0.75) +
    scale_color_viridis_d(begin = 0.2, end = 0.9) +
    facet_wrap(~ species, ncol = 1, scales = "free")






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
           n_pseudo == 9000, wt_vp == "off *Pseudo.*")

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
                   ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") == big_off_df$wt_vp[[1]]) |>
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
           n_pseudo == 9000, wt_vp == "on *Pseudo.*")

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
                   ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") == big_on_df$wt_vp[[1]]) |>
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



