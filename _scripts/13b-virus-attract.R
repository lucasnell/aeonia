
#'
#' Virus attraction plots for larger landscape simulations.
#'

source("_scripts/00-preamble.R")
library(gganimate)
library(magick)


# Plotting and data functions:
source("_scripts/08b-large-plot-funs.R")

.overwrite <- FALSE



# Combine animated plots using magick:
anim_combiner <- function(a, b, fps = 10, width = 400, height = 400,
                          nframes = 101, ...) {

    op <- options(warn = -1)
    on.exit(options(op))

    suppressMessages({a_gif <- animate(plot = a, fps = fps,
                                       width = width, height = height,
                                       renderer = magick_renderer(),
                                       nframes = nframes, ...)})
    suppressMessages({b_gif <- animate(plot = b, fps = fps,
                                       width = width, height = height,
                                       renderer = magick_renderer(),
                                       nframes = nframes, ...)})

    stopifnot(length(a_gif) == length(b_gif))

    new_gif <- image_append(c(a_gif[1], b_gif[1]))
    for(i in 2:length(a_gif)){
        combined <- image_append(c(a_gif[i], b_gif[i]))
        new_gif <- c(new_gif, combined)
    }

    return(new_gif)
}






set.seed(1108067187) # for bootstrapping
va_sim_df <- list.files("_scripts/interm-data/virus-attract", "*.rds",
                        full.names = TRUE) |>
    map(\(x) {
        read_rds(x) |>
            # Remove this vector immediately bc it's quite large
            select(-landscape)
    }) |>
    list_rbind() |>
    mutate(outbreak_size = map_dbl(sim, \(x) mean(x$n_infected[x$n_infected > 1])),
           p_emerge = map_dbl(sim, \(x) mean(x$n_infected > 1)),
           n_infected = map_dbl(sim, \(x) mean(x$n_infected)),
           boots = map(sim, \(x) ci_booter(x$n_infected, "all"))) |>
    mutate(wt_vp = ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") |>
               factor(levels = c("off *Pseudo.*", "on *Pseudo.*")),
           wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
               factor(levels = c("uniform", "clustered")),
           wasp_resp = factor(wasp_resp, levels = levels(wasp_resp_fct))) |>
    add_no_pseudo_points()

va_p <- crossing(fp = c(0.1, 0.5),
         wr = wasp_resp_fct,
         onp = c(FALSE, TRUE)) |>
    # Bc factors behave weird sometimes (but are needed for sorting crossing):
    mutate(across(where(is.factor), paste)) |>
    # Used for tags:
    mutate(tg = LETTERS[1:n()]) |>
    pmap(\(wr, onp, fp, tg) {
        non_defs <- list(p_load = 1)
        if (onp) non_defs <- c(non_defs, list(wt_vp = "on *Pseudo.*"))
        if (fp > 0.1) non_defs <- c(non_defs, list(fly_p = fp))
        p <- baseline_plotter(outcomes = "n_infected", col_fct = "virus_attract",
                              color_vals = c("black",
                                             lighten(par_pal[["virus_attract"]], 0.4),
                                             par_pal[["virus_attract"]]),
                              non_defaults = non_defs,
                              obs_breaks = 0:2 * 5000, obs_max = 10e3,
                              data_df = va_sim_df |> filter(wasp_resp == wr),
                              multiline_col_title = FALSE,
                              p_tag = tg)
        if (wr != levels(wasp_resp_fct)[1] || onp) {
            p <- p + theme(axis.text.y = element_blank())
        }
        if (fp == 0.1) {
            p <- p + theme(axis.text.x = element_blank())
        }
        return(p)
    }) |>
    add_top_labels() |>
    (\(x) {
        fp_labs <- map(c(0.1, 0.5), \(fp) {
            txt <- serify("", sprintf("*p*<sub>fly</sub> = %.01f", fp), "")
            wrap_elements(richtext_grob(txt, gp = gpar(fontsize = 14, lineheight = 0.8)))
        })
        c(x, fp_labs)
    })() |>
    wrap_plots(design = "III#JJJ\nK#L#M#N\nOOOOOOO\nA#B#C#D\nPPPPPPP\nE#F#G#H",
               guides = "collect", axis_titles = "collect",
               widths = c(1, 0.05, 1, 0.2, 1, 0.05, 1),
               heights = c(0.4, 0.3, 0.25, 1, 0.25, 1)) &
    theme(plot.tag.location = "panel",
          plot.tag.position = c(0.1, 1.05),
          legend.position = "top", legend.title.position = "top")

# va_p






# large landscape simmer, just focusing on infections and incoming alates
# through space
large_inf_disp_simmer <- function(n_pseudo,
                                  wasp_resp,
                                  p_load,
                                  fly_p,
                                  virus_attract,
                                  n_sims = 1L) {

    # n_pseudo = 3000L; wasp_resp = "weak"; p_load = 1; fly_p = 0.1
    # virus_attract = 1; n_sims = 1L
    # rm(n_pseudo, wasp_resp, p_load, fly_p, virus_attract, n_sims)
    # rm(landscape, args, sims)

    if (n_pseudo > 0) {
        landscape <- "_scripts/interm-data/large-landscapes.rds" |>
            read_rds() |>
            filter(n_pseudo == .env$n_pseudo,
                   # uniform, virus starts on uninhabited plant:
                   wt_pp == 1, wt_vp < 1) |>
            getElement("landscape") |>
            getElement(1)
        landscape <- landscape[,,1:n_sims,drop=FALSE]
    } else landscape <- array(c(1L, rep(0L, 99999L)), c(100L, 100L, n_sims))

    zeta <- switch(wasp_resp, weak = 0.1, moderate = 0.5, strong = 0.9)

    args <- list(landscape = landscape,
                 virus_attract = virus_attract,
                 fly_p = fly_p,
                 zeta = zeta,
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 sd_N = 0,
                 pseudo_repel = 1,
                 Y0 = 250,
                 N0 = 55,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = dim(landscape)[3],
                 out_dispersals = "in",
                 force_disps = TRUE,
                 summ = "none",
                 out_stages = "two")


    sims <- do.call(big_plantscape, args) |>
        filter(!is.na(x)) |>
        mutate(virus = as.integer(virus),
               aphids = aphids_juv + aphids_adu + alates_juv + alates_adu + parasitized,
               alates = alates_adu) |>
        select(time, x, y, virus, aphids, alates, disps)

    return(sims)

}





# Takes ~1.3 min
set.seed(876826518)
inf_sims <- crossing(wasp_resp = c(levels(wasp_resp_fct), "moderate"),
                     n_pseudo = c(0L, 3000L, 7000L),
                     virus_attract = c(1, 5, 100)) |>
    pmap(\(wasp_resp, n_pseudo, virus_attract) {
        large_inf_disp_simmer(n_pseudo = n_pseudo,
                              wasp_resp = wasp_resp,
                              p_load = 1,
                              fly_p = 0.1,
                              virus_attract = virus_attract,
                              n_sims = 1L) |>
            mutate(wasp_resp = .env$wasp_resp,
                   n_pseudo = .env$n_pseudo,
                   virus_attract = .env$virus_attract)
    }, .progress = .prog_args) |>
    list_rbind() |>
    mutate(wasp_resp = factor(wasp_resp, levels = c(levels(wasp_resp_fct), "moderate")))


MAX_ALATES <- inf_sims |>
    group_by(wasp_resp, n_pseudo, virus_attract, time) |>
    summarize(alates = sum(alates), .groups = "drop") |>
    getElement("alates") |> max()





# Videos ----

one_combined_gif <- function(wasp_resp, n_pseudo, virus_attract, fps = 5, ...) {

    # wasp_resp = "strong"; n_pseudo = 0L; virus_attract = 1; fps = 5
    # rm(wasp_resp, n_pseudo, virus_attract, fps)
    # rm(filt_df, wr, np, va, .title, a, b)

    # Because we're comparing across virus_attract, the point size scale
    # limits are grouped within wasp_resp and n_pseudo:
    max_disps <- inf_sims |>
        filter(wasp_resp == .env$wasp_resp,
               n_pseudo == .env$n_pseudo) |>
        getElement("disps") |>
        max()

    filt_df <- inf_sims |>
        filter(wasp_resp == .env$wasp_resp,
               n_pseudo == .env$n_pseudo,
               virus_attract == .env$virus_attract)

    wr <- paste(filt_df$wasp_resp[[1]])
    np <- filt_df$n_pseudo[[1]]
    va <- filt_df$virus_attract[[1]]
    .title <- serify(paste0(first_cap(wr),", "),
                     sprintf("n<sub>P</sub> = %i, &nu; = %.0f", np, va), "")

    a <- filt_df |>
        filter(virus == 1) |>
        ggplot(aes(x, y)) +
        geom_point(aes(size = disps), color = "red") +
        scale_y_reverse() +
        scale_size_area("Incoming<br>alates", limits = c(0, max_disps), max_size = 4) +
        coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
        labs(title = .title, subtitle = "Time: {frame_time}") +
        transition_time(time) +
        ease_aes("linear")

    # b <- filt_df |>
    #     ggplot(aes(x, y)) +
    #     geom_point(aes(size = alates)) +
    #     scale_y_reverse() +
    #     # scale_size(limits = c(0, 300), range = c(0, 4), guide = "none") +
    #     scale_size_area(max_size = 4) +
    #     coord_cartesian(xlim = c(0, 100), ylim = c(0, 100)) +
    #     labs(title = .title,
    #          subtitle = "Time: {frame_time}") +
    #     transition_time(time) +
    #     ease_aes("linear")

    b <- filt_df |>
        group_by(time) |>
        summarize(alates = sum(alates), .groups = "drop") |>
        ggplot(aes(time, alates)) +
        geom_line(linewidth = 1) +
        geom_point(shape = 21, color = "black", fill = "dodgerblue", size = 5, stroke = 1) +
        coord_cartesian(xlim = c(0, 100), ylim = c(0, MAX_ALATES)) +
        labs(x = "Time", y = "Total adult alates",
             title = .title) +
        transition_reveal(along = time)

    return(anim_combiner(a, b, fps = fps, ...))

}

if (.overwrite || !dir.exists("_plots/virus-attract-gifs")) {

    # Takes ~ 10 min
    inf_gifs <- inf_sims |>
        distinct(wasp_resp, n_pseudo, virus_attract) |>
        mutate(gif = pmap(list(wasp_resp, n_pseudo, virus_attract), one_combined_gif,
                          .progress = .prog_args))
    if (!dir.exists("_plots/virus-attract-gifs")) dir.create("_plots/virus-attract-gifs")
    for (i in 1:nrow(inf_gifs)) {
        x <- sprintf("%s_%04i_%03i", inf_gifs$wasp_resp[[i]], inf_gifs$n_pseudo[[i]],
                     inf_gifs$virus_attract[[i]])
        # fn <- sprintf("_plots/virus-attract-gifs/%s.gif", x)
        # image_write(inf_gifs[[x]], fn, format = "gif")
        fn <- sprintf("_plots/virus-attract-gifs/%s.mp4", x)
        image_write_video(inf_gifs$gif[[i]], fn, framerate = 2)
    }; rm(i, x, fn)

} else {

    inf_gifs <- inf_sims |>
        distinct(wasp_resp, n_pseudo, virus_attract)
    inf_gifs[["gif"]] <- pmap(inf_gifs,
                              \(wasp_resp, n_pseudo, virus_attract) {
                                  x <- sprintf("%s_%04i_%03i", wasp_resp, n_pseudo,
                                               virus_attract)
                                  fn <- sprintf("_plots/virus-attract-gifs/%s.mp4", x)
                                  z <- image_read_video(fn)
                                  return(image_convert(z, "gif"))
                              })

}






# Static figures ----

wr = "moderate"
np = 0L
inf_t = 50L  # time(s) at which new infection explosion occurs


# inf_sims |>
#     filter(wasp_resp == wr, n_pseudo == np) |>
#     filter(time == max(time)) |>
#     group_by(virus_attract) |>
#     summarize(n_infected = sum(virus))

inf_p_list <- c(1, 5) |>
    set_names() |>
    map(
        \(va) {
            max_disps <- inf_sims |>
                filter(wasp_resp == wr, n_pseudo == np, virus_attract < 100) |>
                getElement("disps") |>
                max()
            .title <- serify("No virus attraction (", "<i>&nu;</i> = 1", ")")
            if (va > 1) {
                .title <- .title |>
                    str_replace("No virus", "Virus") |>
                    str_replace("= 1", paste("=", va))
            }
            inf_sims |>
                filter(wasp_resp == wr, n_pseudo == np, virus_attract == va) |>
                filter(virus == 1) |>
                filter(time %in% (inf_t-11L):inf_t) |>
                arrange(disps) |>
                ggplot(aes(x, y)) +
                geom_point(aes(size = disps, color = disps)) +
                scale_y_reverse() +
                scale_size_area("Incoming<br>alates",
                                limits = c(0, max_disps),
                                max_size = 2) +
                scale_color_viridis_c("Incoming<br>alates",
                                      limits = c(0, max_disps),
                                      option = "plasma", direction = -1,
                                      end = 0.9,
                                      guide = guide_colourbar(theme = theme(
                                          # legend.key.size = unit(12, "in"),
                                          legend.key.height = unit(0.24, "in"),
                                          legend.key.width  = unit(1.2, "in")
                                      ))) +
                facet_wrap(~ time, labeller = label_both, nrow = 3, axes = "all") +
                coord_equal(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
                labs(title = .title) +
                # This fixes a very annoying bug where using
                # `legend.text = ggtext::element_markdown(...)` makes the legend
                # really large and doesn't allow me to change its size
                replace_theme(theme_get(),
                              legend.text = element_text(color = "black", size = 9)) +
                theme(legend.position = "top",
                      axis.text.x = element_blank(),
                      axis.text.y = element_blank(),
                      axis.title.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.ticks.y = element_blank())
        })


# inf_p_list[["1"]]
# inf_p_list[["5"]]




infection_plotter <- function(va, np = 0L, wr = "moderate") {

    # va = 1; np = 0L; wr = "moderate"
    # rm(va, np, wr, d_virus_df, inf_t, max_disps, max_bot)
    # rm(.title, spat_p, time_p)

    max_disps <- inf_sims |>
        filter(wasp_resp == wr, n_pseudo == np, virus_attract < 100) |>
        getElement("disps") |>
        max()
    max_bot <- inf_sims |>
        filter(wasp_resp == wr, n_pseudo == np, virus_attract < 100) |>
        group_by(virus_attract, time) |>
        summarize(virus = sum(virus),
                  alates = sum(alates),
                  .groups = "drop") |>
        mutate(d_virus = virus - lag(virus)) |>
        summarize(across(alates:d_virus, \(x) max(x, na.rm = TRUE))) |>
        as.list()
    # alates --> d_virus scale
    f <- function(x) x / max_bot$alates * max_bot$d_virus
    # d_virus --> alates scale
    g <- function(x) x / max_bot$d_virus * max_bot$alates

    d_virus_df <- inf_sims |>
        filter(wasp_resp == wr, n_pseudo == np, virus_attract == va) |>
        select(-wasp_resp, -n_pseudo, -virus_attract) |>
        group_by(time) |>
        summarize(virus = sum(virus),
                  alates = sum(alates),
                  .groups = "drop") |>
        mutate(d_virus = virus - lag(virus),
               alates = f(alates))
    # top three times for max increase in new infected plants
    inf_t <- d_virus_df |>
        mutate(dd_virus = d_virus - lag(d_virus)) |>
        arrange(desc(dd_virus)) |>
        head(n = 3) |>
        getElement("time") |>
        sort()
    .title <- serify("No virus attraction (", "<i>&nu;</i> = 1", ")")
    if (va > 1) {
        .title <- .title |>
            str_replace("No virus", "Virus") |>
            str_replace("= 1", paste("=", va))
    }
    spat_p <- inf_sims |>
        filter(wasp_resp == wr, n_pseudo == np, virus_attract == va) |>
        filter(virus == 1) |>
        filter(time %in% (inf_t-7L)) |>
        arrange(disps) |>
        ggplot(aes(x, y)) +
        geom_point(aes(size = disps, color = disps)) +
        scale_y_reverse() +
        scale_size_area("Incoming<br>alates",
                        limits = c(0, max_disps),
                        max_size = 2) +
        scale_color_viridis_c("Incoming<br>alates",
                              limits = c(0, max_disps),
                              option = "plasma", direction = -1,
                              end = 0.9,
                              guide = guide_colourbar(theme = theme(
                                  # legend.key.size = unit(12, "in"),
                                  legend.key.height = unit(1.2, "in"),
                                  legend.key.width  = unit(0.24, "in")
                              ))) +
        facet_wrap(~ time, labeller = label_both, nrow = 1, axes = "all") +
        coord_equal(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) +
        labs(title = .title) +
        # This fixes a very annoying bug where using
        # `legend.text = ggtext::element_markdown(...)` makes the legend
        # really large and doesn't allow me to change its size
        replace_theme(theme_get(),
                      legend.text = element_text(color = "black", size = 9)) +
        theme(# legend.position = "top",
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              plot.margin = margin(0,0,0,0))
    # spat_p

    time_p <- d_virus_df |>
        ggplot(aes(time, d_virus)) +
        geom_hline(yintercept = 0, color = "gray70") +
        # geom_area(aes(y = alates), fill = "dodgerblue", alpha = 0.5) +
        geom_line(aes(y = alates), color = "dodgerblue", linewidth = 1) +
        geom_line(linewidth = 1, na.rm = TRUE) +
        # geom_vline(xintercept = inf_t -7L, linetype = "solid", color = "gray80",
        #            linewidth = 0.75) +
        geom_vline(xintercept = inf_t, linetype = "22", color = "gray60",
                   linewidth = 0.75) +
        scale_y_continuous(sec.axis = sec_axis(g, "Adult alates")) +
        coord_cartesian(xlim = c(0, 100), ylim = c(0, max_bot$d_virus)) +
        labs(x = "Time (days)", y = "New infected plants") +
        theme(plot.margin = margin(0,0,0,0),
              axis.title.y.right = element_markdown(color = "dodgerblue"),
              axis.text.y.right = element_markdown(color = "dodgerblue"),
              axis.ticks.y.right = element_line(color = "dodgerblue"))
    # time_p
    #
    if (va > 1) {
        time_p <- time_p + theme(axis.text.y = element_blank())
    }
    if (va <= 1) {
        time_p <- time_p + theme(axis.text.y.right = element_blank())
    }

    return(list(spat_p, time_p))
}




inf_p <- map(c(1, 5), infection_plotter) |>
    do.call(what = c) |>
    base::`[`(c(1,3,2,4)) |>
    wrap_plots(design = "A#B\nC#D", axis_titles = "collect", guides = "collect",
               widths = c(1, 0.05, 1)) +
    plot_annotation(tag_levels = "A")

# inf_p




if (.overwrite) {
    save_plot("_plots/virus-attract.pdf", inf_p, width = 7.5, height = 4)
}




# radius ----

Rcpp::cppFunction(
'arma::imat make_neigh_dxdy(const double& radius) {

    typedef arma::uword uint32;
    typedef arma::sword int32;

    // Fill `neigh_dxdy` based on radius:
    uint32 total_rows = 0;
    int32 fl_radius = std::floor(radius);
    std::vector<arma::imat> dxdy_vec;
    dxdy_vec.reserve(fl_radius * 2U + 1U);
    int32 max_dy;
    uint32 rows_x;
    double radius2 = radius * radius;
    for (int32 dx = -fl_radius; dx <= fl_radius; dx++) {
        max_dy = std::floor(std::sqrt(radius2 - static_cast<double>(dx * dx)));
        rows_x = max_dy * 2U;
        if (dx != 0) rows_x++;
        dxdy_vec.emplace_back(rows_x, 2U, arma::fill::none);
        arma::imat& dxdy_i(dxdy_vec.back());
        uint32 i = 0;
        for (int32 dy = -max_dy; dy <= max_dy; dy++) {
            if (dy == 0 && dx == 0) continue;
            dxdy_i.at(i,0) = dx;
            dxdy_i.at(i,1) = dy;
            i++;
        }
        total_rows += dxdy_i.n_rows;
    }

    arma::imat neigh_dxdy(total_rows, 2U);
    uint32 k = 0;
    for (const arma::imat& dxdy : dxdy_vec) {
        for (uint32 j = 0; j < dxdy.n_rows; j++) {
            neigh_dxdy.at(k,0) = dxdy.at(j,0);
            neigh_dxdy.at(k,1) = dxdy.at(j,1);
            k++;
        }
    }

    return neigh_dxdy;
}', depends = "RcppArmadillo")


make_neigh_dxdy(an_environ$radius) |> nrow()
