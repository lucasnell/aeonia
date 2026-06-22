
#'
#' Density plots for larger landscape simulations.
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
                 out_stages = TRUE)


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
va = 100
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
                labs(title = serify(pretty_params("virus_attract",FALSE,TRUE,TRUE),
                                    paste(" =", va), "")) +
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


if (.overwrite) {
    for (x in names(inf_p_list)) {
        fn <- sprintf("_plots/virus-attract-nu-%03i.pdf", as.integer(x))
        save_plot(fn, inf_p_list[[x]], width = 7.5, height = 7.5)
    }; rm(x, fn)
}





landscape <- array(c(1L, rep(0L, 100^2-1)), c(100, 100, 24))

# Takes ~1 min
set.seed(1011497921)
w1_sims <- crossing(va = c(1, 5, 100)) |>
    pmap(\(wr, va) {
        ni <- large_simmer(landscape, wasp_resp = "moderate", virus_attract = va,
                           p_load = 1,
                           w = 1,
                           total_exp_days = 1,
                           fly_p = 0.1,
                           sd_N = 0, pseudo_repel = 1) |>
            getElement("n_infected")
        boots <- booter(ni)
        tibble(virus_attract = va, n_infected = mean(ni),
               lo  = boots[["Lower"]], hi  = boots[["Upper"]])
    }, .progress = .prog_args) |>
    list_rbind()

w1_sims

#   virus_attract n_infected    lo    hi
#           <dbl>      <dbl> <dbl> <dbl>
# 1             1      3888. 3705. 4084.
# 2             5      2563  2382. 2737.
# 3           100       160   140.  179.
