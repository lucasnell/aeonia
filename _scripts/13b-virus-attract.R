
#'
#' Virus attraction plots for larger landscape simulations.
#'

source("_scripts/00-preamble.R")

.overwrite <- FALSE





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

#
# NOTE: check to make sure the legends look okay (not weirdly large).
# Creating `inf_p` again fixes this.
#
# inf_p




if (.overwrite) {
    save_plot("_plots/virus-attract.pdf", inf_p, width = 7.5, height = 4)
}

