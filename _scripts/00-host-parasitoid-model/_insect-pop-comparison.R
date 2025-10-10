
#'
#' Compare aeonia vs gameofclones host-parasitoid simulations
#'

suppressPackageStartupMessages({
    library(tidyverse)
    library(gameofclones)
    library(aeonia)
    library(patchwork)
    library(ggtext)
})

spp_pal <- c(aphids = "dodgerblue", wasps = "firebrick")

# Susceptible aphid line from `gameofclones`:
line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,10), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")

# Leslie matrix (for apterous, non-parasitized):
L0 <- clonal_line("susceptible",
                  density_0 = cbind(c(0,0,0,0,10), rep(0, 5)),
                  surv_juv_apterous = "high",
                  surv_adult_apterous = "high",
                  repro_apterous = "high",
                  p_instar_smooth = 0) |>
    getElement("leslie") |>
    base::`[`(,,1)
# Stable age distribution for this Leslie matrix:
X0 <- gameofclones:::sad_leslie(L0)


run_goc <- function(aphids0, wasps0, max_t = 100) {

    .line_s <- line_s
    .line_s$density_0 <- cbind(X0, 0) * aphids0
    sim_experiments(.line_s, 1, max_t = max_t, alate_b0 = -Inf, alate_b1 = 0,
                    wasp_density_0 = wasps0, wasp_delay = 0, extinct_N = 0) |>
        (\(sims) {
            left_join(sims$aphids |>
                          filter(type != "mummy", type != "parasitized") |>
                          group_by(time) |>
                          summarize(aphids = sum(N)),
                      sims$wasps |>
                          group_by(time) |>
                          summarize(wasps = sum(wasps)),
                      by = join_by(time))
        })() |>
        pivot_longer(aphids:wasps, names_to = "species", values_to = "density")
}



run_aeonia <- function(aphids0, wasps0, max_t = 100) {
    test_insect_pops(max_t = max_t, N0 = aphids0, W0 = 0, Y0 = wasps0, B = 0,
                     alate_0 = -Inf, alate_1 = 0) |>
        mutate(aphids = aphids + alates) |>
        select(time, aphids, wasps) |>
        pivot_longer(-time, names_to = "species", values_to = "density")
}




combos <- crossing(a0 = 10^(0:2), w0 = 10^(-3:-1)) |>
    mutate(w0 = a0 * w0)


sim_goc <- pmap(combos, \(a0, w0) {
    run_goc(aphids0 = a0, wasps0 = w0) |>
        mutate(combo = paste(a0, w0, sep = "_"))
}) |>
    list_rbind() |>
    mutate(combo = factor(combo),
           sim_type = "gameofclones")

sims_aeonia <- pmap(combos, \(a0, w0) {
    run_aeonia(aphids0 = a0, wasps0 = w0) |>
        mutate(combo = paste(a0, w0, sep = "_"))
}) |>
    list_rbind() |>
    mutate(combo = factor(combo),
           sim_type = "aeonia")


sim_ts_plotter <- function(.df, .title = waiver()) {

    sim_types <- unique(.df$sim_type)
    if (length(sim_types) > 2 || ! "gameofclones" %in% sim_types) {
        stop("`sim_type` should have two unique values, one of which is 'gameofclones'")
    }
    other_type <- sim_types[sim_types != "gameofclones"]

    plot_list <- .df |>
        mutate(combo = str_split(combo, "_") |>
                   map_chr(\(x) {
                       sprintf("aphids<sub>0</sub> = %s<br>wasps<sub>0</sub> = %s",
                               x[1], x[2])
                   }),
               sim_type = factor(sim_type, levels = c("gameofclones", other_type))) |>
        split(~ combo) |>
        map(\(dd) {
            .mult <- max(dd$density[dd$species == "aphids"]) /
                max(dd$density[dd$species == "wasps"])


            dd |>
                mutate(density = ifelse(species == "wasps",
                                        density * .mult, density)) |>
                ggplot(aes(time, density, color = species, linetype = sim_type)) +
                geom_hline(yintercept = 0, linewidth = 0.5, color = "gray60") +
                geom_line(linewidth = 0.75) +
                scale_linetype_manual("Model:", values = c("solid", "22")) +
                scale_y_continuous("Aphid density",
                                   sec.axis = sec_axis(\(x) x / .mult,
                                                       "Wasp density")) +
                scale_color_manual(NULL, values = spp_pal, guide = "none") +
                facet_wrap( ~ combo) +
                theme(plot.title = element_markdown(size = 9),
                      strip.text = element_markdown(size = 9),
                      axis.text = element_text(size = 8),
                      axis.text.y.left = element_text(color = spp_pal[["aphids"]]),
                      axis.text.y.right = element_text(color = spp_pal[["wasps"]]),
                      axis.title = element_text(size = 10),
                      axis.title.y.left = element_text(color = spp_pal[["aphids"]]),
                      axis.title.y.right = element_text(color = spp_pal[["wasps"]]),
                      axis.ticks.y.left = element_line(color = spp_pal[["aphids"]]),
                      axis.ticks.y.right = element_line(color = spp_pal[["wasps"]]))
        })

    plot_list |>
        c(list(ncol = 3, guides = "collect", axis_titles = "collect")) |>
        do.call(what = wrap_plots) +
        plot_annotation(title = .title) &
        theme(legend.position = "top")
}


comp_p <- bind_rows(sim_goc, sims_aeonia) |>
    sim_ts_plotter(.title = "aeonia vs gameofclones")

save_plot("_plots/host-paras-model-comp.pdf", comp_p, width = 6, height = 6)
save_plot("_plots/host-paras-model-comp.svg", comp_p, width = 6, height = 6)


