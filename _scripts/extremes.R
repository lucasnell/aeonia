
#'
#' Basic time series simulations
#'

source("_scripts/00-preamble.R")


# ============================================================================*
# Larger sims ----
# ============================================================================*

#'
#' Do a bunch of simulations for sets of parameter values to show how these
#' parameter sets behave over many simulations.
#'



#' Run lil_landscape under two simulation scenarios: large or small,
#' and under two parameter combo types:
#' "low" (Pseudomonas decreases outbreak size)
#' "high" (Pseudomonas increases outbreak size)
#'
run_sim_combos <- function(type, n_pseudo, large_sims = FALSE, ...) {

    stopifnot(length(type) == 1L && type %in% c("low", "high"))
    stopifnot(length(n_pseudo) == 1L && is.numeric(n_pseudo) && n_pseudo >= 0)
    stopifnot(length(large_sims) == 1L && is.logical(large_sims))

    shared_args <- list(Y0 = 2,
                        sd_N = 0,
                        K = 12.5e3,
                        pseudo_surv = 0.85,
                        n_pseudo = n_pseudo)

    if (large_sims) {
        size_args <- list_assign(shared_args, n_sims = 1000,
                                 spat_config = "random", summ = "all")
    } else {
        size_args <- list_assign(shared_args, n_sims = 1,
                                 spat_config = "diagonal", summ = "none")
    }

    if (type == "high") {
        args <- list_assign(size_args,
                            N0 = 15,
                            virus_attract = 1.5,
                            pseudo_repel = 1.5,
                            zeta = 0.1)
    } else {
        args <- list_assign(size_args,
                            N0 = 90,
                            virus_attract = 4,
                            pseudo_repel = 10,
                            zeta = 1)
    }

    args <- list_assign(args, ...)

    return(do.call(lil_plantscape, args))
}



set.seed(259619622)
large_sims <- crossing(type = factor(1:2, labels = c("low", "high")),
                       n_pseudo = c(0L, 3L)) |>
    mutate(outbreak_size = map2(type, n_pseudo, \(type, n_pseudo) {
        run_sim_combos(type, n_pseudo, TRUE) |>
            getElement("outbreak_size")
    }))

large_sims |>
    mutate(outbreak_size = num(map_dbl(outbreak_size, mean), digits = 3))
# # A tibble: 4 × 3
#   type  n_pseudo outbreak_size
#   <fct>    <int>     <num:.3!>
# 1 low          0         5.680
# 2 low          3         3.072
# 3 high         0         2.790
# 4 high         3         5.589

large_sims |>
    mutate(outbreak_size = map_dbl(outbreak_size, mean)) |>
    group_by(type) |>
    summarize(outbreak_size = outbreak_size[n_pseudo != 0] -
                  outbreak_size[n_pseudo == 0]) |>
    mutate(outbreak_size = num(outbreak_size, digits = 3))
# # A tibble: 2 × 2
#   type  outbreak_size
#   <fct>     <num:.3!>
# 1 low          -2.608
# 2 high          2.799



# ============================================================================*
# Histograms ----
# ============================================================================*




low_high_hist_list <- large_sims |>
    split(~ type) |>
    map(\(d) {
        dd <- d |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            unnest(outbreak_size) |>
            group_by(n_pseudo, outbreak_size) |>
            summarize(perc = 100 * n() / 1000, .groups = "drop")
        # dd_means <- dd |>
        #     group_by(n_pseudo) |>
        #     summarize(outbreak_size = sum((perc / 100) * outbreak_size)) |>
        #     mutate(perc = max(dd$perc) * 0.95)
        dd |>
            ggplot(aes(outbreak_size, perc)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_col(aes(fill = n_pseudo), position = "dodge", width = 0.5) +
            # geom_segment(data = dd_means, aes(color = n_pseudo, yend = perc * 1.1), linewidth = 1.5) +
            scale_fill_manual(values = np_pal, aesthetics = c("color", "fill")) +
            scale_x_continuous(breaks = (0:4) * 2 + 1) +
            labs(x = "Outbreak size", y = "Percent of simulations")
    })

wrap_plots(low_high_hist_list, ncol = 1, guides = "collect", axis_titles = "collect")


for (n in names(low_high_hist_list)) {
    save_plot(sprintf("_plots/extremes-histograms-%s.pdf", n),
              low_high_hist_list[[n]] + illustrator_theme, width = 4.5, height = 1)
}; rm(n)





# ============================================================================*
# Time series ----
# ============================================================================*


set.seed(545888412)
low_high_sims <- large_sims |>
    select(-outbreak_size) |>
    mutate(sims = map2(type, n_pseudo, \(type, n_pseudo) {
                run_sim_combos(type, n_pseudo) |>
                    mutate(plant = interaction(y, x),
                           aphids = aphids + parasitized) |>
                    select(plant, time, aphids, alates, wasps) |>
                    pivot_longer(aphids:wasps, names_to = "species",
                                 values_to = "density") |>
                    mutate(species = factor(species, levels = c("aphids", "wasps",
                                                                "alates")))
    })) |>
    unnest(sims)



# # To see that only 3 plants (P- plant & P+ landscape, P+ plant & P+ landscape,
# # P- plant & P- landscape) can represent all in the simulations,
# # run the following:
# low_high_sims |>
#     split(~ type + n_pseudo) |>
#     map(\(d) {
#         d |>
#             filter(!is.na(plant)) |>
#             split(~ species) |>
#             map(\(dd) {
#                 dd |>
#                     ggplot(aes(time, density)) +
#                     geom_line(aes(color = species), linewidth = 1) +
#                     facet_wrap( ~ plant, nrow = 1) +
#                     scale_color_manual(values = spp_pal, guide = "none") +
#                     labs(x = "Time (days)", y = dd[["species"]][[1]])
#             }) |>
#             do.call(what = wrap_plots) +
#             plot_layout(ncol = 1) +
#             plot_annotation(title = sprintf("type = %s, n<sub>P</sub> = %s",
#                                             d$type[[1]], d$n_pseudo[[1]]),
#                             theme = theme(plot.title = element_markdown()))
#     })



# -------------------------------------*
# ... totals ----
# -------------------------------------*

# This is used to get the same axis max values across plots:
tot_empty_data <- tibble(species = c("aphids", "wasps", "alates"),
                         time = 0) |>
    mutate(species = factor(species, levels = species),
           density = map_dbl(species, \(spp) {
               low_high_sims |>
                   filter(species == spp, is.na(plant)) |>
                   getElement("density") |>
                   max()
           }))


tot_low_high_p_list <- low_high_sims |>
    split(~ type) |>
    map(\(x) {
        x |>
            filter(is.na(plant)) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(aes(color = n_pseudo), linewidth = 1) +
            geom_blank(data = tot_empty_data) +
            scale_color_manual(values = np_pal) +
            facet_wrap(~ species, scales = "free")
    })

wrap_plots(tot_low_high_p_list, ncol = 1, guides = "collect", axis_titles = "collect")


# for (n in names(tot_low_high_p_list)) {
#     save_plot(sprintf("_plots/extremes-timeseries-totals-%s.pdf", n),
#               tot_low_high_p_list[[n]] + illustrator_theme, width = 5, height = 1.5)
# }; rm(n)





# -------------------------------------*
# ... by plant ----
# -------------------------------------*

# This is used to get the same axis max values across plots:
empty_data <- tibble(species = c("aphids", "wasps", "alates"),
                     time = 0) |>
    mutate(species = factor(species, levels = species),
           density = map_dbl(species, \(spp) {
               low_high_sims |>
                   filter(species == spp, !is.na(plant)) |>
                   getElement("density") |>
                   max()
           }))

low_high_p_list <- low_high_sims |>
    split(~ type) |>
    map(\(x) {
        x |>
            filter((n_pseudo == 0 &  plant == paste(1, 1, sep = ".")) |
                       (n_pseudo == 3 & plant %in% paste(1, 2:3, sep = "."))) |>
            mutate(plant = factor(plant, labels = c("P_none", "P_land", "P_landplant"))) |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(aes(color = plant), linewidth = 1) +
            geom_blank(data = empty_data) +
            scale_color_manual(values = pt_pal) +
            facet_wrap(~ species, scales = "free")
    })

wrap_plots(low_high_p_list, ncol = 1, guides = "collect", axis_titles = "collect")


# for (n in names(low_high_p_list)) {
#     save_plot(sprintf("_plots/extremes-timeseries-%s.pdf", n),
#               low_high_p_list[[n]] + illustrator_theme, width = 5, height = 1.5)
# }; rm(n)












