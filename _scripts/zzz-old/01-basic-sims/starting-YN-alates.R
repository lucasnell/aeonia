
#'
#' Heatmap showing mean(log10(alates + 1)) for various starting densities of
#' aphids and parasitoids.
#'


source("_scripts/00-preamble.R")
source("_scripts/01-basic-sims/00-time-series-plotter.R")


# sobol_summs <- read_rds("_scripts/interm-data/sobol-sims-summs.rds")
#
# sobol_summs |>
#     filter(alate_dens == 1) |>
#     mutate(n_pseudo = factor(n_pseudo)) |>
#     ggplot(aes(log_alates, outbreak_size)) +
#     geom_point(aes(color = n_pseudo), alpha = 0.05) +
#     scale_color_manual(pretty_params("n_pseudo", TRUE), values = np_pal) +
#     labs(x = "mean(log<sub>10</sub>(alate + 1))", y = "Outbreak size")


alate_bounds <- c(1.7, 2.2)

ob_range_p <- sobol_summs |>
    filter(alate_dens == 1) |>
    mutate(log_alates = cut(log_alates, breaks = 100),
           n_pseudo = factor(n_pseudo)) |>
    group_by(n_pseudo, log_alates) |>
    summarize(ob_min = min(outbreak_size),
              ob_max = max(outbreak_size), .groups = "drop") |>
    mutate(log_alates = log_alates |> paste() |> str_remove_all("\\(|\\]") |>
               str_split(",") |> map_dbl(\(x) mean(as.numeric(x)))) |>
    ggplot(aes(log_alates, fill = n_pseudo, color = n_pseudo)) +
    geom_hline(yintercept = c(1, 9), color ="gray70") +
    geom_ribbon(aes(ymin = ob_min, ymax = ob_max), alpha = 0.25, linewidth = 0.5) +
    geom_vline(xintercept = alate_bounds, linetype = "22", linewidth = 0.75) +
    scale_fill_manual(pretty_params("n_pseudo", TRUE), values = np_pal,
                      aesthetics = c("color", "fill")) +
    labs(x = "mean(log<sub>10</sub>(alate + 1))", y = "Outbreak size")
# save_plot("_plots/outbreak_size-ranges.pdf", ob_range_p, width = 6, height = 4)

ob_mid_prob_p <- sobol_summs |>
    filter(alate_dens == 1) |>
    mutate(log_alates = cut(log_alates, breaks = 50),
           n_pseudo = factor(n_pseudo)) |>
    group_by(n_pseudo, log_alates) |>
    summarize(p_middle = mean(outbreak_size > 2 & outbreak_size < 8), .groups = "drop") |>
    mutate(log_alates = log_alates |> paste() |> str_remove_all("\\(|\\]") |>
               str_split(",") |> map_dbl(\(x) mean(as.numeric(x)))) |>
    ggplot(aes(log_alates, p_middle, color = n_pseudo)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = alate_bounds, linetype = "22", linewidth = 0.75) +
    scale_color_manual(pretty_params("n_pseudo", TRUE), values = np_pal) +
    labs(x = "mean(log<sub>10</sub>(alate + 1))", y = "Prob. medium outbreak size")
# save_plot("_plots/outbreak_size-mid-probs.pdf", ob_mid_prob_p, width = 6, height = 4)



alate_simmer <- function(N0, Y0, K) {
    insect_ptr <- make_insect_ptr(K = K, pseudo_surv = 1, fly_p = 0, zeta = 0)
    sims <- sim_plantscape(landscapes = array(0L, c(3L, 3L, 1L)),
                           max_t = 100L,
                           insect_ptr = insect_ptr,
                           N0 = array(N0, c(3L, 3L, 1L)),
                           Y0 = Y0,
                           W0 = array(0, c(3L, 3L, 1L)),
                           virus_attract = 1,
                           pseudo_repel = 1,
                           epsilon = 1,
                           p_load_alate = 0.5,
                           p_load_plant = 0.5,
                           radius = 1,
                           infect_stop = FALSE,
                           summ = "all")
    return(tibble(N0 = N0, Y0 = Y0, K = K, log_alates = sims$log_alates))
}


n0_sim_list <- c(1L, 10L, 100L) |>
    set_names() |>
    map(\(.n0) {
        sim_plantscape(landscapes = array(0L, c(3L, 3L, 1L)),
                       max_t = 100L,
                       insect_ptr = make_insect_ptr(pseudo_surv = 1, fly_p = 0, zeta = 0),
                       N0 = array(.n0, c(3L, 3L, 1L)),
                       Y0 = 9,
                       W0 = array(0, c(3L, 3L, 1L)),
                       virus_attract = 1,
                       pseudo_repel = 1,
                       epsilon = 1,
                       p_load_alate = 0.5,
                       p_load_plant = 0.5,
                       radius = 1,
                       infect_stop = FALSE,
                       summ = "none")
    })

n0_p_list <- n0_sim_list |>
    imap(\(x, i) sim_plotter(x, zeta = 0, .title = sprintf("N0 = %s", i)))

n0_p_list[[1]]
n0_p_list[[2]]
n0_p_list[[3]]




# Takes ~45 sec
alates_sims <- crossing(N0 = 1:100,
                       Y0 = seq(1, 20, 0.2) |> round(1),
                       K = 12500 * c(0.5, 1, 2)) |>
    pmap(alate_simmer, .progress = .prog_args) |>
    list_rbind()



alates_heatmap <- alates_sims |>
    mutate(K = factor(K, levels = sort(unique(K)),
                      labels = sprintf("<span style=\"font-family: serif\">K</span> = %s",
                                       prettyNum(sort(unique(K)), ",")))) |>
    ggplot(aes(N0, Y0)) +
    geom_raster(aes(fill = log_alates)) +
    # geom_contour(aes(z = log_alates), color = "white") +
    geom_contour(aes(z = log_alates), color = "black", breaks = alate_bounds,
                 linetype = "22", linewidth = 0.75) +
    scale_fill_scico("mean(log<sub>10</sub>(alates + 1))",
                     palette = "vik", midpoint = mean(alate_bounds)) +
    labs(x = pretty_params("N0", cap1 = TRUE),
         y = pretty_params("Y0", cap1 = TRUE)) +
    facet_grid( ~ K) +
    theme(legend.position = "top")

# save_plot("_plots/alates-heatmap.pdf", alates_heatmap, width = 8, height = 4)


alates_sims |>
    ggplot(aes(N0, log_alates, color = Y0, group = interaction(Y0, K))) +
    geom_line() +
    # geom_vline(xintercept = c(1, 10, 100)) +
    labs(x = serify("", "N<sub>0</sub>", ""),
         y = "mean(log<sub>10</sub>(alates + 1))") +
    scale_color_scico(serify("", "Y<sub>0</sub>", ""), palette = "batlow") +
    facet_grid( ~ K)

alates_sims |>
    ggplot(aes(Y0, log_alates, color = N0, group = interaction(N0, K))) +
    geom_line() +
    labs(x = serify("", "Y<sub>0</sub>", ""),
         y = "mean(log<sub>10</sub>(alates + 1))") +
    scale_color_scico(serify("", "N<sub>0</sub>", ""), palette = "glasgow") +
    facet_grid( ~ K)

