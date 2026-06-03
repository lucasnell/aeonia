
#'
#' Density plots for larger landscape simulations.
#'
#' Must run 09-large-densities.R on the cluster first.
#'



source("_scripts/00-preamble.R")


.overwrite <- FALSE


if (!file.exists(interm_files$dens_sims)) stop("Run 09-large-densities.R first!")


dens_sims <- read_csv(interm_files$dens_sims, col_types = "ciidddddddd") |>
    mutate(wasp_resp = factor(wasp_resp, levels = levels(wasp_resp_fct)))


total_dens_blank <- dens_sims |>
    select(wasp_resp, n_pseudo, time, aphids_sum, alates_sum, wasps_sum) |>
    rename_with(\(x) str_remove(x, "_sum")) |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    mutate(species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 50)
max_plant_dens_blank <- dens_sims |>
    select(time, aphids_max, alates_max, wasps_max) |>
    rename_with(\(x) str_remove(x, "_max")) |>
    group_by(time) |>
    summarize(across(aphids:wasps, max), .groups = "drop") |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    mutate(species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 50)


total_dens_p_list <- levels(wasp_resp_fct) |>
    set_names() |>
    map(\(wr) {
        dens_sims |>
            filter(wasp_resp == wr) |>
            select(n_pseudo, time, aphids_sum, alates_sum, wasps_sum) |>
            rename_with(\(x) str_remove(x, "_sum")) |>
            pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
            mutate(species = factor(species, levels = levels(total_dens_blank$species)),
                   n_pseudo = factor(n_pseudo)) |>
            ggplot(aes(time, density / 1e6)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_blank(data = total_dens_blank) +
            geom_line(aes(color = n_pseudo), linewidth = 0.75) +
            scale_color_manual(values = full_np_pal) +
            scale_y_continuous(n.breaks = 4L, expand = expansion(c(0.05, 0.2))) +
            facet_wrap(~ species, ncol = 1, scales = "free_y") +
            labs(x = "Time (days)", y = "Total density (&times; 10<sup>6</sup>)",
                 title = scenario_title(wr, TRUE))
    })





max_plant_dens_p_list <- levels(wasp_resp_fct) |>
    set_names() |>
    map(\(wr) {
        dens_sims |>
            filter(wasp_resp == wr) |>
            select(n_pseudo, time, aphids_max, alates_max, wasps_max) |>
            rename_with(\(x) str_remove(x, "_max")) |>
            pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
            mutate(species = factor(species, levels = levels(max_plant_dens_blank$species)),
                   n_pseudo = factor(n_pseudo)) |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_blank(data = max_plant_dens_blank) +
            geom_line(aes(color = n_pseudo), linewidth = 0.75) +
            scale_color_manual(values = full_np_pal) +
            facet_wrap(~ species, ncol = 1, scales = "free_y") +
            labs(x = "Time (days)", y = "Maximum per-plant density",
                 title = scenario_title(wr, TRUE))
    })



# wrap_plots(total_dens_p_list, nrow = 1, guides = "collect", axes = "collect")

# wrap_plots(max_plant_dens_p_list, nrow = 1, guides = "collect", axes = "collect")



total_dens_p <- total_dens_p_list[[1]] +
    (total_dens_p_list[[2]] + theme(axis.text.y = element_blank())) +
    plot_layout(design = "A#B", widths = c(1, 0.05, 1),
                axis_titles = "collect", guides = "collect") &
    illustrator_theme
# total_dens_p

if (.overwrite) {
    save_plot("_plots/densities-large.pdf", total_dens_p, width = 5, height = 2)
}


