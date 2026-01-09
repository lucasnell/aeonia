
#'
#' Time series and histograms for simulations of two scenarios:
#' one where Pseudomonas decreases outbreaks ("low"),
#' and another where it increases outbreaks ("high")
#'

source("_scripts/00-preamble.R")



# ============================================================================*
# Larger sims ----
# ============================================================================*

#'
#' Do a bunch of simulations for each set of parameter values to show how these
#' parameter sets behave over many simulations.
#'

if (!file.exists(rds_files$extreme_large)) {

    # Takes very little time, but I'm saving RDS to use output in other scripts.
    set.seed(259619622)
    large_sims <- crossing(type = factor(1:2, labels = c("low", "high")),
                           n_pseudo = c(0L, 3L)) |>
        mutate(outbreak_size = map2(type, n_pseudo, \(type, n_pseudo) {
            run_sim_combos(type, n_pseudo, TRUE) |>
                getElement("outbreak_size")
        }))

    write_rds(large_sims, rds_files$extreme_large, compress = "gz")

} else {

    large_sims <- read_rds(rds_files$extreme_large)

}

large_sims |>
    mutate(outbreak_size = num(map_dbl(outbreak_size, mean), digits = 3))
# # A tibble: 4 × 3
#   type  n_pseudo outbreak_size
#   <fct>    <int>     <num:.3!>
# 1 low          0         6.938
# 2 low          3         4.946
# 3 high         0         2.970
# 4 high         3         4.966

large_sims |>
    mutate(outbreak_size = map_dbl(outbreak_size, mean)) |>
    group_by(type) |>
    summarize(outbreak_size = outbreak_size[n_pseudo != 0] -
                  outbreak_size[n_pseudo == 0]) |>
    mutate(outbreak_size = num(outbreak_size, digits = 3))
# # A tibble: 2 × 2
#   type  outbreak_size
#   <fct>     <num:.3!>
# 1 low          -1.992
# 2 high          1.996






# ============================================================================*
# Time series ----
# ============================================================================*


set.seed(1001450726)
low_high_sims <- large_sims |>
    select(-outbreak_size) |>
    mutate(sims = map2(type, n_pseudo, \(type, n_pseudo) {
        # type = "low"; n_pseudo = 3
        # rm(type, n_pseudo, sims, target, .rep)
        sims <- run_sim_combos(type, n_pseudo, n_sims = 100L)
        target <- large_sims |>
            filter(type == .env$type, n_pseudo == .env$n_pseudo) |>
            getElement("outbreak_size") |>
            map_dbl(mean)
        # Choose a representative simulation:
        .rep <- sims |>
            filter(is.na(x)) |>
            group_by(rep) |>
            summarize(outbreak_size = max(virus)) |>
            filter(abs(outbreak_size - target) == min(abs(outbreak_size - target))) |>
            getElement("rep")
        stopifnot(length(.rep) > 0)
        sims |>
            filter(rep == .rep[1]) |>
            mutate(plant = interaction(y, x),
                   aphids = aphids + parasitized) |>
            select(plant, time, virus, aphids, alates, wasps)
    })) |>
    unnest(sims)






# -------------------------------------*
# ... total alates ----
# -------------------------------------*

# This is used to get the same axis max values across plots:
tot_empty_data <- tibble(time = 0,
                         alates = low_high_sims |>
                             filter(is.na(plant)) |>
                             getElement("alates") |>
                             max())


tot_low_high_p_list <- low_high_sims |>
    split(~ type) |>
    map(\(x) {
        x |>
            filter(is.na(plant)) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            ggplot(aes(time, alates)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(aes(color = n_pseudo), linewidth = 1) +
            geom_blank(data = tot_empty_data) +
            scale_color_manual(values = np_pal)
    })

wrap_plots(tot_low_high_p_list, ncol = 1, guides = "collect", axis_titles = "collect")


# for (n in names(tot_low_high_p_list)) {
#     save_plot(sprintf("_plots/extremes-timeseries-total-alates-%s.pdf", n),
#               tot_low_high_p_list[[n]] + illustrator_theme, width = 2.5, height = 1.5)
# }; rm(n)





# -------------------------------------*
# ... by plant ----
# -------------------------------------*

# This is used to get the same axis max values across plots:
empty_data <- low_high_sims |>
    filter(!is.na(plant)) |>
    select(-virus) |>
    pivot_longer(aphids:wasps, names_to = "species",
                 values_to = "density") |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 0,
           species = factor(species, levels = c("aphids", "wasps", "alates")),
           density = ifelse(species == "aphids", density,
                            max(density[species != "aphids"])))

# Object `dens_maxes` depends on the following `empty_data` values:
stopifnot(max(empty_data$density) > 1000 && max(empty_data$density) < 1100)
stopifnot(min(empty_data$density) > 44 && min(empty_data$density) < 45)
dens_maxes <- list(aphids = 1000, alates = 45, wasps = 45)


low_high_p_list <- low_high_sims |>
    mutate(grp = interaction(type, n_pseudo, sep = "-np") |>
               factor(levels = paste(rep(c("low", "high"), each = 2),
                                     rep(c("3", "0"), 2), sep = "-np"))) |>
    split(~ grp) |>
    map(\(x) {
        # x = low_high_sims |> filter(type == "low", n_pseudo == 3)
        # rm(x, vdf, trans, itrans, dd, p)
        vdf <- x |>
            filter(!is.na(plant), virus == 1) |>
            group_by(plant) |>
            summarize(time = min(time)) |>
            mutate(density = max(empty_data$density) * 1.1)
        # wasps / alates --> aphids
        trans <- function(x) {
            x * dens_maxes[["aphids"]] / dens_maxes[["wasps"]]
        }
        # aphids --> wasps / alates
        itrans <- function(x) {
            x * dens_maxes[["wasps"]] / dens_maxes[["aphids"]]
        }
        dd <- x |>
            filter(!is.na(plant)) |>
            select(-virus) |>
            mutate(wasps = trans(wasps),
                   alates = trans(alates))

        p <- dd |>
            # select(-alates) |>
            pivot_longer(aphids:wasps, names_to = "species",
                         values_to = "density") |>
            mutate(species = factor(species, levels = levels(empty_data$species))) |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(aes(color = species, alpha = species),  # , linetype = species),
                      linewidth = 1) +
            geom_blank(data = empty_data, aes(y = density * 1.15)) +
            geom_point(data = vdf, color = color_pal[["virus"]],
                       shape = 8, size = 2) +
            scale_color_manual(values = color_pal) +
            scale_alpha_manual(values = c(1, 1, 1) |>
                                   set_names(c("aphids", "alates", "wasps"))) +
            # scale_linetype_manual(values = c(aphids = "solid",
            #                                  alates = "22",
            #                                  wasps = "33")) +
            labs(x = "Time (days)", y = "Density",
                 title = sprintf("n<sub>P</sub> = %i, type = %s",
                                 x$n_pseudo[[1]], x$type[[1]])) +
            scale_y_continuous(sec.axis = sec_axis(itrans,
                                                   "Density (wasps/alates)",
                                                   breaks = c(0, 20, 40)),
                               breaks = c(0, 500, 1000)) +
            scale_x_continuous(breaks = c(0, 50, 100)) +
            facet_wrap(~ plant) +
            theme(axis.text.x = element_markdown(size = 7),
                  axis.text.y = element_markdown(size = 7),
                  plot.margin = margin(3,3,0,0))
        if (x$n_pseudo[[1]] == 0) {
            p <- p + theme(axis.text.y.left = element_blank(),
                           # axis.ticks.y.left = element_blank(),
                           axis.title.y.left = element_blank())
        } else {
            p <- p + theme(axis.text.y.right = element_blank(),
                           # axis.ticks.y.right = element_blank(),
                           axis.title.y.right = element_blank())
        }
        return(p)
    }) |>
    # Now group by low vs high:
    (\(x) {
        np_vals <- rev(sort(unique(low_high_sims$n_pseudo)))
        out <- levels(low_high_sims$type) |>
            set_names() |>
            map(\(tp) {
                map(np_vals, \(np) x[[paste0(tp, "-np", np)]]) |>
                    wrap_plots(nrow = 1, guides = "collect", design = "A#B",
                               widths = c(1, 0.05, 1))
            })
        return(out)
    })()



wrap_plots(low_high_p_list, nrow = 2, guides = "collect", axis_titles = "collect")


# for (n in names(low_high_p_list)) {
#     save_plot(sprintf("_plots/extremes-timeseries-%s.pdf", n),
#               low_high_p_list[[n]] & illustrator_theme,
#               width = 5, height = 2.5)
# }; rm(n)







# ============================================================================*
# z --> Pr(alates) ----
# ============================================================================*

z_pa_p_list <- low_high_sims |>
    split(~ type) |>
    map(\(x) {
        # x = low_high_sims |> filter(type == "low")
        # rm(x, dd, max_z)
        dd <- x |>
            filter(!is.na(plant)) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            filter(! (plant %in% c("3.1", "2.2", "1.3") & n_pseudo != 0)) |>
            group_by(n_pseudo, plant) |>
            summarize(z = max(aphids + alates), .groups = "drop_last") |>
            summarize(z = max(z)) |>
            mutate(ap = alate_prop(z))
        max_z <- ceiling((max(empty_data$density) * 1.1) / 100) * 100
        tibble(z = 1:max_z, ap = alate_prop(z)) |>
            ggplot(aes(z, ap)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(linewidth = 1) +
            geom_point(data = dd, aes(color = n_pseudo), size = 3) +
            scale_color_manual(values = np_pal) +
            scale_x_continuous(breaks = c(0, 500, 1000))
    })

wrap_plots(z_pa_p_list, ncol = 1, guides = "collect", axis_titles = "collect")


# for (n in names(low_high_p_list)) {
#     save_plot(sprintf("_plots/extremes-z-p_alates-%s.pdf", n),
#               z_pa_p_list[[n]] & illustrator_theme, width = 1.5, height = 1.5)
# }; rm(n)



# Curve for model diagram:
z_pa_curve <- tibble(z = 1:3e3, ap = alate_prop(z)) |>
    ggplot(aes(z, ap)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = 0:3*1000)

# save_plot("_plots/z-p_alates.pdf", z_pa_curve + illustrator_theme,
#           width = 2, height = 2)





# ============================================================================*
# Histograms ----
# ============================================================================*



low_high_hist_list <- levels(large_sims$type) |>
    set_names() |>
    map(\(tp) {
        # tp = "low"
        # rm(tp, d, dd, dd_means)
        d  <- large_sims |>
            filter(type == tp)
        dd <- d |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            unnest(outbreak_size) |>
            mutate(outbreak_size = factor(outbreak_size, levels = 1:9)) |>
            group_by(n_pseudo, outbreak_size) |>
            count(name = "perc", .drop = FALSE) |>
            ungroup() |>
            mutate(outbreak_size = as.integer(paste(outbreak_size)),
                   perc = 100 * perc / 1000)
        dd_means <- d |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            unnest(outbreak_size) |>
            group_by(n_pseudo) |>
            summarize(# boot = list(boot_ci(outbreak_size)),
                      outbreak_size = mean(outbreak_size),
                      .groups = "drop") |>
            # unnest(boot) |>
            mutate(perc = max(dd$perc) * 0.95)
        dd |>
            ggplot(aes(outbreak_size, perc, color = n_pseudo)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_point(size = 2) +
            geom_line(linewidth = 1) +
            # geom_pointrange(data = dd_means, aes(xmin = lo, xmax = hi),
            #                 linewidth = 1.5, size = 0.75, shape = 1) +
            geom_point(data = dd_means, size = 4) +
            scale_color_manual(values = np_pal, aesthetics = c("color", "fill")) +
            scale_x_continuous(breaks = (0:4) * 2 + 1) +
            labs(x = "Outbreak size", y = "Percent of simulations")
    })

wrap_plots(low_high_hist_list, ncol = 1, guides = "collect", axis_titles = "collect")


# for (n in names(low_high_hist_list)) {
#     save_plot(sprintf("_plots/extremes-histograms-%s.pdf", n),
#               low_high_hist_list[[n]] + illustrator_theme, width = 4.5, height = 1.5)
# }; rm(n)

