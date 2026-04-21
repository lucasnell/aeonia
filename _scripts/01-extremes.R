

#'
#' Time series and histograms for simulations of two scenarios:
#' one where Pseudomonas decreases outbreaks ("low"),
#' and another where it increases outbreaks ("high")
#'

source("_scripts/00-preamble.R")

.overwrite <- FALSE

if (! dir.exists("_plots/extremes")) dir.create("_plots/extremes")




one_test <- function(Y0, N0, zeta) {

    # N0 = 55; Y0 = 1; zeta = c(1, 0.3)
    # rm(N0, Y0, zeta)

    stopifnot(length(N0) == 1 && is.numeric(N0) && N0 > 0)
    stopifnot(length(Y0) == 1 && is.numeric(Y0) && Y0 > 0)
    stopifnot(length(zeta) == 2 && is.numeric(zeta))
    stopifnot(all(zeta >= 0) && all(zeta <= 1))

    zeta <- sort(zeta, TRUE) |> set_names(c("low", "high"))

    crossing(np = c(3L, 0L), tp = c("low", "high")) |>
        pmap(\(np, tp) {
            args <- list(type = tp, large_sims = TRUE)
            args[["Y0"]] <- Y0
            args[["N0"]] <- N0
            args[["n_pseudo"]] <- np
            args[["zeta"]] <- zeta[[tp]]
            do.call(run_sim_combos, args) |>
                select(n_infected) |>
                mutate(p_emerge = as.integer(n_infected > 1),
                       outbreak_size = ifelse(n_infected > 1, n_infected, NA)) |>
                summarize(across(everything(), \(x) mean(x, na.rm = TRUE))) |>
                mutate(type = tp, n_pseudo = np)
        }) |>
        list_rbind() |>
        select(type, n_pseudo, n_infected, p_emerge, outbreak_size) |>
        mutate(type = factor(type, levels = c("low", "high"))) |>
        arrange(type, n_pseudo)
}



(sims <- one_test(Y0 = 1, N0 = 55, zeta = c(1, 0.3))); sims |>
    group_by(type) |>
    summarize(across(n_infected:outbreak_size,
                     \(x) x[n_pseudo == 3L] - x[n_pseudo == 0L]))













# ============================================================================*
# Larger sims ----
# ============================================================================*

#'
#' Do a bunch of simulations for each set of parameter values to show how these
#' parameter sets behave over many simulations.
#'



set.seed(259619623)
large_sims <- crossing(type = factor(1:2, labels = c("low", "high")),
                       n_pseudo = c(0L, 3L)) |>
    mutate(sims = map2(type, n_pseudo, \(type, n_pseudo) {
        run_sim_combos(type = type, n_pseudo = n_pseudo, large_sims = TRUE,
                       fly_p = 0.1)
    }))



large_sims |>
    mutate(n_infected = num(map_dbl(sims, \(x) mean(x$n_infected)), digits = 3)) |>
    select(-sims)
#   type  n_pseudo n_infected
#   <fct>    <int>  <num:.3!>
# 1 low          0      5.399
# 2 low          3      1.948
# 3 high         0      2.191
# 4 high         3      5.796


large_sims |>
    mutate(n_infected = map_dbl(sims, \(x) mean(x$n_infected))) |>
    group_by(type) |>
    summarize(n_infected = n_infected[n_pseudo != 0] -
                  n_infected[n_pseudo == 0]) |>
    mutate(n_infected = num(n_infected, digits = 3))
#   type  n_infected
#   <fct>  <num:.3!>
# 1 low       -3.451
# 2 high       3.605


# ============================================================================*
# Time series ----
# ============================================================================*


set.seed(1001450726)
low_high_sims <- large_sims |>
    select(-sims) |>
    mutate(sims = map2(type, n_pseudo, \(type, n_pseudo) {
        # type = "high"; n_pseudo = 3
        # rm(type, n_pseudo, sims, target, .rep)
        sims <- run_sim_combos(type, n_pseudo, n_sims = 100L, out_attack_surv = TRUE)
        target <- large_sims |>
            filter(type == .env$type, n_pseudo == .env$n_pseudo) |>
            getElement("sims") |>
            map_dbl(\(x) mean(x$n_infected))
        # Choose a representative simulation:
        .rep <- sims |>
            filter(is.na(x)) |>
            group_by(rep) |>
            summarize(n_infected = max(virus)) |>
            filter(abs(n_infected - target) ==
                       min(abs(n_infected - target))) |>
            getElement("rep")
        stopifnot(length(.rep) > 0)
        .rep <- sample(.rep, 1)
        sims |>
            filter(rep == .rep) |>
            mutate(plant = interaction(y, x),
                   aphids = aphids + parasitized) |>
            select(plant, time, virus, aphids, alates, wasps, attack_surv)
    })) |>
    unnest(sims)

low_high_sims |>
    filter(is.na(plant)) |>
    group_by(type, n_pseudo) |>
    summarize(n_infected = max(virus), .groups = "drop")
# # A tibble: 4 × 3
#   type  n_pseudo n_infected
#   <fct>    <int>      <dbl>
# 1 low          0          5
# 2 low          3          2
# 3 high         0          2
# 4 high         3          6







# -------------------------------------*
# ... total alates ----
# -------------------------------------*


tot_low_high_p_list <- low_high_sims |>
    split(~ type) |>
    map(\(x) {
        ym <- low_high_sims |>
            filter(is.na(plant)) |>
            getElement("alates") |>
            max()
        x |>
            filter(is.na(plant)) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            ggplot(aes(time, alates)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(aes(linetype = n_pseudo), linewidth = 1, color = "black") +
            scale_linetype_manual(values = np_linetypes) +
            coord_cartesian(ylim = c(0, ym))
    })

wrap_plots(tot_low_high_p_list, ncol = 1, guides = "collect", axis_titles = "collect")

if (.overwrite) {
    for (n in names(tot_low_high_p_list)) {
        save_plot(sprintf("_plots/extremes/timeseries-total-alates-%s.pdf", n),
                  tot_low_high_p_list[[n]] + illustrator_theme,
                  width = 2.5, height = 1.5)
    }; rm(n)
}





# -------------------------------------*
# ... by plant ----
# -------------------------------------*

# This is used to get the same axis max values across plots:
empty_data <- low_high_sims |>
    filter(!is.na(plant)) |>
    select(-virus, -attack_surv) |>
    pivot_longer(aphids:wasps, names_to = "species",
                 values_to = "density") |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 1,
           species = factor(species, levels = c("aphids", "alates", "wasps")))




by_plant_plotter <- function(type, .combine = FALSE) {
    # type == "high"
    # .combine = FALSE
    # rm(type, .combine, vdf, dd, breaks_fun, plant_p_list)

    vdf <- low_high_sims |>
        filter(type == .env$type) |>
        filter(!is.na(plant), virus == 1) |>
        group_by(n_pseudo, plant) |>
        summarize(time = min(time), .groups = "drop") |>
        mutate(density = max(empty_data$density) * 1.1,
               n_pseudo = factor(n_pseudo),
               species = sort(empty_data$species)[1])

    dd <- low_high_sims |>
        filter(type == .env$type) |>
        filter(!is.na(plant)) |>
        select(n_pseudo, plant, time, aphids, alates, wasps)

    breaks_fun <- function(x) {
        if (max(x) >= 1200) {
            return(c(0, 600, 1200))
        } else if (max(x) >= 120) {
            return(c(0, 60, 120))
        } else return(c(0, 30, 60))
    }

    plant_p_list <- levels(dd$plant) |>
        set_names() |>
        map(\(plant) {
            p <- dd |>
                filter(plant == .env$plant) |>
                pivot_longer(aphids:wasps, names_to = "species",
                             values_to = "density") |>
                mutate(n_pseudo = factor(n_pseudo),
                       species = factor(species, levels =
                                            levels(empty_data$species))) |>
                ggplot(aes(time, density)) +
                geom_hline(yintercept = 0, color = "gray70") +
                geom_line(aes(color = species, linetype = n_pseudo,
                              linewidth = n_pseudo)) +
                geom_blank(data = empty_data, aes(y = density * 1.15))
            if (plant %in% vdf$plant) {
                p <- p +
                    geom_point(aes(shape = n_pseudo),
                               data = filter(vdf, plant == .env$plant),
                               color = color_pal[["virus"]],
                               size = 2, stroke = 1) +
                    scale_shape_manual(values = np_shapes)
            }
            p +
                scale_color_manual(values = color_pal) +
                scale_linetype_manual(values = np_linetypes) +
                scale_linewidth_manual(values = np_linewidths) +
                labs(x = "Time (days)", y = "Density", tag = plant) +
                scale_x_continuous(breaks = c(0, 50, 100)) +
                scale_y_continuous(breaks = breaks_fun) +
                coord_cartesian(clip = "off") +
                facet_wrap(~ species, ncol = 1, scales = "free_y") +
                guides(linetype = guide_legend(
                    override.aes = list(color = "black"))) +
                theme(axis.text.x = element_markdown(size = 7),
                      axis.text.y = element_markdown(size = 7),
                      strip.text.x = element_blank(),
                      plot.margin = margin(3,3,0,0))
        })

    if (!.combine) return(plant_p_list)

    wrap_plots(plant_p_list, ncol = 3, guides = "collect",
               axis_titles = "collect")

}

# by_plant_plotter("low", .combine = TRUE)
# by_plant_plotter("high", .combine = TRUE)


if (.overwrite) {
    for (ty in levels(low_high_sims$type)) {
        plot_list <- by_plant_plotter(ty)
        for (pl in names(plot_list)) {
            fn <- sprintf("_plots/extremes/timeseries-%s-%s.pdf", ty, pl)
            p <- plot_list[[pl]] &
                illustrator_theme &
                theme(axis.text.y = element_blank(),
                      axis.text.x = element_blank())
            save_plot(fn, p, width = 1.5, height = 1.5)
        }
    }; rm(ty, plot_list, pl, fn, p)
}






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
            geom_point(data = dd, aes(shape = n_pseudo), size = 4, stroke = 1) +
            scale_shape_manual(values = np_shapes) +
            scale_x_continuous(breaks = c(0, 500, 1000, 1500))
    })

wrap_plots(z_pa_p_list, ncol = 1, guides = "collect", axis_titles = "collect")

if (.overwrite) {
    for (n in names(z_pa_p_list)) {
        save_plot(sprintf("_plots/extremes/z-p_alates-%s.pdf", n),
                  z_pa_p_list[[n]] & illustrator_theme, width = 1.5, height = 1.5)
    }; rm(n)
}




# Curve for model diagram:
z_pa_curve <- tibble(z = 1:3e3, ap = alate_prop(z)) |>
    ggplot(aes(z, ap)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(linewidth = 1) +
    scale_x_continuous(breaks = 0:3*1000)

if (.overwrite) {
    save_plot("_plots/z-p_alates.pdf", z_pa_curve + illustrator_theme,
              width = 2, height = 2)
}






# ============================================================================*
# Outbreak sizes - freq. polygons ----
# ============================================================================*



low_high_hist_list <- levels(large_sims$type) |>
    set_names() |>
    map(\(tp) {
        # tp = "low"
        # rm(tp, d, dd, dd_means)
        d  <- large_sims |>
            filter(type == tp)|>
            mutate(n_pseudo = factor(n_pseudo)) |>
            mutate(outbreak_size = map(sims, \(x) x$n_infected)) |>
            select(-sims) |>
            unnest(outbreak_size) |>
            filter(outbreak_size > 1)
        dd <- d |>
            mutate(outbreak_size = factor(outbreak_size, levels = 2:9)) |>
            group_by(n_pseudo, outbreak_size) |>
            count(name = "n_obs", .drop = FALSE) |>
            group_by(n_pseudo) |>
            mutate(perc = 100 * n_obs / sum(n_obs)) |>
            ungroup() |>
            mutate(outbreak_size = as.integer(paste(outbreak_size)))
        dd_means <- d |>
            group_by(n_pseudo) |>
            summarize(outbreak_size = mean(outbreak_size),
                      .groups = "drop")
        dd |>
            ggplot(aes(outbreak_size, perc)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_point(aes(shape = n_pseudo), size = 2, stroke = 1) +
            geom_line(aes(linetype = n_pseudo), linewidth = 1) +
            # geom_pointrange(data = dd_means, aes(xmin = lo, xmax = hi),
            #                 linewidth = 1.5, size = 0.75, shape = 1) +
            geom_vline(data = dd_means,
                       aes(xintercept = outbreak_size, linetype = n_pseudo),
                       linewidth = 1) +
            scale_shape_manual(values = np_shapes) +
            scale_linetype_manual(values = np_linetypes) +
            scale_x_continuous(breaks = (0:4) * 2 + 1) +
            labs(x = "Outbreak size", y = "Percent of simulations")
    })

wrap_plots(low_high_hist_list, ncol = 1, guides = "collect", axis_titles = "collect")

if (.overwrite) {
    for (n in names(low_high_hist_list)) {
        save_plot(sprintf("_plots/extremes/histograms-%s.pdf", n),
                  low_high_hist_list[[n]] + illustrator_theme,
                  width = 2.5, height = 1.5)
    }; rm(n)
}





# ============================================================================*
# Pr(emergence) bar graphs ----
# ============================================================================*


low_high_bar_list <- levels(large_sims$type) |>
    set_names() |>
    map(\(tp) {
        large_sims |>
            filter(type == tp)|>
            mutate(n_pseudo = factor(n_pseudo)) |>
            mutate(p_emerge = map_dbl(sims, \(x) mean(x$n_infected > 1))) |>
            ggplot(aes(p_emerge, n_pseudo)) +
            geom_vline(xintercept = 0, color = "gray70") +
            geom_col(aes(fill = n_pseudo), color = "black", width = 0.45,
                     linewidth = 0.75, linejoin = "mitre") +
            labs(x = "Prob. of emergence", y = "Number of Pseudo. patches") +
            coord_cartesian(xlim  = c(0, 1)) +
            scale_fill_manual(values = c(`0` = "white", `3` = "black")) +
            theme(axis.text.y = element_markdown(color = NA))
    })

# wrap_plots(low_high_bar_list, ncol = 1, guides = "collect", axis_titles = "collect")

if (.overwrite) {
    for (n in names(low_high_bar_list)) {
        save_plot(sprintf("_plots/extremes/barplots-%s.pdf", n),
                  low_high_bar_list[[n]] + illustrator_theme,
                  width = 1.75, height = 1.5)
    }; rm(n)
}





# ============================================================================*
# Empirical zeta ----
# ============================================================================*


# In Ives et al. (1999), they found that parasitoids spent ~3.76 more time
# foraging at plants where they encountered an aphid.
# Below simulates our model with varying zeta, then compares the observed
# wasp abundances to those predicted when parasitoids never encounter an aphid
# on a Pseudomonas-inhabited plant but always do on plants without Pseudomonas.

set.seed(1180209329)
zeta_low_high_sims <- map(c(0.5, 0.6, 0.7, 0.8, 0.9), \(zeta) {
        sims <- run_sim_combos(type = "high", n_pseudo = 3L, n_sims = 12L,
                               out_attack_surv = TRUE, zeta = zeta) |>
            filter(!is.na(x)) |>
            mutate(plant = interaction(y, x),
                   aphids = aphids + parasitized) |>
            mutate(zeta = .env$zeta) |>
            select(zeta, rep, plant, time, virus, aphids, alates, wasps,
                   attack_surv) |>
            split(~ rep + time) |>
            map(\(x) {
                Y <- sum(x$wasps)
                p_plants <- x$plant %in% c("3.1", "2.2", "1.3")
                pred_wasps <- numeric(nrow(x))
                pred_wasps[p_plants] <- 1
                pred_wasps[!p_plants] <- 3.76
                pred_wasps <- Y * pred_wasps / sum(pred_wasps)
                x[["pred_wasps"]] <- pred_wasps
                return(x)
            }) |>
            list_rbind()
    }) |>
    list_rbind()



empir_zeta_p <- zeta_low_high_sims |>
    split(~ zeta + rep) |>
    map(\(x) {
        max_N_t <- x |>
            group_by(time) |>
            summarize(aphids = max(aphids + alates)) |>
            filter(aphids == max(aphids)) |>
            getElement("time") |> getElement(1)
        return(x |> filter(time == max_N_t))
    }) |>
    list_rbind() |>
    select(zeta, rep, plant, wasps, pred_wasps) |>
    mutate(pseudo = factor(plant %in% c("3.1", "2.2", "1.3"),
                           levels = c(TRUE, FALSE),
                           labels = c("*Pseudomonas*", "No *Pseudomonas*")),
           zeta = factor(zeta),
           rep = factor(rep),
           rel = wasps / pred_wasps) |>
    ggplot(aes(zeta, rel)) +
    geom_jitter(aes(color = zeta), height = 0, width = 0.2, shape = 1, size = 2) +
    geom_hline(yintercept = 1, color = "black", linewidth = 1.5, linetype = "22") +
    # facet_wrap(~ plant, nrow = 3, scales = "free") +
    facet_wrap(~ pseudo, nrow = 1, scales = "free") +
    labs(x = pretty_params("zeta", cap1 = TRUE),
         y = "Observed / predicted parasitoid density") +
    scale_color_viridis_d(guide = "none")



# empir_zeta_p

if (.overwrite) {
    save_plot("_plots/empirical-zeta.pdf", empir_zeta_p, width = 6, height = 3)
}



