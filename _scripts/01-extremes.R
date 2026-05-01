

#'
#' Time series and histograms for simulations of two levels of wasp
#' responsiveness to aphid densities:
#' "strong" (so Pseudomonas decreases outbreak size)
#' "weak" (so Pseudomonas increases outbreak size)
#'

source("_scripts/00-preamble.R")

.overwrite <- FALSE

if (! dir.exists("_plots/extremes")) dir.create("_plots/extremes")




# ============================================================================*
# Larger sims ----
# ============================================================================*

#'
#' Do a bunch of simulations for each set of parameter values to show how these
#' parameter sets behave over many simulations.
#'



set.seed(259619623)
large_sims <- crossing(wasp_resp = factor(1:2, labels = c("weak", "strong")),
                       n_pseudo = c(0L, 3L)) |>
    mutate(sims = map2(wasp_resp, n_pseudo, \(wasp_resp, n_pseudo) {
        run_sim_combos(wasp_resp = wasp_resp, n_pseudo = n_pseudo, large_sims = TRUE)
    }))



large_sims |>
    mutate(n_infected = num(map_dbl(sims, \(x) mean(x$n_infected)), digits = 3)) |>
    select(-sims)
#   wasp_resp n_pseudo n_infected
#   <fct>        <int>  <num:.3!>
# 1 weak             0      3.686
# 2 weak             3      6.491
# 3 strong           0      3.671
# 4 strong           3      2.141


large_sims |>
    mutate(n_infected = map_dbl(sims, \(x) mean(x$n_infected))) |>
    group_by(wasp_resp) |>
    summarize(n_infected = n_infected[n_pseudo != 0] -
                  n_infected[n_pseudo == 0]) |>
    mutate(n_infected = num(n_infected, digits = 3))
#   wasp_resp n_infected
#   <fct>      <num:.3!>
# 1 weak           2.805
# 2 strong        -1.530


# ============================================================================*
# Time series ----
# ============================================================================*


set.seed(1001450726)
weak_strong_sims <- large_sims |>
    select(-sims) |>
    mutate(sims = map2(wasp_resp, n_pseudo, \(wasp_resp, n_pseudo) {
        # wasp_resp = "weak"; n_pseudo = 3
        # rm(wasp_resp, n_pseudo, sims, target, .rep, half_inf_times, hit_adiffs)
        sims <- run_sim_combos(wasp_resp, n_pseudo, n_sims = 1000L,
                               out_attack_surv = TRUE, out_stages = TRUE)
        # target <- large_sims |>
        #     filter(wasp_resp == .env$wasp_resp, n_pseudo == .env$n_pseudo) |>
        #     getElement("sims") |>
        #     map_dbl(\(x) mean(x$n_infected))
        # # Choose a representative simulation:
        # .rep <- sims |>
        #     filter(is.na(x)) |>
        #     group_by(rep) |>
        #     summarize(n_infected = max(virus)) |>
        #     filter(abs(n_infected - target) ==
        #                min(abs(n_infected - target))) |>
        #     getElement("rep")
        # stopifnot(length(.rep) > 0)
        # sims <- sims |>
        #     filter(rep %in% .rep)
        # # Now find rep with median time to halfway max infected:
        # if (length(.rep) > 1) {
        #     half_inf_times <- sims |>
        #         split(~ rep) |>
        #         map(\(x) {
        #             max_inf <- max(x$virus)
        #             tibble(rep = x$rep[[1]],
        #                    time = min(x$time[is.na(x$x) & x$virus >= (max_inf / 2)]))
        #         }) |>
        #         list_rbind()
        #     hit_adiffs <- abs(half_inf_times$time - median(half_inf_times$time))
        #     .rep <- half_inf_times$rep[hit_adiffs == min(hit_adiffs)]
        #     if (length(.rep) > 1) .rep <- sample(.rep, 1)
        #     sims <- sims |>
        #         filter(rep %in% .rep)
        # }
        sims |>
            # filter(rep == .rep) |>
            mutate(plant = interaction(y, x),
                   aphids = aphids_juv + aphids_adu + alates_juv + alates_adu + parasitized,
                   alates = alates_adu) |>
            select(rep, plant, time, virus, starts_with("aphids"),
                   starts_with("alates"), wasps, attack_surv) |>
            group_by(time, plant) |>
            summarize(across(virus:attack_surv, mean), .groups = "drop")
    })) |>
    unnest(sims)

weak_strong_sims |>
    filter(is.na(plant)) |>
    group_by(wasp_resp, n_pseudo) |>
    summarize(n_infected = max(virus), .groups = "drop")
#   wasp_resp n_pseudo n_infected
#   <fct>        <int>      <dbl>
# 1 weak             0       3.69
# 2 weak             3       6.50
# 3 strong           0       3.66
# 4 strong           3       2.04







# -------------------------------------*
# ... totals ----
# -------------------------------------*

# This is used to get the same axis max values across plots:
tot_empty_data <- weak_strong_sims |>
    filter(is.na(plant)) |>
    select(virus, aphids, alates, wasps) |>
    mutate(aphids = aphids / 100) |>
    pivot_longer(virus:wasps, names_to = "species", values_to = "density") |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 1,
           species = factor(species, levels = c("aphids", "alates",
                                                "wasps", "virus")),
           density = ifelse(species == "virus", 9, density)) |>
    arrange(desc(density))


total_plotter <- function(wasp_resp, delay_virus = FALSE) {
    # wasp_resp = "weak"; delay_virus = FALSE
    # rm(wasp_resp, delay_virus, dd)

    dd <- weak_strong_sims |>
        filter(wasp_resp == .env$wasp_resp) |>
        filter(is.na(plant)) |>
        select(n_pseudo, time, virus, aphids, alates, wasps) |>
        mutate(aphids = aphids / 100) |>
        pivot_longer(virus:wasps, names_to = "species",
                     values_to = "density") |>
        mutate(species = factor(species, levels =
                                    levels(tot_empty_data$species)))
    if (delay_virus) {
        dd <- dd |>
            mutate(time = ifelse(species == "virus", time - 7L, time)) |>
            filter(time >= 0L)
    }

    dd |>
        # filter(n_pseudo > 0) |>
        mutate(n_pseudo = factor(n_pseudo)) |>
        ggplot(aes(time, density)) +
        geom_hline(yintercept = 0, color = "gray70") +
        # geom_area(data = np_df, aes(fill = species), alpha = 0.4) +
        # geom_line(aes(color = species), linewidth = 0.75) +
        # geom_line(aes(color = species, linewidth = n_pseudo, linetype = n_pseudo)) +
        geom_line(aes(color = n_pseudo), linewidth = 1.25) +
        geom_blank(data = tot_empty_data, aes(y = density * 1.15)) +
        # scale_color_manual(values = color_pal, aesthetics = c("color", "fill")) +
        scale_color_manual(values = np_pal, aesthetics = c("color", "fill")) +
        # scale_linetype_manual(values = np_linetypes) +
        # scale_linewidth_manual(values = np_linewidths) +
        labs(x = "Time (days)", y = "Total density") +
        scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
        scale_y_continuous(breaks = \(x) {
            if (max(x) >= 400) {# wasps
                return(c(0, 200, 400))
            } else if (max(x) >= 120) {# aphids (total)
                return(c(0, 60, 120))
            } else if (max(x) >= 50) {# alates (adult)
                return(c(0, 25, 50))
            } else return(c(1, 5, 9))# virus
        }) +
        coord_cartesian(clip = "off") +
        facet_wrap(~ species, ncol = 1, scales = "free_y") +
        theme(strip.text.x = element_blank())

}


# total_plotter("weak")
# total_plotter("strong")







# -------------------------------------*
# ... by plant ----
# -------------------------------------*

# This is used to get the same axis max values across plots:
bp_empty_data <- weak_strong_sims |>
    filter(!is.na(plant)) |>
    select(aphids, alates, wasps) |>
    mutate(aphids = aphids / 100) |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 1,
           species = factor(species, levels = c("aphids", "alates", "wasps"))) |>
    arrange(desc(density))


by_plant_plotter <- function(wasp_resp) {
    # wasp_resp = "weak"
    # rm(wasp_resp, dd)

    dd <- weak_strong_sims |>
        filter(wasp_resp == .env$wasp_resp) |>
        filter(!is.na(plant)) |>
        select(n_pseudo, plant, time, aphids, alates, wasps) |>
        mutate(aphids = aphids / 100) |>
        pivot_longer(aphids:wasps, names_to = "species",
                     values_to = "density") |>
        mutate(species = factor(species, levels =
                                    levels(bp_empty_data$species))) |>
        mutate(plant_type = case_when(n_pseudo == 0L ~ "noP",
                                      plant %in% c("3.1", "2.2", "1.3") ~ "allP",
                                      TRUE ~ "PnoP")) |>
        group_by(n_pseudo, plant_type, time, species) |>
        summarize(density = mean(density), .groups = "drop") |>
        mutate(n_pseudo = factor(n_pseudo)) |>
        mutate(grp = interaction(plant_type, species, drop = TRUE))

    dd |>
        ggplot(aes(time, density)) +
        geom_hline(yintercept = 0, color = "gray70", linewidth = 0.5) +
        geom_line(aes(color = n_pseudo, linetype = plant_type, group = grp),
                  linewidth = 1.25) +
        geom_blank(data = bp_empty_data, aes(y = density * 1.15)) +
        scale_color_manual(values = np_pal, aesthetics = c("color", "fill")) +
        scale_linetype_manual(values = c("solid", "solid", "22") |>
                                   set_names(c("noP", "PnoP", "allP"))) +
        labs(x = "Time (days)", y = "Density per plant") +
        scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
        scale_y_continuous(breaks = \(x) {
            if (max(x) >= 50) {# wasps
                return(c(0, 25, 50))
            } else if (max(x) >= 14) {# aphids (total)
                return(c(0, 7, 14))
            } else return(c(0, 4, 8))# alates (adult)
        }) +
        coord_cartesian(clip = "off") +
        facet_wrap(~ species, ncol = 1, scales = "free_y") +
        theme(strip.text.x = element_blank())

}





# by_plant_plotter("weak")
# by_plant_plotter("strong")



# -------------------------------------*
# ... stitch together ----
# -------------------------------------*


timeseries_p <- crossing(wasp_resp = sort(unique(weak_strong_sims$wasp_resp)),
                         fxn = c(total_plotter, by_plant_plotter)) |>
    pmap(\(wasp_resp, fxn) {
        p <- fxn(wasp_resp) +
            theme(axis.text.x = element_markdown(size = 10),
                  axis.text.y = element_markdown(size = 10))
        if (wasp_resp == "weak") p <- p + theme(axis.text.y = element_blank())
        if (identical(fxn, total_plotter)) p <- p +
                theme(axis.text.x = element_blank())
        return(p)
    }) |>
    wrap_plots(design = "A#C\n###\nB#D", widths = c(1, 0.05, 1),
               heights = c(4, 0.35, 3), guides = "collect", axis_titles = "collect")
# timeseries_p


if (.overwrite) {
    save_plot("_plots/extremes/timeseries.pdf",
              timeseries_p & illustrator_theme, width = 6.9, height = 6)
}






# ============================================================================*
# z --> Pr(alates) ----
# ============================================================================*

z_pa_p_list <- weak_strong_sims |>
    split(~ wasp_resp) |>
    map(\(x) {
        # x = weak_strong_sims |> filter(wasp_resp == "strong")
        # rm(x, dd, max_z)
        dd <- x |>
            filter(!is.na(plant)) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            filter(! (plant %in% c("3.1", "2.2", "1.3") & n_pseudo != 0)) |>
            group_by(n_pseudo, plant) |>
            summarize(z = max(aphids + alates), .groups = "drop_last") |>
            summarize(z = max(z)) |>
            mutate(ap = alate_prop(z))
        max_z <- ceiling((100 * bp_empty_data$density[bp_empty_data$species == "aphids"] * 1.15) / 100) * 100
        tibble(z = 1:max_z, ap = alate_prop(z)) |>
            ggplot(aes(z, ap)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_line(linewidth = 1) +
            geom_point(data = dd, aes(color = n_pseudo), size = 4, stroke = 1) +
            scale_color_manual(values = np_pal) +
            scale_x_continuous(breaks = c(0, 500, 1000, 1500)) +
            theme(axis.text.x = element_markdown(size = 8),
                  axis.text.y = element_markdown(size = 8))
    })

# wrap_plots(z_pa_p_list, ncol = 1, guides = "collect", axis_titles = "collect")

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



weak_strong_hist_list <- levels(large_sims$wasp_resp) |>
    set_names() |>
    map(\(tp) {
        # tp = "strong"
        # rm(tp, d, dd, dd_means)
        d  <- large_sims |>
            filter(wasp_resp == tp)|>
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
            geom_point(aes(color = n_pseudo), size = 2, stroke = 1) +
            geom_line(aes(color = n_pseudo), linewidth = 1) +
            # geom_pointrange(data = dd_means, aes(xmin = lo, xmax = hi),
            #                 linewidth = 1.5, size = 0.75, shape = 1) +
            geom_vline(data = dd_means,
                       aes(xintercept = outbreak_size, color = n_pseudo),
                       linewidth = 1, linetype = "22") +
            scale_color_manual(values = np_pal) +
            # scale_linetype_manual(values = np_linetypes) +
            scale_x_continuous(breaks = (0:4) * 2 + 1) +
            labs(x = "Outbreak size", y = "Percent of simulations")
    })

# wrap_plots(weak_strong_hist_list, ncol = 1, guides = "collect", axis_titles = "collect")

if (.overwrite) {
    for (n in names(weak_strong_hist_list)) {
        save_plot(sprintf("_plots/extremes/histograms-%s.pdf", n),
                  weak_strong_hist_list[[n]] + illustrator_theme,
                  width = 3.2, height = 1.5)
    }; rm(n)
}





# ============================================================================*
# Pr(emergence) bar graphs ----
# ============================================================================*


weak_strong_bar_list <- levels(large_sims$wasp_resp) |>
    set_names() |>
    map(\(tp) {
        large_sims |>
            filter(wasp_resp == tp)|>
            mutate(n_pseudo = factor(n_pseudo)) |>
            mutate(p_emerge = map_dbl(sims, \(x) mean(x$n_infected > 1))) |>
            ggplot(aes(p_emerge, n_pseudo)) +
            geom_vline(xintercept = 0, color = "gray70") +
            geom_col(aes(fill = n_pseudo), color = NA, width = 0.45,
                     linewidth = 0.75, linejoin = "mitre") +
            labs(x = "Prob. of emergence", y = "Number of Pseudo. plants") +
            coord_cartesian(xlim  = c(0, 1)) +
            scale_fill_manual(values = np_pal) +
            theme(axis.text.y = element_markdown(color = NA))
    })

# wrap_plots(weak_strong_bar_list, ncol = 1, guides = "collect", axis_titles = "collect")

if (.overwrite) {
    for (n in names(weak_strong_bar_list)) {
        save_plot(sprintf("_plots/extremes/barplots-%s.pdf", n),
                  weak_strong_bar_list[[n]] + illustrator_theme,
                  width = 3.15, height = 1.5)
    }; rm(n)
}





# ============================================================================*
# Empirical zeta ----
# ============================================================================*


# In Ives et al. (1999), they found that parasitoids spent ~3.76 more time
# foraging at plants where they encountered an aphid.
# Bestrong simulates our model with varying zeta, then compares the observed
# wasp abundances to those predicted when parasitoids never encounter an aphid
# on a Pseudomonas-inhabited plant but always do on plants without Pseudomonas.

set.seed(1180209329)
zeta_weak_strong_sims <- map(c(0.5, 0.6, 0.7, 0.8, 0.9), \(zeta) {
        sims <- run_sim_combos(wasp_resp = "weak", n_pseudo = 3L, n_sims = 12L,
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



empir_zeta_p <- zeta_weak_strong_sims |>
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



