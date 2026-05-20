
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
large_sims <- crossing(wasp_resp = factor(1:2, labels = levels(wasp_resp_fct)),
                       n_pseudo = c(0L, 3L)) |>
    mutate(sims = map2(wasp_resp, n_pseudo, \(wasp_resp, n_pseudo) {
        run_sim_combos(wasp_resp = wasp_resp, n_pseudo = n_pseudo, large_sims = TRUE)
    })) |>
    mutate(n_infected = map_dbl(sims, \(x) mean(x$n_infected)),
           p_emerge = map_dbl(sims, \(x) mean(x$n_infected > 1)),
           outbreak_size = map_dbl(sims, \(x) mean(x$n_infected[x$n_infected > 1]))) |>
    mutate(pe_boots = map(sims, \(x) {
        ci <- booter(as.integer(x$n_infected > 1L))
        return(ci[c("Lower", "Upper")] |> as.list() |> as_tibble())
    }),
    ob_boots = map(sims, \(x) {
        z <- x$n_infected[x$n_infected > 1]
        stopifnot(length(z) > 1)
        ci <- booter(z)
        return(ci[c("Lower", "Upper")] |> as.list() |> as_tibble())
    }))



large_sims |>
    select(wasp_resp, n_pseudo, n_infected:outbreak_size) |>
    mutate(across(n_infected:outbreak_size, \(x) num(x, digits = 3)))
#   wasp_resp n_pseudo n_infected  p_emerge outbreak_size
#   <fct>        <int>  <num:.3!> <num:.3!>     <num:.3!>
# 1 weak             0      3.686     0.796         4.374
# 2 weak             3      8.273     1.000         8.273
# 3 strong           0      3.671     0.800         4.339
# 4 strong           3      2.141     0.441         3.587


large_sims |>
    group_by(wasp_resp) |>
    summarize(n_infected = n_infected[n_pseudo != 0] - n_infected[n_pseudo == 0],
              p_emerge = p_emerge[n_pseudo != 0] - p_emerge[n_pseudo == 0],
              outbreak_size = outbreak_size[n_pseudo != 0] - outbreak_size[n_pseudo == 0]) |>
    mutate(across(n_infected:outbreak_size, \(x) num(x, digits = 3)))
#   wasp_resp n_infected  p_emerge outbreak_size
#   <fct>      <num:.3!> <num:.3!>     <num:.3!>
# 1 weak           4.587     0.204         3.899
# 2 strong        -1.530    -0.359        -0.751


# ============================================================================*
# Time series ----
# ============================================================================*


set.seed(1001450726)
weak_strong_sims <- large_sims |>
    select(wasp_resp, n_pseudo) |>
    mutate(sims = map2(wasp_resp, n_pseudo, \(wasp_resp, n_pseudo) {
        # wasp_resp = "weak"; n_pseudo = 3
        # rm(wasp_resp, n_pseudo, sims, target, .rep, half_inf_times, hit_adiffs)
        sims <- run_sim_combos(wasp_resp, n_pseudo, n_sims = 1000L,
                               out_attack_surv = TRUE, out_stages = TRUE)
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
# 2 weak             3       8.28
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


total_plotter <- function(wasp_resp,
                          delay_virus = FALSE,
                          .data_df = weak_strong_sims) {
    # wasp_resp = "weak"; delay_virus = FALSE
    # rm(wasp_resp, delay_virus, dd)

    dd <- .data_df
    if ("wasp_resp" %in% colnames(dd)) dd <- dd |> filter(wasp_resp == .env$wasp_resp)
    if ("plant" %in% colnames(dd)) dd <- dd |> filter(is.na(plant))
    dd <- dd |>
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
        geom_hline(data = filter(tot_empty_data, species == "virus"),
                   aes(yintercept = density), color = "gray70") +
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
            } else if (max(x) >= 80) {# alates (adult)
                return(c(0, 40, 80))
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


by_plant_plotter <- function(wasp_resp, .data_df = weak_strong_sims) {
    # wasp_resp = "weak"
    # rm(wasp_resp, dd)

    dd <- .data_df |>
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
        if (wasp_resp == levels(wasp_resp_fct)[2]) p <- p + theme(axis.text.y = element_blank())
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
# Disease outcomes ----
# ============================================================================*



# -------------------------------------*
# ... Pr(emergence) bar graphs ----
# -------------------------------------*



weak_strong_bar_list <- levels(wasp_resp_fct) |>
    set_names() |>
    map(\(wr) {
        large_sims |>
            filter(wasp_resp == wr)|>
            select(n_pseudo, p_emerge, pe_boots) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            unnest(pe_boots) |>
            ggplot(aes(p_emerge, n_pseudo)) +
            geom_vline(xintercept = 0, color = "gray70") +
            geom_col(aes(fill = n_pseudo), color = NA, width = 0.45,
                     linewidth = 0.75, linejoin = "mitre") +
            geom_errorbar(aes(xmin = Lower, xmax = Upper),
                          width = 0.2, linewidth = 1, orientation = "y") +
            labs(x = "Prob. of emergence", y = "Number of *Pseudo.* plants") +
            coord_cartesian(xlim  = c(0, 1)) +
            scale_fill_manual(values = np_pal) +
            theme(axis.text.y = element_blank())
    })

# wrap_plots(weak_strong_bar_list, nrow = 1, guides = "collect", axis_titles = "collect")


# -------------------------------------*
# ... Outbreak sizes - freq. line graphs ----
# -------------------------------------*

hist_empty_data <- large_sims |>
    mutate(outbreak_size = map(sims, \(x) x$n_infected)) |>
    select(-sims) |>
    unnest(outbreak_size) |>
    filter(outbreak_size > 1) |>
    mutate(outbreak_size = factor(outbreak_size, levels = 2:9)) |>
    group_by(wasp_resp, n_pseudo, outbreak_size) |>
    count(name = "n_obs", .drop = FALSE) |>
    group_by(wasp_resp, n_pseudo) |>
    mutate(perc = 100 * n_obs / sum(n_obs)) |>
    ungroup() |>
    mutate(outbreak_size = as.integer(paste(outbreak_size))) |>
    filter(perc == max(perc)) |>
    select(outbreak_size, perc) |>
    mutate(n_pseudo = list(factor(unique(large_sims$n_pseudo)))) |>
    unnest(n_pseudo)


weak_strong_hist_list <- levels(wasp_resp_fct) |>
    set_names() |>
    map(\(wr) {
        # wr = "strong"
        # rm(wr, dd, dd_means)
        dd  <- large_sims |>
            filter(wasp_resp == wr) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            mutate(outbreak_size = map(sims, \(x) x$n_infected)) |>
            select(n_pseudo, outbreak_size) |>
            unnest(outbreak_size) |>
            filter(outbreak_size > 1) |>
            mutate(outbreak_size = factor(outbreak_size, levels = 2:9)) |>
            group_by(n_pseudo, outbreak_size) |>
            count(name = "n_obs", .drop = FALSE) |>
            group_by(n_pseudo) |>
            mutate(perc = 100 * n_obs / sum(n_obs)) |>
            ungroup() |>
            mutate(outbreak_size = as.integer(paste(outbreak_size)))
        dd_means <- large_sims |>
            filter(wasp_resp == wr) |>
            mutate(n_pseudo = factor(n_pseudo)) |>
            select(n_pseudo, outbreak_size, ob_boots) |>
            unnest(ob_boots)
        p <- dd |>
            ggplot(aes(outbreak_size, perc)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_blank(data = hist_empty_data) +
            geom_point(aes(color = n_pseudo), size = 2, stroke = 1) +
            geom_line(aes(color = n_pseudo), linewidth = 1) +
            geom_vline(data = dd_means,
                       aes(xintercept = outbreak_size, color = n_pseudo),
                       linewidth = 1, linetype = "22") +
            geom_ribbon(data = dd_means |> crossing(perc = c(-Inf, Inf)),
                       aes(y = perc,
                           xmin = Lower, xmax = Upper, fill = n_pseudo),
                       color = NA, alpha = 0.25) +
            scale_color_manual(values = np_pal, aesthetics = c("color", "fill")) +
            scale_x_continuous(breaks = (0:4) * 2 + 1) +
            scale_y_continuous(breaks = \(y) {
                if (max(y) < 65 && max(y) > 50) {
                    return(c(0, 25, 50))
                } else return(scales::breaks_extended(n = 5)(y))
            }) +
            labs(x = "Outbreak size", y = "Percent of simulations")
        if (wr == levels(large_sims$wasp_resp)[[2]])
            p <- p + theme(axis.text.y = element_blank())
        return(p)
    })

# wrap_plots(weak_strong_hist_list, nrow = 1, guides = "collect", axis_titles = "collect")



# -------------------------------------*
# ... stitch together ----
# -------------------------------------*


if (.overwrite) {
    .p <- c(weak_strong_bar_list, weak_strong_hist_list) |>
        wrap_plots(design = "A#B\n###\nC#D", widths = c(1, 0.1, 1),
                   heights = c(1, 0.4, 1), guides = "collect",
                   axis_titles = "collect") & illustrator_theme &
        theme(plot.margin = margin(0,0,0,0))
    save_plot("_plots/extremes/disease-outcomes.pdf", .p,
              width = 6.2, height = 3.825)
    rm(.p)
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
                           labels = c("Plants with<br>*Pseudomonas*",
                                      "Plants without<br>*Pseudomonas*")),
           zeta = factor(zeta),
           rep = factor(rep),
           rel = wasps / pred_wasps) |>
    ggplot(aes(zeta, rel)) +
    geom_jitter(aes(color = zeta), height = 0, width = 0.2, shape = 1, size = 2, alpha = 0.5) +
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



# ============================================================================*
# Density-dep. dispersal ----
# ============================================================================*



# Time series showing alate_p = 0.05 and 0.08 work best for high and low values
# (when comparing to Figure 2 A,B):
#
# Takes just a couple of seconds
#
set.seed(1145518329)
dd_disp_ts_sims <- crossing(alate_p = c(0.05, 0.08),
                         wasp_resp = wasp_resp_fct,
                         n_pseudo = c(0L, 3L)) |>
    pmap(\(alate_p, wasp_resp, n_pseudo) {
        run_sim_combos(wasp_resp = wasp_resp, n_pseudo = n_pseudo,
                       large_sims = TRUE, summ = "time", out_stages = TRUE,
                       alate_b1 = 0, alate_b0 = logit(alate_p)) |>
            mutate(aphids = aphids_juv + aphids_adu + alates_juv + alates_adu + parasitized,
                   alates = alates_adu) |>
            select(rep, time, aphids, alates, wasps, virus) |>
            group_by(time) |>
            summarize(across(aphids:virus, mean), .groups = "drop") |>
            mutate(alate_p = .env$alate_p,
                   wasp_resp = .env$wasp_resp,
                   n_pseudo = .env$n_pseudo) |>
            select(alate_p, wasp_resp, n_pseudo, everything())
    }) |>
    list_rbind() |>
    mutate(alate_p = factor(alate_p))



dd_disp_ts_p <- dd_disp_ts_sims |>
    pivot_longer(aphids:virus, names_to = "species",
                 values_to = "density") |>
    mutate(species = factor(species, levels =
                                levels(tot_empty_data$species))) |>
    left_join(weak_strong_sims |>
                  filter(is.na(plant)) |>
                  select(wasp_resp, n_pseudo, time, aphids, alates, wasps, virus) |>
                  pivot_longer(aphids:virus, names_to = "species",
                               values_to = "density") |>
                  mutate(species = factor(species, levels =
                                              levels(tot_empty_data$species))) |>
                  rename(density_dd = density),
              by = c("wasp_resp", "n_pseudo", "time", "species"),
              relationship = "many-to-one") |>
    pivot_longer(density:density_dd, names_to = "dd_disp",
                 values_to = "density") |>
    mutate(dd_disp = factor(dd_disp, levels = c("density", "density_dd"),
                            labels = c("fixed", "dens.-dep.")),
           wasp_resp = factor(wasp_resp, levels = levels(wasp_resp_fct),
                              labels = scenario_title(levels(wasp_resp_fct),
                                                      TRUE, TRUE)),
           n_pseudo = factor(n_pseudo)) |>
    # mutate(density = ifelse(species == "aphids", density / 100, density),
    #        density_dd = ifelse(species == "aphids", density_dd / 100, density_dd)) |>
    split(~ alate_p) |>
    map(\(xx) {
        dd_empty_df <- xx |>
            group_by(species) |>
            summarize(density = max(density)) |>
            mutate(time = 1,
                   density = ifelse(species == "virus", 9, density))
        p <- xx |>
            ggplot(aes(time, density)) +
            geom_hline(yintercept = 0, color = "gray70") +
            geom_hline(data = filter(dd_empty_df, species == "virus"),
                       aes(yintercept = density), color = "gray70") +
            geom_blank(data = dd_empty_df) +
            geom_line(aes(color = n_pseudo, linetype = dd_disp,
                          linewidth = dd_disp)) +
            labs(x = "Time (days)", y = "Total density",
                 title = sprintf("Fixed alate proportion = %s",
                                 paste(xx$alate_p[[1]]))) +
            scale_x_continuous(breaks = c(0, 25, 50, 75, 100)) +
            scale_y_continuous(breaks = \(x) {
                if (max(x) < 10) return(c(1, 5, 9))
                scales::extended_breaks(n = 4)(x)}) +
            scale_color_manual("*Pseudo.*<br>plants:", values = np_pal) +
            scale_linetype_manual("Alate<br>production:",
                                  values = c("solid", "22")) +
            scale_linewidth_manual("Alate<br>production:",
                                  values = c(0.75, 1)) +
            facet_grid(species ~ wasp_resp, scales = "free_y",
                       switch = "y") +
            theme(strip.text.y.left = element_markdown(
                size = 9, angle = 0, hjust = 1),
                strip.placement = "outside") +
            guides(color = guide_legend(override.aes = list(linewidth = 1)))
        # if (paste(xx$alate_p[[1]]) == levels(xx$alate_p)[1]) {
        #     p <- p +
        #         geom_richtext(data = dd_empty_df[1,c("species", "density")] |>
        #                           crossing(tibble(n_pseudo = sort(unique(xx$n_pseudo)))) |>
        #                           mutate(time = 100, density = density * c(1, 0.2),
        #                                  wasp_resp = tail(sort(unique(xx$wasp_resp)),1),
        #                                  np_fct = fct_recode(n_pseudo,
        #                                                      "0 *Pseudo.*" = "0",
        #                                                      "3 *Pseudo.*" = "3")),
        #                       aes(color = n_pseudo, label = np_fct),
        #                       hjust = 1, vjust = c(1, 0), label.colour = NA, fill = NA,
        #                       fontface = "bold")
        # }
        return(p)
    }) |>
    wrap_plots(nrow = 2, axis_titles = "collect", guides = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"),
          plot.tag.location = c("panel", "plot", "margin")[1],
          plot.tag.position = c(0.05, 0.95),
          legend.position = "bottom", legend.direction = "vertical",
          legend.title.position = "left",
          legend.title = element_markdown(hjust = 1))

# dd_disp_ts_p

if (.overwrite) {
    save_plot("_plots/dens-dep-dispersal-densities.pdf", dd_disp_ts_p,
              width = 6, height = 9)
}



if (.overwrite || !file.exists(rds_files$dd_disp_sims)) {

    # Takes ~1 min
    set.seed(1512618759)
    dd_disp_sims <- crossing(disp = c(NA_real_, 0.05, 0.08),
                             zeta = seq(0, 1, 0.05),
                             n_pseudo = c(0L, 3L)) |>
        pmap(\(disp, zeta, n_pseudo) {
            b0 <- eval(formals(make_insects_ptr)[["alate_b0"]])
            b1 <- eval(formals(make_insects_ptr)[["alate_b1"]])
            if (!is.na(disp)) {
                b0 <- disp |> logit()
                b1 <- 0
            } else disp <- 0 # << used for output
            sims <- run_sim_combos(wasp_resp = "weak", # << this doesn't matter bc of zeta
                                   n_pseudo = n_pseudo,
                                   zeta = zeta,
                                   large_sims = TRUE,
                                   alate_b1 = b1, alate_b0 = b0)
            pe <- as.integer(sims$n_infected > 1)
            ob <- sims$n_infected[sims$n_infected > 1]
            boots <- list(p_emerge = booter(pe),
                          outbreak_size = booter(ob))
            tibble(disp = .env$disp, zeta = .env$zeta,
                   n_pseudo = .env$n_pseudo,
                   p_emerge = mean(sims$n_infected > 1),
                   outbreak_size = mean(sims$n_infected[sims$n_infected > 1]),
                   ci = list(boots))
        }, .progress = .prog_args) |>
        list_rbind()

    write_rds(dd_disp_sims, rds_files$dd_disp_sims, compress = "gz")

} else {

    dd_disp_sims <- read_rds(rds_files$dd_disp_sims)

}




dd_disp_p <- dd_disp_sims |>
    mutate(n_pseudo = factor(n_pseudo),
           disp = factor(disp,
                         labels = c("density<br>dependent",
                                    sprintf("fixed<br>prop = %.2f",
                                            sort(unique(disp))[-1])))) |>
    pivot_longer(p_emerge:outbreak_size, names_to = "outcome") |>
    mutate(outcome = factor(outcome, levels = c("p_emerge", "outbreak_size"))) |>
    filter(!is.na(value)) |>
    mutate(lower = map2_dbl(ci, outcome, \(.c, .o) .c[[.o]][["Lower"]]),
           upper = map2_dbl(ci, outcome, \(.c, .o) .c[[.o]][["Upper"]])) |>
    select(-ci) |>
    split(~ disp + outcome) |>
    map(\(ddd) {
        yvar <- paste(ddd$outcome[[1]])
        .disp <- paste(ddd$disp[[1]])
        yl <- switch(yvar, p_emerge = c(0, 1), outbreak_size = c(2, 9))
        yb <- switch(yvar, p_emerge = 0:2/2, outbreak_size = c(3,6,9))
        # .strip <- switch(yvar, p_emerge = element_markdown(),
        #                  outbreak_size = element_blank())
        .title <- switch(yvar, p_emerge = .disp, outbreak_size = waiver())
        y_lab <- yvar_desc[[yvar]] |>
            (\(x) ifelse(yvar == "outbreak_size", paste("mean", x), x))() |>
            first_cap()
        leg_pos <- ifelse(yvar == "p_emerge" && .disp == tail(levels(ddd$disp),1),
                          list(c(0.5, 0.1)), list("none"))[[1]]
        p <- ddd |>
            ggplot(aes(zeta, value, color = n_pseudo)) +
            geom_hline(yintercept = yl, color = "gray70") +
            geom_vline(xintercept = c(0.1, 0.9), color = "gray70",
                       linetype = "33") +
            geom_ribbon(aes(ymin = lower, ymax = upper, fill = n_pseudo),
                        color = NA, alpha = 0.25) +
            geom_line(linewidth = 1) +
            labs(x = pretty_params("zeta", cap1 = TRUE), y = y_lab,
                 title = .title) +
            scale_x_continuous(breaks = c(0, 0.5, 1)) +
            scale_y_continuous(breaks = yb) +
            scale_color_manual("*Pseudo.*<br>plants", values = np_pal,
                               aesthetics = c("color","fill")) +
            # facet_wrap(~ disp, nrow = 1) +
            coord_cartesian(xlim = range(dd_disp_sims$zeta),
                            ylim = yl * c(1, 1.1)) +
            # theme(strip.text.x = .strip, panel.spacing = unit(1, "lines"))
            theme(legend.position = leg_pos, legend.justification = c(0.5, 0))
        if (.disp != levels(ddd$disp)[1]) {
            p <- p + theme(axis.title.y = element_blank(),
                           axis.text.y = element_blank())
        }
        if (yvar == "p_emerge") {
            p <- p + theme(axis.title.x = element_blank(),
                           axis.text.x = element_blank())
        }
        return(p)
    }) |>
    wrap_plots(nrow = 2, axis_titles = "collect") +
    plot_annotation(tag_levels = "A") &
    theme(plot.tag = element_markdown(face = "bold"),
          plot.tag.location = c("panel", "plot", "margin")[1],
          plot.tag.position = c(0.1, 0.95))

# dd_disp_p

if (.overwrite) {
    save_plot("_plots/dens-dep-dispersal.pdf", dd_disp_p, width = 7.5, height = 6)
}

