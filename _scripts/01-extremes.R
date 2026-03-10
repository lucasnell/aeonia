
#'
#' Time series and histograms for simulations of two scenarios:
#' one where Pseudomonas decreases outbreaks ("low"),
#' and another where it increases outbreaks ("high")
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

if (!file.exists(rds_files$extreme_large) || .overwrite) {

    # Takes very little time, but I'm saving RDS to use output in other scripts.
    set.seed(259619623)
    large_sims <- crossing(type = factor(1:2, labels = c("low", "high")),
                           n_pseudo = c(0L, 3L)) |>
        mutate(sims = map2(type, n_pseudo, \(type, n_pseudo) {
            run_sim_combos(type, n_pseudo, TRUE)
        }))

    write_rds(large_sims, rds_files$extreme_large, compress = "gz")

} else {

    large_sims <- read_rds(rds_files$extreme_large)

}

large_sims |>
    mutate(outbreak_size = num(map_dbl(sims, \(x) mean(x$outbreak_size)), digits = 3)) |>
    select(-sims)
# # A tibble: 4 × 3
#   type  n_pseudo outbreak_size
#   <fct>    <int>     <num:.3!>
# 1 low          0         4.986
# 2 low          3         3.200
# 3 high         0         3.020
# 4 high         3         4.863

large_sims |>
    mutate(outbreak_size = map_dbl(sims, \(x) mean(x$outbreak_size))) |>
    group_by(type) |>
    summarize(outbreak_size = outbreak_size[n_pseudo != 0] -
                  outbreak_size[n_pseudo == 0]) |>
    mutate(outbreak_size = num(outbreak_size, digits = 3))
# # A tibble: 2 × 2
#   type  outbreak_size
#   <fct>     <num:.3!>
# 1 low          -1.786
# 2 high          1.843



# ============================================================================*
# Time series ----
# ============================================================================*


set.seed(1001450726)
low_high_sims <- large_sims |>
    select(-sims) |>
    mutate(sims = map2(type, n_pseudo, \(type, n_pseudo) {
        # type = "low"; n_pseudo = 3
        # rm(type, n_pseudo, sims, target, .rep)
        sims <- run_sim_combos(type, n_pseudo, n_sims = 100L, out_stages = TRUE)
        target <- large_sims |>
            filter(type == .env$type, n_pseudo == .env$n_pseudo) |>
            getElement("sims") |>
            map_dbl(\(x) mean(x$outbreak_size))
        # Choose a representative simulation:
        .rep <- sims |>
            filter(is.na(x)) |>
            group_by(rep) |>
            summarize(outbreak_size = max(virus)) |>
            filter(abs(outbreak_size - target) ==
                       min(abs(outbreak_size - target))) |>
            getElement("rep")
        stopifnot(length(.rep) > 0)
        sims |>
            filter(rep == .rep[1]) |>
            mutate(plant = interaction(y, x),
                   alates = alates_juv + alates_adu,
                   aphids = aphids_juv + aphids_adu + parasitized) |>
            select(plant, time, virus, aphids, alates, wasps,
                   aphids_juv, aphids_adu, alates_juv, alates_adu)
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
            geom_line(aes(linetype = n_pseudo), linewidth = 1, color = "black") +
            geom_blank(data = tot_empty_data) +
            scale_linetype_manual(values = np_linetypes)
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
    select(-virus, -ends_with("_juv"), -ends_with("_adu")) |>
    pivot_longer(aphids:wasps, names_to = "species",
                 values_to = "density") |>
    group_by(species) |>
    summarize(density = max(density)) |>
    mutate(time = 0,
           species = factor(species, levels = c("aphids", "wasps", "alates")),
           density = ifelse(species == "aphids", density,
                            max(density[species != "aphids"])))

dens_maxes <- list(aphids = (max(empty_data$density) %/% 100L) * 100.0,
                   alates = ceiling(min(empty_data$density)),
                   wasps = ceiling(min(empty_data$density)))




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

if (.overwrite) {
    for (n in names(low_high_p_list)) {
        save_plot(sprintf("_plots/extremes/timeseries-%s.pdf", n),
                  low_high_p_list[[n]] & illustrator_theme,
                  width = 5, height = 2.5)
    }; rm(n)
}







by_plant_plotter <- function(x) {
    # x = low_high_sims |> filter(type == "high")
    # rm(x, trans, itrans, vdf, dd, p)
    # wasps / alates --> aphids
    trans <- function(x) {
        x * dens_maxes[["aphids"]] / dens_maxes[["wasps"]]
    }
    # aphids --> wasps / alates
    itrans <- function(x) {
        x * dens_maxes[["wasps"]] / dens_maxes[["aphids"]]
    }

    vdf <- x |>
        filter(!is.na(plant), virus == 1) |>
        group_by(n_pseudo, plant) |>
        summarize(time = min(time), .groups = "drop") |>
        mutate(density = max(empty_data$density) * 1.1,
               n_pseudo = factor(n_pseudo))

    dd <- x |>
        filter(!is.na(plant)) |>
        select(n_pseudo, plant, time, aphids, alates, wasps) |>
        mutate(wasps = trans(wasps),
               alates = trans(alates))

    dd |>
        pivot_longer(aphids:wasps, names_to = "species",
                     values_to = "density") |>
        mutate(n_pseudo = factor(n_pseudo),
               species = factor(species, levels = levels(empty_data$species))) |>
        ggplot(aes(time, density)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(aes(color = species, linetype = n_pseudo), linewidth = 1) +
        geom_blank(data = empty_data, aes(y = density * 1.15)) +
        geom_point(aes(shape = n_pseudo), data = vdf, color = color_pal[["virus"]],
                   size = 2) +
        scale_color_manual(values = color_pal) +
        scale_linetype_manual(values = np_linetypes) +
        scale_shape_manual(values = c("0" = 4, "3" = 13)) +
        labs(x = "Time (days)", y = "Density") +
        scale_y_continuous(sec.axis = sec_axis(itrans,
                                               "Density (wasps/alates)",
                                               breaks = c(0, 20, 40)),
                           breaks = c(0, 500, 1000)) +
        scale_x_continuous(breaks = c(0, 50, 100)) +
        facet_wrap(~ plant) +
        guides(linetype = guide_legend(override.aes = list(color = "black"))) +
        theme(axis.text.x = element_markdown(size = 7),
              axis.text.y = element_markdown(size = 7),
              plot.margin = margin(3,3,0,0))

}










A_plotter <- function(x) {
    # x = low_high_sims |> filter(type == "high")
    # rm(x, A_calc)
    A_calc <- function(wasps, aphids_juv, aphids_adu, alates_juv, alates_adu) {
        # n = 10; wasps = runif(n) * 2; aphids_juv = runif(n) * 10;
        # aphids_adu = runif(n) * 40; alates_juv = runif(n) * 2;
        # alates_adu = runif(n) * 6
        # rm(n, wasps, aphids_juv, aphids_adu, alates_juv, alates_adu)
        # rm(a, h, k, R, X, xt, A)
        a <- pop_info$a
        h <- pop_info$h
        k <- pop_info$k
        R <- pop_info$R
        X <- cbind(aphids_juv, aphids_adu, alates_juv, alates_adu)
        xt <- rowSums(X)
        A <- lapply(1:nrow(X), \(i) {
            .A <- (1 + R * a * wasps[i] / (k * (h * xt[i] + 1)))^(-k)
            return(sum(.A * (X[i,] / xt[i])))
        }) |>
            do.call(what = c)
        return(A)
    }
    x |>
        filter(!is.na(plant)) |>
        mutate(attack_surv = A_calc(wasps, aphids_juv, aphids_adu, alates_juv, alates_adu)) |>
        mutate(n_pseudo = factor(n_pseudo)) |>
        ggplot(aes(time, attack_surv)) +
        geom_line(aes(linetype = n_pseudo), linewidth = 1) +
        scale_linetype_manual(values = np_linetypes) +
        labs(x = "Time (days)", y = "Attack survival") +
        scale_x_continuous(breaks = c(0, 50, 100)) +
        facet_wrap(~ plant)
}


#


c(0L, 3L) |>
    map(\(np) {
        run_sim_combos(type = "high", n_pseudo = np, n_sims = 10L, pseudo_surv = 0.98) |>
            mutate(n_pseudo = np)
    }) |>
    list_rbind() |>
    mutate(plant = interaction(y, x),
           aphids = aphids + parasitized) |>
    filter(rep == sample.int(10L)) |>
    by_plant_plotter()









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
            geom_point(data = dd, aes(shape = n_pseudo), size = 4) +
            scale_shape_manual(values = np_shapes) +
            scale_x_continuous(breaks = c(0, 500, 1000))
    })

wrap_plots(z_pa_p_list, ncol = 1, guides = "collect", axis_titles = "collect")

if (.overwrite) {
    for (n in names(low_high_p_list)) {
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
            mutate(outbreak_size = map(sims, \(x) x$outbreak_size)) |>
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
            geom_point(aes(shape = n_pseudo), size = 2) +
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
            mutate(p_emerge = map_dbl(sims, \(x) mean(x$outbreak_size > 1))) |>
            select(-outbreak_size) |>
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


