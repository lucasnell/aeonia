#'
#' Plots for larger landscape simulations.
#' This includes results from `01-large-main.sh` and `02-large-interior.sh`
#' inside `_scripts/01-large-hpc`.
#' It also includes the plot of total aphid abundance through time,
#' with and without variation among plants in initial aphid abundances.
#'

source("_scripts/00-preamble.R")


# Plotting and data functions:
source("_scripts/03-large/00-large-plot-funs.R")

.overwrite <- FALSE


#' Because below returns an empty tibble, it shows that when we set
#' `wt_vp = 1e-6` it does indeed result in no overlap between
#' *Pseudomonas* and virus.
#' Unsurprisingly, when we set `wt_vp = 100` and manually have them overlap,
#' there is always overlap.
#'
# list.files("_scripts/interm-data", "large-main-.?.?.rds",
#            full.names = TRUE) |>
#     map(read_rds) |>
#     list_rbind() |>
#     select(n_pseudo:wt_pp, landscape) |>
#     filter(!is.na(wt_pp)) |>
#     # Predicted number of instances of both Pseudomonas and virus
#     # across 100 sims:
#     mutate(n_both_pred = ifelse(wt_vp < 1, 0L, 100L)) |>
#     # Observed:
#     mutate(n_both_obs = map_int(landscape, \(x) sum(x == 3L))) |>
#     select(n_pseudo:wt_pp, starts_with("n_both")) |>
#     filter(n_both_pred != n_both_obs)




# Process groups of large landscape sim files:
process_large_sim_files <- function(type, .seed) {
    set.seed(.seed) # for bootstrapping
    list.files(paste0("_scripts/interm-data/large-", type),
               sprintf("large-%s-.?.?.rds", type),
               full.names = TRUE) |>
        map(\(x) {
            read_rds(x) |>
                # Remove this vector immediately bc it's quite large
                select(-landscape)
        }) |>
        list_rbind() |>
        mutate(outbreak_size = map_dbl(sim, \(x) mean(x$n_infected[x$n_infected > 1])),
               log_outbreak_size = map_dbl(sim, \(x) mean(log10(x$n_infected[x$n_infected > 1]))),
               p_emerge = map_dbl(sim, \(x) mean(x$n_infected > 1)),
               n_infected = map_dbl(sim, \(x) mean(x$n_infected)),
               boots = map(sim, \(x) ci_booter(x$n_infected, "all"))) |>
        mutate(wt_vp = ifelse(wt_vp < 1, "off *Pseudo.*", "on *Pseudo.*") |>
                   factor(levels = c("off *Pseudo.*", "on *Pseudo.*")),
               wt_pp = ifelse(wt_pp > 1, "clustered", "uniform") |>
                   factor(levels = c("uniform", "clustered")),
               wasp_resp = factor(wasp_resp, levels = levels(wasp_resp_fct))) |>
        add_no_pseudo_points()
}



# Takes a few seconds
# from 01-large-main.sh:
sim_df <- process_large_sim_files("main", 2120927824)




# Factor for levels of parameters to show in main text:
par_lvls <- c("sd_N", "pseudo_repel", "virus_attract") |>
    (\(x) factor(x, levels = x))()




#' Basic check to make sure that the two scenarios still do what we expect
#' them to: Pseudomonas increases # infected plants when wasp response is
#' weak and decreases when strong

sim_df |>
    filter(sd_N == 0 & virus_attract == 1 & pseudo_repel == 1 &
               wt_vp == "off *Pseudo.*" & wt_pp == "uniform" &
               n_pseudo %in% c(0L, 7e3L)) |>
    select(wasp_resp, p_load, n_pseudo, n_infected) |>
    arrange(wasp_resp, p_load, n_pseudo)
#   wasp_resp p_load n_pseudo n_infected
#   <fct>      <dbl>    <int>      <dbl>
# 1 weak        0.05        0       2.58
# 2 weak        0.05     7000       4.23
# 3 weak        0.5         0    2151.
# 4 weak        0.5      7000    4574.
# 5 strong      0.05        0       2.55
# 6 strong      0.05     7000       1.46
# 7 strong      0.5         0    2100.
# 8 strong      0.5      7000     454.







# =============================================================================*
# =============================================================================*
# Basic plot ----
# =============================================================================*
# =============================================================================*

large_outcomes_p <- baseline_plotter(outcomes = "n_infected", col_fct = "wt_vp",
                                     color_vals = c("black", "#4EEE94"),
                                     obs_breaks = 0:2 * 2000)
# large_outcomes_p





if (.overwrite) {
    save_plot("_plots/large-outcomes.pdf",
              large_outcomes_p & illustrator_theme, width = 4.88, height = 1.5)
}




# =============================================================================*
# =============================================================================*
# Manips ----
# =============================================================================*
# =============================================================================*

#' - `wt_vp`
#' - `pseudo_repel`
#' - `virus_attract`
#' - `sd_N`


# Effects (or lack thereof) of clustering vs uniform Pseudomonas:
large_manip_clust_p <- crossing(wr = wasp_resp_fct,
         onp = c(FALSE, TRUE)) |>
    # Bc factors behave weird sometimes (but are needed for sorting crossing):
    mutate(across(where(is.factor), paste)) |>
    # Used for tags:
    mutate(tg = LETTERS[1:n()]) |>
    pmap(\(wr, onp, tg) {
        non_defs <- list()
        if (onp) {
            non_defs <- list(wt_vp = "on *Pseudo.*")
        }
        p <- baseline_plotter(outcomes = "n_infected", col_fct = "wt_pp",
                              color_vals = c("black", "dodgerblue"),
                              incl_vals = TRUE,
                              non_defaults = non_defs,
                              obs_breaks = 0:2 * 2000, obs_max = 4900,
                              data_df = sim_df |> filter(wasp_resp == wr),
                              multiline_col_title = FALSE,
                              p_tag = tg)
        if (wr != levels(wasp_resp_fct)[1] || onp) {
            p <- p + theme(axis.text.y = element_blank())
        }
        return(p)
    }) |>
    add_top_labels() |>
    wrap_plots(design = "EEE#FFF\nG#H#I#J\nA#B#C#D",
               guides = "collect", axis_titles = "collect",
               widths = c(1, 0.05, 1, 0.2, 1, 0.05, 1), heights = c(0.4, 0.3, 1)) &
    theme(plot.tag.location = "panel",
          plot.tag.position = c(0.05, 1.05),
          legend.position = "top", legend.title.position = "top")

# large_manip_clust_p

if (.overwrite) {
    save_plot("_plots/large-baseline-wt_pp.pdf",
              large_manip_clust_p, width = 6.5, height = 4)
}



large_manip_plots <- crossing(cf = par_lvls,
                              wr = wasp_resp_fct,
                              onp = c(FALSE, TRUE)) |>
    # Bc factors behave weird sometimes (but are needed for sorting crossing):
    mutate(across(where(is.factor), paste)) |>
    pmap(\(cf, wr, onp) {
        col_pal <- c("black", par_pal[[cf]])
        non_defs <- list()
        if (onp) {
            non_defs <- list(wt_vp = "on *Pseudo.*")
        }
        .data_df <- sim_df |> filter(wasp_resp == wr)
        if (cf == "virus_attract") {
            .data_df <- .data_df |> filter(.data[[cf]] < max(.data[[cf]]))
        }
        p <- baseline_plotter(outcomes = "n_infected", col_fct = cf,
                              color_vals = col_pal, non_defaults = non_defs,
                              obs_breaks = 0:2 * 3000, obs_max = 7400,
                              data_df = .data_df)
        if (wr != levels(wasp_resp_fct)[1] || onp) {
            p <- p + theme(axis.text.y = element_blank())
        }
        if (cf != tail(levels(par_lvls), 1)) {
            p <- p + theme(axis.text.x = element_blank())
        }
        if (wr == tail(levels(wasp_resp_fct), 1) && onp) {
            if (cf == tail(levels(par_lvls), 1)) out <- list(p)
            else out <- c(list(p), rep(list(plot_spacer()), 7))
        } else out <- list(p, plot_spacer())
        return(out)
    }) |>
    do.call(what = c)

# wrap_plots(large_manip_plots, ncol = 7, guides = "collect", axes = "collect",
#            widths = c(1, 0.05, 1, 0.2, 1, 0.05, 1),
#            heights = c(1, 0.05, 1, 0.05, 1)) &
#     illustrator_theme




if (.overwrite) {
    save_plot("_plots/large-baselines.pdf",
              wrap_plots(large_manip_plots, ncol = 7, guides = "collect",
                         axis_titles = "collect",
                         widths = c(1, 0.05, 1, 0.2, 1, 0.05, 1),
                         heights = c(1, 0.05, 1, 0.05, 1)) &
                  illustrator_theme,
              width = 6.5, height = 4)
}








# =============================================================================*
# =============================================================================*
# Small outbreaks ----
# =============================================================================*
# =============================================================================*


large_small_outcomes_p <- baseline_plotter(outcomes = "all",
                                           non_defaults = list(p_load = 0.05),
                                           multiline_ylab = FALSE,
                                           p_tag = as.list(LETTERS[1:6]),
                                           p_title = c(map(levels(wasp_resp_fct),
                                                           \(wr) scenario_title(wr)),
                                                       rep(list(waiver()), 4)),
                                           return_list = TRUE) |>
    wrap_plots(design = "A#B\n###\nC#D\n###\nE#F",
               guides = "collect", axis_titles = "collect",
               widths = c(1, 0.05, 1),
               heights = c(1, 0.1, 1, 0.1, 1)) &
    theme(plot.tag.location = "panel",
          plot.tag.position = c(0.05, 0.95))
# large_small_outcomes_p

if (.overwrite) {
    save_plot("_plots/large-outcomes-small-outbreaks.pdf",
              large_small_outcomes_p, width = 5, height = 5)
}


large_small_outs_main_p <- crossing(cf = par_lvls,
                                    ou = c("p_emerge", "outbreak_size", "n_infected") |>
                                        (\(x) factor(x, levels = x))(),
                                    wr = wasp_resp_fct,
                                    onp = c(FALSE, TRUE)) |>
    # Bc factors behave weird sometimes (but are needed for sorting crossing):
    mutate(across(where(is.factor), paste)) |>
    # Used for tags:
    mutate(tg = map(1:n(), \(i) {
        if (((i-1L) %% 12L) <= 3L) {
            j <- i %% 12L
            k <- i %/% 12L * 4L
            return(LETTERS[j+k])
        }
        return(waiver())
    })) |>
    pmap(\(cf, ou, wr, onp, tg) {
        col_pal <- c("black", par_pal[[cf]])
        non_defs <- list(p_load = 0.05)
        if (onp) non_defs <- c(non_defs, list(wt_vp = "on *Pseudo.*"))
        obs_max <- case_when(cf == "pseudo_repel" ~ 11.5,
                             cf == "virus_attract" ~ 36.7,
                             cf == "sd_N" ~ 7)
        obs_breaks <- case_when(cf == "pseudo_repel" ~ list(c(1, 6, 11)),
                                cf == "virus_attract" ~ list(c(1, 18, 35)),
                                cf == "sd_N" ~ list(c(1, 4, 7)))[[1]]
        if (ou == "outbreak_size") obs_breaks[1] <- 2
        .data_df <- sim_df |> filter(wasp_resp == wr)
        if (cf == "virus_attract") {
            .data_df <- .data_df |> filter(.data[[cf]] < max(.data[[cf]]))
        }
        p <- baseline_plotter(outcomes = ou, col_fct = cf,
                              color_vals = col_pal, non_defaults = non_defs,
                              multiline_ylab = TRUE,
                              obs_breaks = obs_breaks,
                              obs_max = obs_max,
                              data_df = .data_df,
                              p_tag = list(tg))
        stopifnot(is_ggplot(p))
        if (wr != levels(wasp_resp_fct)[1] || onp) {
            p <- p + theme(axis.text.y = element_blank())
        }
        if (!(cf == tail(levels(par_lvls),1) && ou == "n_infected")) {
            p <- p + theme(axis.text.x = element_blank())
        }
        # if (cf == levels(par_lvls)[1] && ou == "p_emerge") {
        #     p <- p + labs(title = paste0("Virus starts on<br>",
        #                                  ifelse(onp, "*Pseudomonas*", "uninhabited"),
        #                                  "<br>plant")) +
        #         theme(plot.title = element_markdown(size = 10, lineheight = 0.8))
        # }
        pl <- list(p)
        out <- pl
        # Last column of last outcome (except for last row) gets vertical spacers:
        if (wr == tail(levels(wasp_resp_fct), 1) && onp) {
            if (ou == "n_infected") {
                if (cf == tail(levels(par_lvls),1)) out <- pl
                else out <- c(pl, rep(list(plot_spacer()), 7))
            }
        } else {
            # First two columns get horizontal spacers:
            out <- c(pl, list(plot_spacer()))
        }
        return(out)
    }) |>
    do.call(what = c) |>
    # Add space for titles before first row:
    (\(x) c(rep(list(plot_spacer()), 7), x))() |>
    wrap_plots(ncol = 7, guides = "collect", axis_titles = "collect",
               widths = c(1, 0.01, 1, 0.05, 1, 0.01, 1),
               heights = c(1, head(rep(c(1, 1, 1, 0.05), 3), -1))) &
    theme(# axis.title.y = element_markdown(hjust = 1, vjust = 0.5, angle = 0),
          plot.margin = margin(0,0,0,0))




large_small_outs_p <- function() {

    grid.newpage()

    # Main plot
    pushViewport(viewport(x = 0, y = 0, width = 1, height = 0.95, name = "main",
                          just = c("left", "bottom")))
    grid.draw(patchworkGrob(large_small_outs_main_p))
    popViewport()

    # Top labels
    labs <- scenario_title(levels(wasp_resp_fct))
    top_x <- c(0.315, 0.7)
    for (i in 1:length(top_x)) {
        pushViewport(viewport(x = top_x[i], y = 1, width = 1/3, height = 0.1,
                              name = paste("top", i), just = c("center", "top")))
        grid.draw(richtext_grob(labs[i], gp = gpar(fontsize = 13, lineheight = 0.8)))
        popViewport()
    }

    # Middle labels
    labs <- rep(paste0("Virus starts on<br>", c("uninhabited", "*Pseudomonas*"),
                       "<br>plant"), 2)
    mid_x <- c(top_x[1] + c(-1,1) * 0.09, top_x[2] + c(-1,1) * 0.09)
    for (i in 1:length(mid_x)) {
        pushViewport(viewport(x = mid_x[i], y = 0.9, width = 1/6, height = 0.05,
                              name = paste0("mid", i), just = c("center", "top")))
        grid.draw(richtext_grob(labs[i], gp = gpar(fontsize = 10, lineheight = 0.8)))
        popViewport()
    }

}

# large_small_outs_p()


if (.overwrite) {
    save_plot("_plots/large-small-manips.pdf", large_small_outs_p,
              width = 8, height = 9)
}








# =============================================================================*
# =============================================================================*
# Interior virus locations ----
# =============================================================================*
# =============================================================================*


# from 02-large-interior.sh:
int_sim_df <- process_large_sim_files("interior", 1309766214)


large_int_manip_p <- crossing(cf = par_lvls,
                              wr = wasp_resp_fct,
                              onp = c(FALSE, TRUE)) |>
    # Bc factors behave weird sometimes (but are needed for sorting crossing):
    mutate(across(where(is.factor), paste)) |>
    # Used for tags:
    mutate(tg = LETTERS[1:n()]) |>
    pmap(\(cf, wr, onp, tg) {
        non_defs <- list()
        col_pal <- c("black", par_pal[[cf]])
        if (onp) {
            non_defs <- list(wt_vp = "on *Pseudo.*")
        }
        p <- baseline_plotter(outcomes = "n_infected", col_fct = cf,
                              color_vals = col_pal, non_defaults = non_defs,
                              obs_breaks = 0:2 * 5000, obs_max = 10e3,
                              data_df = int_sim_df |> filter(wasp_resp == wr),
                              p_tag = tg)
        if (wr != levels(wasp_resp_fct)[1] || onp) {
            p <- p + theme(axis.text.y = element_blank())
        }
        if (cf != tail(levels(par_lvls),1)) {
            p <- p + theme(axis.text.x = element_blank())
        }
        return(p)
    }) |>
    do.call(what = c) |>
    add_top_labels() |>
    wrap_plots(design = "MMM#NNN\nO#P#Q#R\nA#B#C#D\n#######\nE#F#G#H\n#######\nI#J#K#L",
               guides = "collect", axis_titles = "collect",
               widths = c(1, 0.05, 1, 0.2, 1, 0.05, 1),
               heights = c(0.5, 0.4, 1, 0.1, 1, 0.1, 1)) &
    theme(plot.tag.location = "panel",
          plot.tag.position = c(0.1, 1))

# large_int_manip_p


if (.overwrite) {
    save_plot("_plots/large-manips-interior.pdf", large_int_manip_p,
              width = 7.5, height = 6)
}







# =============================================================================*
# =============================================================================*
# sd_N ----
# =============================================================================*
# =============================================================================*



N0 <- 55
sd_N <- 50
# Converting to mu and sigma for underlying normal distribution:
mu_N <- log(N0^2 / sqrt(N0^2 + sd_N^2))
sigma_N <- sqrt(log(1 + sd_N^2 / N0^2))


sd_sims_args <- list(max_t = 30L,
                     n_reps = 100L,
                     insects_ptr = make_insects_ptr(pseudo_surv = 1, zeta = 0.5))


set.seed(1573430855)
variable_N_df <- map(1:sd_sims_args$n_reps, \(i) {
    test_insect_pops(N0 = rlnorm(1, mu_N, sigma_N), max_t = sd_sims_args$max_t,
                     W0 = 0, M0 = 0, Y0 = 0,
                     sd_sims_args$insects_ptr) |>
        mutate(rep = i)
}) |>
    list_rbind() |>
    mutate(aphids = aphids + alates + parasitized) |>
    select(rep, time, aphids, wasps)

fixed_N_df <- test_insect_pops(N0 = 55, max_t = sd_sims_args$max_t,
                               W0 = 0, M0 = 0, Y0 = 0,
                               sd_sims_args$insects_ptr) |>
    mutate(aphids = aphids + alates + parasitized) |>
    select(time, aphids, wasps)

sd_N_dens_p <- variable_N_df |>
    ggplot(aes(time, aphids)) +
    geom_line(aes(group = rep), color = "dodgerblue", alpha = 0.1) +
    geom_line(data = variable_N_df |>
                  group_by(time) |>
                  summarize(aphids = mean(aphids)),
              color = "dodgerblue", linewidth = 1) +
    geom_line(data = fixed_N_df, linewidth = 1) +
    labs(x = "Time (days)", y = "Total aphids")

if (.overwrite) {
    save_plot("_plots/sd_N-densities.pdf", sd_N_dens_p, width = 6, height = 4)
}

