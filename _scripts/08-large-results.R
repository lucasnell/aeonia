#'
#' Plots for larger landscape simulations
#'

source("_scripts/00-preamble.R")


# Plotting and data functions:
source("_scripts/08b-large-plot-funs.R")

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
# from 06-large-main.sh:
sim_df <- process_large_sim_files("main", 2120927824)









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

large_outcomes_p <- baseline_plotter(outcomes = "n_infected", obs_breaks = 0:2 * 2000)
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

# # If you want to combine by `wasp_resp`:
#
# # Effects (or lack thereof) of clustering vs uniform Pseudomonas:
# large_manip_clust_p <- tibble(wr = levels(wasp_resp_fct)) |>
#     # Used for tags:
#     mutate(tg = LETTERS[1:n()]) |>
#     pmap(\(wr, tg) {
#         p <- baseline_plotter(outcomes = "n_infected", col_fct = "wt_pp",
#                               lty_fct = "wt_vp",
#                               color_vals = c("black", "dodgerblue"),
#                               incl_vals = TRUE,
#                               obs_breaks = 0:2 * 2000, obs_max = 4900,
#                               data_df = sim_df |> filter(wasp_resp == wr),
#                               multiline_col_title = FALSE,
#                               p_tag = tg,
#                               p_title = scenario_title(wr, TRUE, TRUE))
#         if (wr != levels(wasp_resp_fct)[1]) {
#             p <- p + theme(axis.text.y = element_blank())
#         }
#         return(p)
#     }) |>
#     # add_top_labels(add_bot_labs = FALSE) |>
#     (\(x) c(guide_area(), x))() |>
#     wrap_plots(design = "AAA\nB#C",
#                guides = "collect", axis_titles = "collect",
#                widths = c(1, 0.05, 1), heights = c(0.2, 1)) &
#     theme(plot.tag.location = "panel",
#           plot.tag.position = c(0.05, 1.05),
#           legend.position = "top", legend.title.position = "top")



# large_manip_clust_p

if (.overwrite) {
    save_plot("_plots/large-baseline-wt_pp.pdf",
              large_manip_clust_p, width = 6.5, height = 4)
}



large_manip_plots <- crossing(cf = c("pseudo_repel", "virus_attract", "sd_N") |>
                                  (\(x) factor(x, levels = x))(),
                              wr = wasp_resp_fct,
                              onp = c(FALSE, TRUE)) |>
    # Bc factors behave weird sometimes (but are needed for sorting crossing):
    mutate(across(where(is.factor), paste)) |>
    pmap(\(cf, wr, onp) {
        col_pal <- c("black", par_pal[[cf]])
        if (cf == "virus_attract") {
            col_pal <- c("black", lighten(par_pal[[cf]], 0.4), par_pal[[cf]])
        }
        non_defs <- list()
        if (onp) {
            non_defs <- list(wt_vp = "on *Pseudo.*")
        }
        p <- baseline_plotter(outcomes = "n_infected", col_fct = cf,
                              color_vals = col_pal, non_defaults = non_defs,
                              obs_breaks = 0:2 * 3000, obs_max = 7400,
                              data_df = sim_df |> filter(wasp_resp == wr))
        if (wr != levels(wasp_resp_fct)[1] || onp) {
            p <- p + theme(axis.text.y = element_blank())
        }
        if (cf != "sd_N") {
            p <- p + theme(axis.text.x = element_blank())
        }
        if (wr == tail(levels(wasp_resp_fct), 1) && onp) {
            if (cf == "sd_N") out <- list(p)
            else out <- c(list(p), rep(list(plot_spacer()), 7))
        } else out <- list(p, plot_spacer())
        return(out)
    }) |>
    do.call(what = c)

# wrap_plots(large_manip_plots, ncol = 7, guides = "collect", axes = "collect",
#            widths = c(1, 0.05, 1, 0.2, 1, 0.05, 1),
#            heights = c(1, 0.05, 1, 0.05, 1)) &
#     illustrator_theme



# # If you want to combine by `wasp_resp`:
#
# large_manip_plots <- crossing(cf = c("pseudo_repel", "virus_attract", "sd_N") |>
#                                   (\(x) factor(x, levels = x))(),
#                               wr = wasp_resp_fct) |>
#     # Bc factors behave weird sometimes (but are needed for sorting crossing):
#     mutate(across(where(is.factor), paste)) |>
#     pmap(\(cf, wr) {
#         col_pal <- c("black", par_pal[[cf]])
#         obs_max <- ifelse(wr == levels(wasp_resp_fct)[1], 7400, 4200)
#         obs_breaks <- ifelse(wr == levels(wasp_resp_fct)[1],
#                              list(0:2 * 3000), list(0:2 * 2000))[[1]]
#         p <- baseline_plotter(outcomes = "n_infected", col_fct = cf,
#                               lty_fct = "wt_vp",
#                               color_vals = col_pal,
#                               obs_breaks = obs_breaks, obs_max = obs_max,
#                               data_df = sim_df |> filter(wasp_resp == wr))
#         # if (wr != levels(wasp_resp_fct)[1]) {
#         #     p <- p + theme(axis.text.y = element_blank())
#         # }
#         if (cf != "sd_N") {
#             p <- p + theme(axis.text.x = element_blank())
#         }
#         if (wr == tail(levels(wasp_resp_fct), 1)) {
#             if (cf == "sd_N") out <- list(p)
#             else out <- c(list(p), rep(list(plot_spacer()), 3))
#         } else out <- list(p, plot_spacer())
#         return(out)
#     }) |>
#     do.call(what = c)
#
# wrap_plots(large_manip_plots, ncol = 3, guides = "collect", axis_titles = "collect",
#            widths = c(1, 0.2, 1),
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


# ----------------------------------*
## Plots of interactions, if interested:
# ----------------------------------*

# baseline_plotter(outcomes = "n_infected",
#                  col_fct = c("wt_vp", "pseudo_repel"),
#                  color_vals = c("black", "blue", "red", "purple"),
#                  inters = TRUE)
# baseline_plotter(outcomes = "n_infected",
#                  col_fct = c("wt_vp", "virus_attract"),
#                  color_vals = c("black", "blue", "red", "purple"),
#                  inters = TRUE)







# =============================================================================*
# =============================================================================*
# Small outbreaks ----
# =============================================================================*
# =============================================================================*


large_small_outcomes_p <- baseline_plotter(outcomes = "all",
                                           non_defaults = list(p_load = 0.05),
                                           multiline_ylab = TRUE,
                                           p_tag = as.list(LETTERS[1:6]),
                                           p_title = c(map(levels(wasp_resp_fct),
                                                           \(wr) scenario_title(wr, TRUE, TRUE)),
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


large_small_outs_main_p <- crossing(cf = c("pseudo_repel", "virus_attract", "sd_N") |>
             (\(x) factor(x, levels = x))(),
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
                             cf == "virus_attract" ~ 36.5,
                             cf == "sd_N" ~ 7)
        obs_breaks <- case_when(cf == "pseudo_repel" ~ list(c(1, 6, 11)),
                                cf == "virus_attract" ~ list(c(1, 18, 35)),
                                cf == "sd_N" ~ list(c(1, 4, 7)))[[1]]
        if (ou == "outbreak_size") obs_breaks[1] <- 2
        p <- baseline_plotter(outcomes = ou, col_fct = cf,
                              color_vals = col_pal, non_defaults = non_defs,
                              multiline_ylab = TRUE,
                              obs_breaks = obs_breaks,
                              obs_max = obs_max,
                              data_df = sim_df |> filter(wasp_resp == wr),
                              p_tag = list(tg))
        stopifnot(is_ggplot(p))
        if (wr != levels(wasp_resp_fct)[1] || onp) {
            p <- p + theme(axis.text.y = element_blank())
        }
        if (!(cf == "sd_N" && ou == "n_infected")) {
            p <- p + theme(axis.text.x = element_blank())
        }
        if (cf == "pseudo_repel" && ou == "p_emerge") {
            p <- p + labs(title = paste0("Virus starts on<br>",
                                         ifelse(onp, "*Pseudomonas*", "uninhabited"),
                                         "<br>plant")) +
                theme(plot.title = element_markdown(size = 10, lineheight = 0.8))
        }
        pl <- list(p)
        out <- pl
        # Last column of last outcome (except for last row) gets vertical spacers:
        if (wr == tail(levels(wasp_resp_fct), 1) && onp) {
            if (ou == "n_infected") {
                if (cf == "sd_N") out <- pl
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
    theme(axis.title.y = element_markdown(hjust = 1, vjust = 0.5, angle = 0),
          plot.margin = margin(0,0,0,0))

large_small_outs_p <- function() {
    grid.newpage()
    grid.draw(patchworkGrob(large_small_outs_main_p))
    pushViewport(viewport(x = 0.17, y = 1, width = 1/3, height = 0.1, name = "top-left",
                          just = c("left", "top")))
    grid.draw(richtext_grob(scenario_title(levels(wasp_resp_fct)[1], TRUE, TRUE),
                            gp = gpar(fontsize = 13, lineheight = 0.8)))
    popViewport()
    pushViewport(viewport(x = 0.54, y = 1, width = 1/3, height = 0.1, name = "top-right",
                          just = c("left", "top")))
    grid.draw(richtext_grob(scenario_title(levels(wasp_resp_fct)[2], TRUE, TRUE),
                            gp = gpar(fontsize = 13, lineheight = 0.8)))
    popViewport()
}


if (.overwrite) {
    save_plot("_plots/large-small-manips.pdf", large_small_outs_p,
              width = 8, height = 9)
}








# =============================================================================*
# =============================================================================*
# Interior virus locations ----
# =============================================================================*
# =============================================================================*


# from 07b-large-interior.sh:
int_sim_df <- process_large_sim_files("interior", 1309766214)


large_int_manip_p <- crossing(cf = c("pseudo_repel", "virus_attract", "sd_N") |>
                                  (\(x) factor(x, levels = x))(),
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
        if (cf != "sd_N") {
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


# =============================================================================*
# =============================================================================*
# Empirical zeta ----
# =============================================================================*
# =============================================================================*





# In Ives et al. (1999), they found that parasitoids spent ~3.76 more time
# foraging at plants where they encountered an aphid.
# Bestrong simulates our model with varying zeta, then compares the observed
# wasp abundances to those predicted when parasitoids never encounter an aphid
# on a Pseudomonas-inhabited plant but always do on plants without Pseudomonas.

if (!file.exists(interm_files$large_zeta_sims)) stop("Run 11-large-emp-zeta.R first!")

zeta_sims <- read_csv(interm_files$large_zeta_sims, col_types = "diicidddddld")

add_n_pseudo_fct <- function(d, label_both = FALSE) {
    fmt <- ifelse(label_both, "%.0f%% *Pseudo.*", "%.0f%%")
    n_pseudo_labs <- n_pseudo_lvls |>
        (\(x) x / 10e3 * 100)() |>
        (\(x) sprintf(fmt, x))()
    d |>
        mutate(n_pseudo_fct = factor(n_pseudo, levels = n_pseudo_lvls,
                                     labels = n_pseudo_labs))
}

# Takes ~ 2 min
set.seed(1540192361)
zeta_sim_summs <- zeta_sims |>
    mutate(rel = wasps / pred_wasps) |>
    group_by(zeta, n_pseudo, pseudo) |>
    summarize(ci = list(booter(rel)[c("Lower", "Upper")] |>
                            as.list() |> as_tibble()),
              min = min(rel), max = max(rel),
              rel = mean(rel), .groups = "drop") |>
    unnest(ci) |>
    mutate(pseudo = factor(pseudo, levels = c(FALSE, TRUE),
                           labels = c("Plants without<br>*Pseudomonas*",
                                      "Plants with<br>*Pseudomonas*"))) |>
    add_n_pseudo_fct(TRUE)


# approximate using linear interpolation where lines equal 1:
zeta_one_ests <- zeta_sim_summs |>
    filter(n_pseudo > 0) |>
    group_by(n_pseudo_fct, n_pseudo, pseudo) |>
    summarize(zeta = (\(z, r) {
        # Linear interpolation:
        f <- approxfun(z, r)
        # Find where it equals 1:
        approx_z <- uniroot(function(x) f(x) - 1, range(z))$root
        return(approx_z)
    })(zeta, rel), .groups = "drop_last") |>
    summarize(zeta_d = abs(diff(range(zeta))),
              zeta = mean(zeta), .groups = "drop")

# They are the same across whether Pseudomonas is on plant:
zeta_one_ests$zeta_d |> max()
# [1] 1.332268e-15

zeta_one_ests |> select(n_pseudo, zeta)
#   n_pseudo  zeta
#      <int> <dbl>
# 1     1000 0.749
# 2     3000 0.698
# 3     5000 0.625
# 4     7000 0.509
# 5     9000 0.296


zeta_p <- zeta_sim_summs |>
    filter(n_pseudo > 0) |>
    ggplot(aes(zeta, rel, color = pseudo)) +
    geom_hline(yintercept = 1, color = "gray70", linewidth = 1) +
    geom_vline(data = zeta_one_ests, aes(xintercept = zeta),
               color = "black", linetype = "22") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = pseudo),
    # geom_ribbon(aes(ymin = min, ymax = max, fill = pseudo),
                color = NA, alpha = 0.25) +
    geom_line(linewidth = 0.75) +
    facet_wrap(~ n_pseudo_fct, nrow = 2, axes = "all", axis.labels = "margins") +
    labs(x = pretty_params("zeta", cap1 = TRUE),
         y = "Observed / predicted parasitoid density") +
    scale_color_manual(NULL, values = np_pal |> set_names(nm = NULL),
                          aesthetics = c("color", "fill")) +
    theme(legend.position = c(5/6, 1/4),
          legend.justification = c(0.5, 0.5),
          legend.key.spacing.y = unit(1, "lines"))

# zeta_p

if (.overwrite) {
    save_plot("_plots/large-empirical-zeta.pdf", zeta_p, width = 6, height = 5)
}



# Now plot for zeta ~ n_infected plants
set.seed(19813000)
zeta_ninf_sim_summs <- zeta_sims |>
    filter(n_pseudo > 0) |>
    group_by(zeta, n_pseudo, rep) |>
    summarize(n_infected = sum(max_virus), .groups = "drop_last") |>
    summarize(ci = list(booter(n_infected)[c("Lower", "Upper")] |>
                            as.list() |> as_tibble()),
              n_infected = mean(n_infected), .groups = "drop") |>
    unnest(ci)

set.seed(694769807)
no_pseudo_df <- zeta_sims |>
    filter(n_pseudo == 0) |>
    group_by(zeta, rep) |>
    summarize(n_infected = sum(max_virus), .groups = "drop") |>
    summarize(ci = list(booter(n_infected)[c("Lower", "Upper")] |>
                            as.list() |> as_tibble()),
              n_infected = mean(n_infected), .groups = "drop") |>
    unnest(ci)





zeta_ninf_p <- zeta_ninf_sim_summs |>
    filter(n_pseudo > 0) |>
    add_n_pseudo_fct() |>
    ggplot(aes(zeta, n_infected, color = n_pseudo_fct)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = n_pseudo_fct), alpha = 0.25, color = NA) +
    geom_line() +
    geom_hline(yintercept = no_pseudo_df$n_infected) +
    geom_hline(yintercept = as.numeric(no_pseudo_df[,c("Lower", "Upper")]),
               linetype = "22") +
    geom_point(data = zeta_one_ests |>
                   # Use linear interpolation to estimate value of n_infected for each zeta:
                   mutate(n_infected = map2_dbl(n_pseudo, zeta, \(np, z) {
                       d <- zeta_ninf_sim_summs |> filter(n_pseudo == np)
                       f <- approxfun(d$zeta, d$n_infected)
                       return(f(z))
                   })) |>
                   add_n_pseudo_fct(),
               size = 3) +
    scale_color_manual("*Pseudo.*<br>plants",
                       values = full_np_pal[paste(n_pseudo_lvls)] |>
                           set_names(sprintf("%.0f%%", n_pseudo_lvls / 10e3 * 100)),
                       aesthetics = c("color", "fill")) +
    labs(x = pretty_params("zeta", cap1 = TRUE),
         y = "Mean peak infected plants")

# zeta_ninf_p


# # What proportion are above # infected when no Pseudomonas?
# zeta_ninf_sim_summs |>
#     getElement("n_infected") |>
#     (\(x) mean(x > no_pseudo_df$n_infected))()
# # [1] 0.4222222
#
# zeta_ninf_sim_summs |>
#     filter(n_pseudo < 9000) |>
#     getElement("n_infected") |>
#     (\(x) mean(x > no_pseudo_df$n_infected))()
# # [1] 0.5
#
# zeta_ninf_sim_summs |>
#     filter(n_pseudo == 5000) |>
#     getElement("n_infected") |>
#     (\(x) mean(x > no_pseudo_df$n_infected))()
# # [1] 0.4444444

if (.overwrite) {
    save_plot("_plots/large-zeta-n_infected.pdf", zeta_ninf_p, width = 6, height = 4)
}
