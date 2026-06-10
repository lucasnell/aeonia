#'
#' Plots for larger landscape simulations
#'

source("_scripts/00-preamble.R")
library("grid")
library("gridtext")

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
# 1 weak        0.05        0       2.42
# 2 weak        0.05     7000       4.46
# 3 weak        0.5         0    2020.
# 4 weak        0.5      7000    4612.
# 5 strong      0.05        0       2.55
# 6 strong      0.05     7000       1.51
# 7 strong      0.5         0    2041.
# 8 strong      0.5      7000     452.









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
        non_defs <- list()
        col_pal <- c("black", par_pal[[cf]])
        if (onp) {
            non_defs <- list(wt_vp = "on *Pseudo.*")
            # col_pal[1] <- par_pal[["wt_vp"]]
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
    }) |>
    do.call(what = c)


# wrap_plots(large_manip_plots, ncol = 7, guides = "collect", axes = "collect",
#            widths = c(1, 0.05, 1, 0.2, 1, 0.05, 1),
#            heights = c(1, 0.05, 1, 0.05, 1)) &
#     illustrator_theme

if (.overwrite) {
    save_plot("_plots/large-baselines.pdf",
              wrap_plots(large_manip_plots, ncol = 7, guides = "collect", axes = "collect",
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


large_small_outs_plots <- c("off", "on") |>
    set_names() |>
    map(\(pos) {
        wt_vp <- sprintf("%s *Pseudo.*", pos)
        title <- paste("Virus starts on",
                       ifelse(pos == "off", "uninhabited", "*Pseudomonas*"),
                       "plant")
        c(`1` = "pseudo_repel", `3` = "virus_attract", `5` = "sd_N") |>
            imap(\(n, i) {
                i <- as.integer(i)
                titles <- list(waiver())
                if (i == 1) {
                    titles <- c(map(levels(wasp_resp_fct),
                                    \(wr) scenario_title(wr, TRUE, TRUE)),
                                rep(list(waiver()), 4))
                }
                non_defs <- list(p_load = 0.05)
                if (pos == "on") non_defs <- c(non_defs, list("wt_vp" = wt_vp))
                pl <- baseline_plotter(outcomes = "all", col_fct = n,
                                       color_vals = c("black", par_pal[[n]]),
                                       non_defaults = non_defs,
                                       multiline_ylab = TRUE,
                                       p_tag = c(as.list(LETTERS[i:(i+1)]),
                                                 rep(list(waiver()), 4)),
                                       p_title = titles,
                                       obs_breaks = scales::breaks_extended(n = 5))
                if (i < 5) {
                    pl <- pl & theme(axis.title.x.bottom = element_blank(),
                                     axis.text.x = element_blank())
                    pl <- list(pl, plot_spacer())
                } else pl <- list(pl)
                return(pl)
            }) |>
            do.call(what = c) |>
            wrap_plots(ncol = 1, guides = "collect", axes = "collect",
                       heights = c(1, 0.05, 1, 0.05, 1))  +
            plot_annotation(title = title,
                            theme = theme(plot.title = element_markdown(
                                face = "bold", size = 16))) &
            theme(plot.tag.location = "panel",
                  plot.tag.position = c(0.05, 1.05),
                  axis.title.y = element_markdown(size = 10, angle = 0,
                                                  hjust = 1, vjust = 0.5))
    })


# large_small_outs_plots



if (.overwrite) {
    for (n in names(large_small_outs_plots)) {
        f <- sprintf("_plots/large-small-manips-%s-pseudo.pdf", n)
        save_plot(f, large_small_outs_plots[[n]], width = 6.5, height = 9)
    }; rm(n, f)
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
# w=1 ----
# =============================================================================*
# =============================================================================*


# from 07-large-w1.sh:
w1_sim_df <- process_large_sim_files("w1", 1627402402)


large_w1_manip_p <- crossing(wr = wasp_resp_fct,
                             onp = c(FALSE, TRUE)) |>
    # Bc factors behave weird sometimes (but are needed for sorting crossing):
    mutate(across(where(is.factor), paste)) |>
    # Used for tags:
    mutate(tg = LETTERS[1:n()]) |>
    pmap(\(wr, onp, tg) {
        cf <- "virus_attract"
        non_defs <- list()
        col_pal <- c("black", par_pal[[cf]], lighten(par_pal[[cf]], 0.4))
        if (onp) {
            non_defs <- list(wt_vp = "on *Pseudo.*")
        }
        p <- baseline_plotter(outcomes = "n_infected", col_fct = cf,
                              color_vals = col_pal, non_defaults = non_defs,
                              obs_breaks = 0:2 * 3000, obs_max = 7780,
                              data_df = w1_sim_df |> filter(wasp_resp == wr),
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


if (.overwrite) {
    save_plot("_plots/large-manips-w1.pdf", large_w1_manip_p,
              width = 6.5, height = 4)
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

zeta_sims <- read_csv(interm_files$large_zeta_sims, col_types = "diicidddddl",
                      progress = FALSE)

# Takes ~ 10 sec
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
                                      "Plants with<br>*Pseudomonas*")),
           n_pseudo_fct = factor(n_pseudo,
                                 levels = n_pseudo_lvls[n_pseudo_lvls > 0],
                                 labels = n_pseudo_lvls[n_pseudo_lvls > 0] |>
                                     (\(x) x / 10e3 * 100)() |>
                                     (\(x) sprintf("%.0f%% *Pseudo.*", x))()))


# approximate using linear interpolation where lines equal 1:
zeta_one_ests <- zeta_sim_summs |>
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

zeta_one_ests |> select(n_pseudo, zeta)
#   n_pseudo  zeta
#      <int> <dbl>
# 1     1000 0.749
# 2     3000 0.698
# 3     5000 0.626
# 4     7000 0.509
# 5     9000 0.296


zeta_p <- zeta_sim_summs |>
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
