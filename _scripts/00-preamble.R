
suppressPackageStartupMessages({
    library(pillar) # num() for setting sig figs in tibble printing
    library(tidyverse)
    library(aeonia)
    library(patchwork)
    library(ggtext)
    library(scico)
    library(viridisLite)
    library(RColorBrewer)
})

# RDS files with simulation output:
rds_files <- list(extreme_manip = "extremes-manip-sims.rds",
                  extreme_manip2 = "extremes-manip2-sims.rds",
                  dens_sims = "density-sims.rds") |>
    map(\(x) paste0("_scripts/interm-data/", x))



# Set threads for simulations:
options("mc.cores" = max(1L, parallel::detectCores()-2L))

# For purrr progress bar:
.prog_args <- list(clear = FALSE,
                   format = paste("{cli::pb_bar}",
                                  "{cli::pb_percent}",
                                  "[{cli::pb_elapsed}] |",
                                  "ETA: {cli::pb_eta}"))

#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists(".LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}



color_pal <- c("#999999", "#0046D2", "#FF64FF",
               viridisLite::inferno(100)[c(50, 25, 80)]) |>
    set_names(c("no_pseudo", "pseudo", "virus", "aphids", "alates", "wasps"))
# scales::show_col(color_pal, labels = FALSE)

# for numbers of pseudomonas patches:
np_pal <- color_pal[c("no_pseudo", "pseudo")] |>
    set_names(c("0", "3"))
# scales::show_col(np_pal, labels = FALSE)

np_shapes <- c(`0` = 1, `3` = 8)
np_linetypes <- c(`0` = "22", `3` = "solid")
np_linewidths <- c(`0` = 1, `3` = 0.75)


serify <- function(prefix, x, suffix) {
    xx <- paste0("<span style=\"font-family: serif\">", x, "</span>")
    paste0(prefix, xx, suffix)
}

# Only first letter capitalized (others not forced to lowercase, as in `str_to_sentence`):
first_cap <- \(str) paste(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)), sep="")




pretty_params <- function(x, short = FALSE, cap1 = FALSE, serif = FALSE) {
    if (short) {
        out <- case_when(x == "pseudo_surv" ~ "<i>&psi;</i>",
                         x == "virus_attract" ~ "<i>&nu;</i>",
                         x == "pseudo_repel" ~ "<i>&rho;</i>",
                         x == "epsilon" ~ "<i>&epsilon;</i>",
                         x == "zeta" ~ "<i>&zeta;</i>",
                         x == "sd_N" ~ "<i>&sigma;</i><sub>N</sub>",
                         x == "Y0" ~ "<i>Y</i><sub>0</sub>",
                         x == "mean_N" ~ "<i>&mu;</i><sub>N</sub>",
                         x == "N0" ~ "<i>N</i><sub>0</sub>",
                         x == "alate_slope" ~ "<i>b</i><sub>slope</sub>",
                         x == "alate_max" ~ "<i>b</i><sub>max</sub>",
                         x == "n_pseudo" ~ "<i>n</i><sub>P</sub>",
                         x == "spat_config" ~ "spat. config.",
                         x == "wt_vp" ~ "virus place.",
                         x == "wt_pp" ~ "*Pseudo.* spacing",
                         .default = x)
        if (serif && any(x != "spat_config")) {
            out[x != "spat_config"] <- serify("", out[x != "spat_config"], "")
        }
    } else {
        out <- case_when(x == "pseudo_surv" ~ serify("*Pseudomonas* survival (", "<i>&psi;</i>", ")"),
                         x == "virus_attract" ~ serify("virus attraction (", "<i>&nu;</i>", ")"),
                         x == "pseudo_repel" ~ serify("*Pseudomonas* repellence (", "<i>&rho;</i>", ")"),
                         x == "epsilon" ~ serify("virus effect on staying (", "<i>&epsilon;</i>", ")"),
                         x == "zeta" ~ serify("parasitoid responsiveness to aphid density (", "<i>&zeta;</i>", ")"),
                         x == "sd_N" ~ serify("initial aphid density SD (", "<i>&sigma;</i><sub>N</sub>", ")"),
                         x == "Y0" ~ serify("initial wasp density (", "<i>Y</i><sub>0</sub>", ")"),
                         x == "mean_N" ~ serify("initial aphid density mean (", "<i>&mu;</i><sub>N</sub>", ")"),
                         x == "N0" ~ serify("initial aphid density (", "<i>N</i><sub>0</sub>", ")"),
                         x == "alate_slope" ~ serify("slope for aphid density ~ alate offspring (", "<i>b</i><sub>slope</sub>", ")"),
                         x == "alate_max" ~ serify("max alate proportion (", "<i>b</i><sub>max</sub>", ")"),
                         x == "n_pseudo" ~ serify("number of *Pseudomonas* patches (", "<i>n</i><sub>P</sub>", ")"),
                         x == "K" ~ serify("aphid density dependence (", "<i>K</i>", ")"),
                         x == "spat_config" ~ "spatial configuration",
                         x == "wt_vp" ~ "virus placement",
                         x == "wt_pp" ~ "*Pseudomonas* spacing",
                         .default = x)
    }
    if (cap1) out <- first_cap(out)
    return(out)
}

# Descriptions for each y variable:
yvar_desc <- list(infect_time = "days to 5 plants infected",
                 infect_time_Inf = "percent where ≥ 5 plants were infected",
                 outbreak_size = "outbreak size",
                 log_outbreak_size = "log10(outbreak size)",
                 sd_outbreak_size = "outbreak size SD",
                 p_emerge = "prob. emergence",
                 p_outbreak = "outbreak probability",
                 p_alates = "mean alate proportion",
                 log_aphids = "mean log aphid abundance",
                 aphids = "mean aphid abundance",
                 log_alates = "mean log alate abundance",
                 alates = "mean alate abundance",
                 log_parasitized = "mean log parasitized aphid abundance",
                 parasitized = "mean parasitized aphid abundance",
                 log_mummies = "mean log mummy abundance",
                 mummies = "mean mummy abundance",
                 log_wasps = "mean log wasp abundance",
                 wasps = "mean wasp abundance")


spat_config_lvls <- c("random", "no virus", "diagonal",
                      "near virus", "far virus", "over virus")
spat_config_abbrevs <- c("random" = "rnd+v",
                         "no virus" = "rnd&minus;v",
                         "diagonal" = "diag",
                         "near virus" = "nr v",
                         "far virus" = "fr v",
                         "over virus" = "on v")


#'
#' I add this to ggplot objects when I want to use them inside a figure
#' I'm stitching together in Adobe Illustrator, where I'll add all the titles
#' and annotations.
#'
illustrator_theme <- theme(plot.title = element_blank(),
                           plot.tag = element_blank(),
                           legend.position = "none",
                           axis.title.y = element_blank(),
                           axis.title.y.right = element_blank(),
                           axis.title.y.left = element_blank(),
                           axis.title.x = element_blank(),
                           strip.text = element_blank(),
                           strip.text.x = element_blank(),
                           strip.text.y = element_blank(),
                           panel.background = element_rect(fill="transparent", color=NA),
                           plot.background = element_rect(fill="transparent", color=NA),
                           strip.background = element_blank())

# Nonparametric bootstrap CI:
boot_ci <- function(x, alpha = 0.01, R = 2000L) {
    b <- sapply(1:R, \(i) mean(sample(x, replace = TRUE), na.rm = TRUE))
    ci <- tibble(lo = quantile(b, alpha/2),
                 hi = quantile(b, 1-alpha/2))
    return(ci)
}



#' Run lil_landscape under two simulation scenarios: large or small,
#' and under two parameter combo types:
#' "low" (Pseudomonas decreases outbreak size)
#' "high" (Pseudomonas increases outbreak size)
#'
#' Use of `...` allows you to adjust parameter values.
#'
run_sim_combos <- function(type,
                           n_pseudo,
                           large_sims = FALSE,
                           return_args = FALSE,
                           ...) {

    # type = "low"; n_pseudo = 3L; large_sims = TRUE; return_args = FALSE

    stopifnot(length(type) == 1L && type %in% c("low", "high"))
    stopifnot(length(n_pseudo) == 1L && is.numeric(n_pseudo) && n_pseudo >= 0)
    stopifnot(length(large_sims) == 1L && is.logical(large_sims))

    shared_args <- list(sd_N = 0,
                        Y0 = 1,
                        pseudo_surv = 0.85,
                        n_pseudo = n_pseudo,
                        spat_config = "diagonal")

    if (large_sims) {
        size_args <- list_assign(shared_args, n_sims = 1000, summ = "all")
    } else {
        size_args <- list_assign(shared_args, n_sims = 1, summ = "none")
    }

    if (type == "low") {
        args <- list_assign(size_args,
                            N0 = 60,
                            zeta = 1)
    } else {
        args <- list_assign(size_args,
                            N0 = 45,
                            zeta = 0.3)
    }

    args <- list_assign(args, ...)

    if (return_args) return(args)

    return(do.call(lil_plantscape, args))
}


# # Used PlotDigitizer to get endpoints of regression lines of
# # minutes parasitoids foraged on plant ~ number of aphids * whether parasitoid
# # encountered an aphid from Ives et al. (1999)
# wasp_forage_lm <- "encounter, n_aphids, min_foraged
# yes, 3.625149424824349, 6.419469008502544
# yes, 54.704106691221476, 4.229949925046568
# no, 3.5784584360456937, 1.8839591285895756
# no, 54.463777356970255, 1.2101207112360788" |>
#     I() |>
#     read_csv(col_types = "fdd") |>
#     mutate(across(n_aphids:min_foraged, log10)) |>
#     arrange(encounter, n_aphids) |>
#     split(~ encounter) |>
#     map(\(x) {
#         lm(min_foraged ~ n_aphids, x)
#     })
#
#
# # From same paper, probability of parasitoid encountering an aphid based
# # on number of aphids on plant:
# pr_encounter <- function(x) {
#     p <- numeric(length(x))
#     p[x > 0] <- aeonia::inv_logit(-1.99 + 0.82 * log(x[x > 0]))
#     return(p)
# }
#
#
# # Put both above together for expected minutes for aphid abundances:
# wasp_forage_minutes <- function(n_aphids) {
#     pr_e <- pr_encounter(n_aphids)
#     pred_df <- data.frame(n_aphids = log10(n_aphids))
#     min_enc <- 10^predict(wasp_forage_lm$yes, pred_df) |> unname()
#     min_noenc <- 10^predict(wasp_forage_lm$no, pred_df) |> unname()
#     return(pr_e * min_enc + (1 - pr_e) * min_noenc)
# }
# # ratio of foraging time for with / without encounter
# foraging_ratio <- function(n_aphids) {
#     pred_df <- data.frame(n_aphids = log10(n_aphids))
#     min_enc <- 10^predict(wasp_forage_lm$yes, pred_df) |> unname()
#     min_noenc <- 10^predict(wasp_forage_lm$no, pred_df) |> unname()
#     return(min_enc / min_noenc)
# }

#
#
#
#
# require(gtools)
# require(tidyverse)
#
# p_load <- 0.5
#
# # Probability of at least one infected plant visited:
# p_inf_visit <- function(p_infected, plant_visits) {
#     1 - (1 - p_infected)^plant_visits
# }
# # Probability of at least one UNinfected plant visited:
# p_uninf_visit <- function(p_infected, plant_visits) {
#     1 - p_infected^plant_visits
# }
# # Expected number of events where uninfected visited after infected:
# expected_visits <- function(p_infected, plant_visits) {
#     perms <- permutations(n = 2, r = plant_visits + 1L, v = 0:1,
#                           repeats.allowed = TRUE)
#     n_events <- sapply(1:nrow(perms), \(i) sum(diff(perms[i,]) == -1))
#     out <- numeric(length(p_infected))
#     for (j in 1:length(p_infected)) {
#         p_inf <- p_infected[j]
#         probs <- sapply(1:nrow(perms), \(i) {
#             (1 - p_inf)^sum(perms[i,] == 0) * p_inf^sum(perms[i,] == 1)
#         })
#         out[j] <- sum(n_events * probs)
#     }
#     return(out)
# }
#
#
# tibble(plant_visits = 1:5) |>
#     mutate(n_visits = map(plant_visits, \(pv) {
#         ev <- expected_visits(1:100 / 100, pv)
#         return(tibble(p_infected = 1:100 / 100, n_visits = ev))
#     })) |>
#     unnest(n_visits) |>
#     ggplot(aes(p_infected, n_visits, color = factor(plant_visits))) +
#     geom_hline(yintercept = 0, color = "gray70") +
#     geom_line() +
#     scale_color_viridis_d(begin = 0.1, end = 0.9)















