
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
rds_files <- list(extreme_large = "extremes-large-sims.rds",
                  extreme_manip = "extremes-manip-sims.rds",
                  extreme_manip2 = "extremes-manip2-sims.rds") |>
    map(\(x) paste0("_scripts/interm-data/", x))



# Set threads for simulations:
options("mc.cores" = max(1L, parallel::detectCores()-2))

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

# Carrying capacity with no alates or natural enemies:
CC <- function(K, pseudo_surv, pred_surv){
    L  <- rbind(c(pop_info$surv_j, pop_info$fecund),
                c(pop_info$recruit, pop_info$surv_a))
    K * (max(abs(eigen(L)[["values"]])) * pseudo_surv * pred_surv - 1)
}


color_pal <- c("#999999", "#0046D2", "#FF64FF",
               viridisLite::inferno(100)[c(50, 25, 80)]) |>
    set_names(c("no_pseudo", "pseudo", "virus", "aphids", "alates", "wasps"))
# scales::show_col(color_pal, labels = FALSE)

# for numbers of pseudomonas patches:
np_pal <- color_pal[c("no_pseudo", "pseudo")] |>
    set_names(c("0", "3"))
# scales::show_col(np_pal, labels = FALSE)

np_shapes <- c(`0` = 1, `3` = 19)
np_linetypes <- c(`0` = "22", `3` = "solid")


serify <- function(prefix, x, suffix) {
    xx <- paste0("<span style=\"font-family: serif\">", x, "</span>")
    paste0(prefix, xx, suffix)
}

# Only first letter capitalized (others not forced to lowercase, as in `str_to_sentence`):
first_cap <- \(str) paste(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)), sep="")




pretty_params <- function(x, short = FALSE, cap1 = FALSE, serif = FALSE) {
    if (short) {
        out <- case_when(x == "pseudo_surv" ~ "&psi;",
                         x == "virus_attract" ~ "&nu;",
                         x == "pseudo_repel" ~ "&rho;",
                         x == "epsilon" ~ "&epsilon;",
                         x == "zeta" ~ "&zeta;",
                         x == "sd_N" ~ "&sigma;<sub>N</sub>",
                         x == "Y0" ~ "Y<sub>0</sub>",
                         x == "mean_N" ~ "&mu;<sub>N</sub>",
                         x == "N0" ~ "N<sub>0</sub>",
                         x == "alate_slope" ~ "b<sub>slope</sub>",
                         x == "alate_max" ~ "b<sub>max</sub>",
                         x == "n_pseudo" ~ "n<sub>P</sub>",
                         x == "spat_config" ~ "spat. config.",
                         x == "wt_vp" ~ "*Pseudo.* place.",
                         x == "wt_pp" ~ "*Pseudo.* spacing",
                         .default = x)
        if (serif && any(x != "spat_config")) {
            out[x != "spat_config"] <- serify("", out[x != "spat_config"], "")
        }
    } else {
        out <- case_when(x == "pseudo_surv" ~ serify("*Pseudomonas* survival (", "&psi;", ")"),
                         x == "virus_attract" ~ serify("virus attraction (", "&nu;", ")"),
                         x == "pseudo_repel" ~ serify("*Pseudomonas* repellence (", "&rho;", ")"),
                         x == "epsilon" ~ serify("virus effect on staying (", "&epsilon;", ")"),
                         x == "zeta" ~ serify("wasp density response (", "&zeta;", ")"),
                         x == "sd_N" ~ serify("initial aphid density SD (", "&sigma;<sub>N</sub>", ")"),
                         x == "Y0" ~ serify("initial wasp density (", "Y<sub>0</sub>", ")"),
                         x == "mean_N" ~ serify("initial aphid density mean (", "&mu;<sub>N</sub>", ")"),
                         x == "N0" ~ serify("initial aphid density (", "N<sub>0</sub>", ")"),
                         x == "alate_slope" ~ serify("slope for aphid density ~ alate offspring (", "b<sub>slope</sub>", ")"),
                         x == "alate_max" ~ serify("max alate proportion (", "b<sub>max</sub>", ")"),
                         x == "n_pseudo" ~ serify("number of *Pseudomonas* patches (", "n<sub>P</sub>", ")"),
                         x == "K" ~ serify("aphid density dependence (", "K", ")"),
                         x == "spat_config" ~ "spatial configuration",
                         x == "wt_vp" ~ "*Pseudomonas* placement",
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

    stopifnot(length(type) == 1L && type %in% c("low", "high"))
    stopifnot(length(n_pseudo) == 1L && is.numeric(n_pseudo) && n_pseudo >= 0)
    stopifnot(length(large_sims) == 1L && is.logical(large_sims))

    shared_args <- list(sd_N = 0,
                        K = 12.5e3,
                        virus_attract = 1,
                        pseudo_repel = 1,
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
                            Y0 = 3,
                            N0 = 110,
                            zeta = 1)
    } else {
        args <- list_assign(size_args,
                            Y0 = 2,
                            N0 = 10,
                            zeta = 0.37)
    }

    args <- list_assign(args, ...)

    if (return_args) return(args)

    return(do.call(lil_plantscape, args))
}
