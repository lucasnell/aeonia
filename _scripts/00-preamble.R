
suppressPackageStartupMessages({
    library(pillar) # num() for setting sig figs in tibble printing
    library(tidyverse)
    library(aeonia)
    library(patchwork)
    library(ggtext)
    library(scico)
})

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

serify <- function(prefix, x, suffix) {
    xx <- paste0("<span style=\"font-family: serif\">", x, "</span>")
    paste0(prefix, xx, suffix)
}

# Only first letter capitalized (others not forced to lowercase, as in `str_to_sentence`):
first_cap <- \(str) paste(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)), sep="")




pretty_params <- function(x, short = FALSE, cap1 = FALSE) {
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
                         .default = x)
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



#'
#' I add this to ggplot objects when I want to use them inside a figure
#' I'm stitching together in Adobe Illustrator, where I'll add all the titles
#' and annotations.
#'
illustrator_theme <- theme(plot.title = element_blank(),
                           legend.position = "none",
                           axis.title.y = element_blank(),
                           axis.title.x = element_blank(),
                           strip.text = element_blank(),
                           panel.background = element_rect(fill="transparent", color=NA),
                           plot.background = element_rect(fill="transparent", color=NA),
                           strip.background = element_blank())


#' Run lil_landscape under two simulation scenarios: large or small,
#' and under two parameter combo types:
#' "low" (Pseudomonas decreases outbreak size)
#' "high" (Pseudomonas increases outbreak size)
#'
run_sim_combos <- function(type, n_pseudo, large_sims = FALSE, ...) {

    stopifnot(length(type) == 1L && type %in% c("low", "high"))
    stopifnot(length(n_pseudo) == 1L && is.numeric(n_pseudo) && n_pseudo >= 0)
    stopifnot(length(large_sims) == 1L && is.logical(large_sims))

    shared_args <- list(Y0 = 2,
                        sd_N = 0,
                        K = 12.5e3,
                        pseudo_surv = 0.9,
                        virus_attract = 1,
                        n_pseudo = n_pseudo,
                        spat_config = "diagonal")

    if (large_sims) {
        size_args <- list_assign(shared_args, n_sims = 1000, summ = "all")
    } else {
        size_args <- list_assign(shared_args, n_sims = 1, summ = "none")
    }

    if (type == "low") {
        args <- list_assign(size_args,
                            N0 = 100,
                            pseudo_repel = 10,
                            zeta = 1)
    } else {
        args <- list_assign(size_args,
                            N0 = 10,
                            pseudo_repel = 1.5,
                            zeta = 0.2)
    }

    args <- list_assign(args, ...)

    return(do.call(lil_plantscape, args))
}
