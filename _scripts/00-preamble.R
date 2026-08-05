

suppressPackageStartupMessages({
    library(pillar) # num() for setting sig figs in tibble printing
    library(tidyverse)
    library(aeonia)
    library(patchwork)
    library(ggtext)
    library(scico)
    library(viridisLite)
    library(RColorBrewer)
    library(colorspace)
})

# RDS / CSV files with simulation output:
interm_files <- list(extreme_manip = "extremes-manip-sims.rds",
                     extreme_manip2 = "extremes-manip2-sims.rds",
                     dens_sims = "density-sims.csv",
                     dd_disp_sims = "dd_disp-sims.csv",
                     large_zeta_sims = "large-zeta-sims.csv.gz") |>
    map(\(x) paste0("_scripts/interm-data/", x))



# Set threads for simulations:
options("mc.cores" = max(1L, parallel::detectCores()-2L))
# And for readr functions
options("readr.num_threads" = max(1L, parallel::detectCores()-2L))


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




# For numbers of pseudomonas plants:
# First, for number of Pseudomonas plants, from 0 to 9000, by increments of 1000
full_np_pal <- viridisLite::plasma(10, begin = 0.1, end = 0.9, direction = -1L) |>
    set_names(paste(0:9 * 1000L))
# scales::show_col(full_np_pal, labels = FALSE)
# scales::show_col(full_np_pal[c(1,2,4,6,8,10)], labels = FALSE)

np_pal <- full_np_pal[c(1, 7)] |>
    set_names(c("0", "3"))
# scales::show_col(np_pal, labels = FALSE)



# For parameters:
#
# par_pal <- viridisLite::turbo(100)[c(10, 30, 90, 70)] |>
#     set_names("wt_vp", "pseudo_repel", "virus_attract", "sd_N")
# # scales::show_col(par_pal, labels = FALSE)
#
# To get the colors below (that play well in cmyk or rgb color models),
# I inserted the colors from above into Illustrator as new swatches,
# then opened the (modified) swatches back up and pasted the hex codes below.
par_pal <- c(wt_vp = "#4c5aa8", pseudo_repel = "#48c2c6",
             virus_attract = "#c22c26", sd_N = "#faa733")




# Factor for parasitoid wasp responsiveness to aphid densities:
wasp_resp_fct <- factor(1:2, labels = c("weak", "strong"))

# Levels of Pseudomonas densities for large landscapes:
n_pseudo_lvls <- as.integer(c(0, 0:4 * 2000 + 1000))


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
                         x == "N0" ~ "<i>N</i><sub>0</sub>",
                         x == "n_pseudo" ~ "<i>n</i><sub>P</sub>",
                         x == "wt_vp" ~ "virus place.",
                         x == "wt_pp" ~ "<i>&tau;</i>",
                         .default = x)
        if (serif) {
            out <- serify("", out, "")
        }
    } else {
        out <- case_when(x == "pseudo_surv" ~ serify("aphid *Pseudomonas* survival (", "<i>&psi;</i>", ")"),
                         x == "virus_attract" ~ serify("virus attraction (", "<i>&nu;</i>", ")"),
                         x == "pseudo_repel" ~ serify("*Pseudomonas* repellence (", "<i>&rho;</i>", ")"),
                         x == "epsilon" ~ serify("virus effect on staying (", "<i>&epsilon;</i>", ")"),
                         x == "zeta" ~ serify("parasitoid responsiveness to aphid density (", "<i>&zeta;</i>", ")"),
                         x == "sd_N" ~ serify("starting aphid density SD (", "<i>&sigma;</i><sub>N</sub>", ")"),
                         x == "Y0" ~ serify("starting parasitoid density (", "<i>Y</i><sub>0</sub>", ")"),
                         x == "N0" ~ serify("starting aphid density (", "<i>N</i><sub>0</sub>", ")"),
                         x == "n_pseudo" ~ serify("number of *Pseudomonas* plants (", "<i>n</i><sub>P</sub>", ")"),
                         x == "K" ~ serify("aphid density dependence (", "<i>K</i>", ")"),
                         x == "wt_vp" ~ "virus placement",
                         x == "wt_pp" ~ "*Pseudomonas* spacing",
                         .default = x)
    }
    if (cap1) out <- first_cap(out)
    return(out)
}

# yvar_desc ----
# Descriptions for each y variable:
yvar_desc <- list(infect_time = "days to 5 plants infected",
                 infect_time_Inf = "percent where ≥ 5 plants were infected",
                 outbreak_size = "outbreak size",
                 log_outbreak_size = "log10(outbreak size)",
                 sd_outbreak_size = "outbreak size SD",
                 p_emerge = "prob. emergence",
                 p_outbreak = "outbreak probability",
                 n_infected = "peak infected plants",
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
                           plot.tag = element_blank(),
                           legend.position = "none",
                           axis.title.y = element_blank(),
                           axis.title.y.right = element_blank(),
                           axis.title.y.left = element_blank(),
                           axis.title.x = element_blank(),
                           axis.title.x.top = element_blank(),
                           axis.title.x.bottom = element_blank(),
                           strip.text = element_blank(),
                           strip.text.x = element_blank(),
                           strip.text.y = element_blank(),
                           panel.background = element_rect(fill="transparent", color=NA),
                           plot.background = element_rect(fill="transparent", color=NA),
                           strip.background = element_blank(),
                           plot.margin = margin(0,0,0,0))

# Nonparametric bootstrap CI:
ci_booter <- function(n_infected, .outcomes) {
    if (length(.outcomes) == 1L && .outcomes == "all")
        .outcomes <- c("p_emerge", "outbreak_size", "n_infected")
    out <- list()
    if ("p_emerge" %in% .outcomes) {
        out[["p_emerge"]] <- booter(as.integer(n_infected > 1))
    }
    if ("outbreak_size" %in% .outcomes) {
        obs <- n_infected[n_infected > 1]
        if (length(obs) < 2L) {
            if (length(obs) == 0L) obs <- NA_real_
            out[["outbreak_size"]] <- c(Lower = obs, Median = obs, Upper = obs)
        } else out[["outbreak_size"]] <- booter(obs)
    }
    if ("n_infected" %in% .outcomes) {
        out[["n_infected"]] <- booter(n_infected)
    }
    return(out)
}








#' Run lil_landscape under two simulation scenarios: large or small,
#' and under two levels of wasp responsiveness to aphid densities:
#' "strong" (so Pseudomonas decreases outbreak size)
#' "weak" (so Pseudomonas increases outbreak size)
#'
#' Use of `...` allows you to adjust parameter values.
#'
run_sim_combos <- function(wasp_resp,
                           n_pseudo,
                           large_sims = FALSE,
                           return_args = FALSE,
                           ...) {

    stopifnot(length(wasp_resp) == 1L && wasp_resp %in% c("strong", "weak"))
    stopifnot(length(n_pseudo) == 1L && is.numeric(n_pseudo) && n_pseudo >= 0)
    stopifnot(length(large_sims) == 1L && is.logical(large_sims))

    args <- list(sd_N = 0,
                 Y0 = 1,
                 N0 = 55,
                 pseudo_surv = 0.85,
                 n_pseudo = n_pseudo,
                 spat_config = "diagonal")

    if (large_sims) {
        args[["n_sims"]]  <- 1000
        args[["summ"]]  <- "all"
    } else {
        args[["n_sims"]]  <- 1
        args[["summ"]]  <- "none"
    }

    args[["zeta"]] <- ifelse(wasp_resp == "strong", 0.9, 0.1)

    args <- list_assign(args, ...)

    if (return_args) return(args)

    return(do.call(lil_plantscape, args))
}




scenario_title <- function(.wasp_resp,
                           .value = TRUE,
                           .break = TRUE,
                           .cap1 = TRUE,
                           .nl = "<br>") {
    .fmt <- ifelse(.break, "%s parasitoid<br>response to aphids",
                   "%s parasitoid response to aphids")
    if (.break && .nl != "<br>") .fmt <- .fmt |>
            str_replace_all("<br>", .nl)
    out <- sprintf(.fmt, .wasp_resp)
    if (.cap1) out <- first_cap(out)
    if (.value) {
        val <- sprintf("<i>&zeta;</i> = %.1f",
                       ifelse(.wasp_resp == "weak", 0.1, 0.9))
        out <- paste(out, serify("(", val, ")"), sep = ifelse(.break, .nl, " "))
    }
    return(out)
}




large_simmer <- function(landscape, wasp_resp, p_load, sd_N,
                         virus_attract, pseudo_repel, ...) {

    args <- list(landscape = landscape,
                 sd_N = sd_N,
                 virus_attract = virus_attract,
                 pseudo_repel = pseudo_repel,
                 Y0 = 250,
                 N0 = 55,
                 zeta = ifelse(wasp_resp == "weak", 0.1, 0.9),
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = dim(landscape)[3],
                 summ = "all")

    args <- list_assign(args, ...)

    return(do.call(big_plantscape, args))

}
