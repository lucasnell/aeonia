
suppressPackageStartupMessages({
    library(tidyverse)
    library(aeonia)
    library(patchwork)
    library(ggtext)
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


spp_pal <- viridisLite::plasma(101)[c(75, 40, 10)] |>
    set_names(c("aphids", "alates", "wasps"))
# Uncomment below to see (in RStudio):
# "#F79143FF" "#AD2793FF" "#3E049CFF"


# viridisLite::plasma(101)[c(75, 60, 30, 10)]
# "#F79143FF" "#DE6263FF" "#8D0BA5FF" "#3E049CFF"


# for numbers of pseudomonas patches:
np_pal <- c(`0` = "goldenrod", `3` = "dodgerblue")



# pseudo_pal <- c(`3` = "#1E90FF", `0` = "gray60")

serify <- function(prefix, x, suffix) {
    xx <- paste0("<span style=\"font-family: serif\">", x, "</span>")
    paste0(prefix, xx, suffix)
}

# Only first letter capitalized (others not forced to lowercase, as in `str_to_sentence`):
first_cap <- \(str) paste(toupper(substr(str, 1, 1)), substr(str, 2, nchar(str)), sep="")


pretty_params <- function(x, short = FALSE) {
    if (short) {
        case_when(x == "pseudo_surv" ~ "&psi;",
                  x == "virus_attract" ~ "&nu;",
                  x == "pseudo_repel" ~ "&rho;",
                  x == "epsilon" ~ "&epsilon;",
                  x == "zeta" ~ "&zeta;",
                  x == "sd_N" ~ "&sigma;<sub>N</sub>",
                  x == "Y0" ~ "Y<sub>0</sub>",
                  x == "mean_N" ~ "&mu;<sub>N</sub>",
                  x == "alate_slope" ~ "b<sub>slope</sub>",
                  x == "alate_max" ~ "b<sub>max</sub>",
                  x == "n_pseudo" ~ "n<sub>P</sub>",
                  x == "spat_config" ~ "spat. conf.",
                  .default = x)
    } else {
        case_when(x == "pseudo_surv" ~ serify("*Pseudomonas* survival (", "&psi;", ")"),
                  x == "virus_attract" ~ serify("virus attraction (", "&nu;", ")"),
                  x == "pseudo_repel" ~ serify("*Pseudomonas* repellence (", "&rho;", ")"),
                  x == "epsilon" ~ serify("virus effect on staying (", "&epsilon;", ")"),
                  x == "zeta" ~ serify("wasp density response (", "&zeta;", ")"),
                  x == "sd_N" ~ serify("initial aphid density SD (", "&sigma;<sub>N</sub>", ")"),
                  x == "Y0" ~ serify("initial wasp density (", "Y<sub>0</sub>", ")"),
                  x == "mean_N" ~ serify("initial aphid density mean (", "&mu;<sub>N</sub>", ")"),
                  x == "alate_slope" ~ serify("slope for aphid density ~ alate offspring (", "b<sub>slope</sub>", ")"),
                  x == "alate_max" ~ serify("max alate proportion (", "b<sub>max</sub>", ")"),
                  x == "n_pseudo" ~ serify("number of *Pseudomonas* patches (", "n<sub>P</sub>", ")"),
                  x == "K" ~ serify("aphid density dependence (", "K", ")"),
                  x == "spat_config" ~ "spatial configuration",
                  .default = x)
    }
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
