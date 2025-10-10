
suppressPackageStartupMessages({
    library(tidyverse)
    library(aeonia)
    library(future.apply)
    library(progressr)
    library(patchwork)
    library(ggtext)
})

# Set threads for simulations:
options("mc.cores" = max(1L, parallel::detectCores()-2))

plan(multisession, workers = options()[["mc.cores"]])
handlers(global = TRUE)
handlers("progress")

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
CC <- function(K, B, m){
    L  <- rbind(c(pop_info$surv_j, pop_info$fecund),
                c(pop_info$recruit, pop_info$surv_a))
    K * (max(abs(eigen(L)[["values"]])) * (1 - B) * (1 - m) - 1)
}


spp_pal <- viridisLite::plasma(101)[c(10, 50, 80)] |>
    set_names(c("aphids", "alates", "wasps"))
# pseudo_pal <- c(`3` = "#1E90FF", `0` = "gray60")


pretty_params <- function(x) {
    case_when(x == "alpha" ~ "&alpha;",
              x == "beta" ~ "&beta;",
              x == "epsilon" ~ "&epsilon;",
              x == "wasp_disp_m0" ~ "&delta;<sub>0</sub>",
              x == "wasp_disp_m1" ~ "&delta;<sub>z</sub>",
              x == "mu_Y" ~ "&mu;<sub>Y</sub>",
              x == "mu_N" ~ "&mu;<sub>N</sub>",
              x == "n_pseudo" ~ "n<sub>P</sub>",
              .default = x)
}
