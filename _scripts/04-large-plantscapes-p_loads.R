
#'
#' cd /home2/lan68
#' srun -N 1 -n 1 -c 50 --mem=25G --time=20:00:00 --job-name="large-test" --pty bash -l
#' R --vanilla
#'

source("03-large-preamble.R")


large_simmer <- function(landscape, sd_N, virus_attract, pseudo_repel,
                         Y0, N0, zeta, p_load) {

    args <- list(landscape = landscape,
                 sd_N = sd_N,
                 virus_attract = virus_attract,
                 pseudo_repel = pseudo_repel,
                 Y0 = Y0,
                 N0 = N0,
                 zeta = zeta,
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = 100L,
                 summ = "all")

    return(do.call(big_plantscape, args))

}




# =============================================================================*
# First test for starting conditions and zeta ----
# =============================================================================*

landscape1 <- sim_df |>
    filter(n_pseudo == 5000, wt_vp == 1e-6, wt_pp == 1,
           # These do not affect landscape:
           type == "low", sd_N == 0, virus_attract == 1,
           pseudo_repel == 1) |>
    getElement("landscape") |> getElement(1)
landscape0 <- array(c(1L, rep(0L, 99999L)), c(100L, 100L, 1L))



one_test <- function(p_load, Y0, N0, zeta) {

    n_inf1 <- large_simmer(landscape = landscape1,
                          sd_N = 0, virus_attract = 1, pseudo_repel = 1,
                          Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
        getElement("outbreak_size")
    n_inf0 <- large_simmer(landscape = landscape0,
                          sd_N = 0, virus_attract = 1, pseudo_repel = 1,
                          Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
        getElement("outbreak_size")

    tibble(n_pseudo = c(5000L, 0L),
           p_load = .env$p_load, Y0 = .env$Y0, N0 = .env$N0, zeta = .env$zeta,
           outbreak_size = c(mean(n_inf1[n_inf1 > 1]), mean(n_inf0[n_inf0 > 1])),
           p_emerge = c(mean(n_inf1 > 1), mean(n_inf0 > 1)))

}

# Takes ~ min
test_sims <- crossing(p_load = round(seq(0.05, 0.2, 0.05), 2),
                      Y0 = c(90, 100, 110),
                      N0 = c(9:11, 100, 110, 120),
                      zeta = c(0.08, 0.1, 0.12, 0.9, 0.95, 1)) |>
    # Filter for two scenarios: Pseudomonas promotes or inhibits viruses, resp.:
    filter((zeta < 0.2 & N0 < 20) | (zeta > 0.8 & N0 > 100)) |>
    pmap(one_test, .progress = list(clear = FALSE,
                                    format = paste("{cli::pb_bar}",
                                                   "{cli::pb_percent}",
                                                   "[{cli::pb_elapsed}] |",
                                                   "ETA: {cli::pb_eta}"))) |>
    list_rbind()

write_rds(test_sims, "TEST-large-plantscapes.rds", compress = "gz")

# Then run
# scp lan68@cbsugreischar.biohpc.cornell.edu:/home2/lan68/TEST-large-plantscapes.rds ~/GitHub/Cornell/aeonia/_scripts/interm-data/


# test_sims <- read_rds("_scripts/interm-data/TEST-large-plantscapes.rds")
#
#
# test_sims |>
#     mutate(outbreak_size = log10(outbreak_size)) |>
#     pivot_longer(outbreak_size:p_emerge, names_to = "outcome") |>
#     filter(!is.na(value)) |>
#     mutate(outcome = factor(outcome,
#                             levels = c("outbreak_size", "p_emerge"),
#                             labels = c("log10(outbreak size)",
#                                        "prob. emerge"))) |>
#     ggplot(aes(p_load, value, color = scenario)) +
#     geom_hline(yintercept = 0) +
#     geom_point() +
#     geom_line() +
#     facet_wrap(~ outcome, ncol = 1, scales = "free")
#

