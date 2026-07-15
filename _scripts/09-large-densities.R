
#'
#' Density plots for larger landscape simulations.
#' Do NOT run this locally!
#'
#' Instead run on cluster using an interactive job:
#' cd /home2/lan68/09-large-densities
#' srun -N 1 -n 1 -c 50 --mem=200G --time=20:00:00 --job-name="lp-dens" --pty R --vanilla
#'
#' Then, when the jobs are done (assuming you're back on your desktop in
#' directory `~/GitHub/Cornell/aeonia/_scripts`):
#'
#' export RDS_DIR="/home2/lan68/09-large-densities"
#' scp lan68@cbsugreischar.biohpc.cornell.edu:${RDS_DIR}/density-sims.csv \
#'     ./interm-data/


source("../03-large-preamble.R")

# Factor for parasitoid wasp responsiveness to aphid densities:
wasp_resp_fct <- factor(1:2, labels = c("weak", "strong"))

# Levels of Pseudomonas densities for large landscapes:
n_pseudo_lvls <- as.integer(c(0, 0:4 * 2000 + 1000))



large_simmer <- function(landscape, wasp_resp, p_load = 0.5, sd_N = 0,
                         virus_attract = 1, pseudo_repel = 1, ...) {

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
                 n_sims = 100L,
                 summ = "all")

    args <- list_assign(args, ...)

    sims <- do.call(big_plantscape, args)

    return(sims)

}


# Takes ~26 sec
set.seed(314679353)
dens_sims <- crossing(wr = wasp_resp_fct,
                      np = n_pseudo_lvls) |>
    pmap(\(wr, np) {
        if (np > 0) {
            .landscape <- sim_df |>
                filter(n_pseudo == np,
                       # uniform, virus starts on uninhabited plant:
                       wt_pp == 1, wt_vp < 1) |>
                getElement("landscape") |>
                getElement(1)
        } else {
            .landscape <- array(c(1L, rep(0L, 100L*100L-1L)), rep(100L, 3))
        }
        summer <- function(x, p) x[is.na(p)]       # total across plants
        maxxer <- function(x, p) max(x[!is.na(p)]) # max per plant

        large_simmer(.landscape, wasp_resp = paste(wr),
                     summ = "none", out_stages = "two",
                     n_sims = dim(.landscape)[3]) |>
            mutate(plant = interaction(y, x),
                   aphids = aphids_juv + aphids_adu + alates_juv + alates_adu + parasitized,
                   alates = alates_adu) |>
            select(rep, plant, time, aphids, alates, wasps, virus) |>
            group_by(time, rep) |>
            summarize(across(aphids:virus, list(sum = \(x) summer(x, p = plant),
                                                max = \(x) maxxer(x, p = plant))),
                      .groups = "drop_last") |>
            summarize(across(aphids_sum:virus_max, mean), .groups = "drop") |>
            mutate(wasp_resp = wr, n_pseudo = np) |>
            select(wasp_resp, n_pseudo, everything())
    }, .progress = list(clear = FALSE,
                        format = paste("{cli::pb_bar}",
                                       "{cli::pb_percent}",
                                       "[{cli::pb_elapsed}] |",
                                       "ETA: {cli::pb_eta}"))) |>
    list_rbind()


write_csv(dens_sims, "density-sims.csv")


