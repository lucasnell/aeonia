
#'
#' Dispersal patterns with virus attraction
#' Do NOT run this locally!
#'
#' Instead run on cluster using an interactive job:
#' cd /home2/lan68/13-virus_attract
#' srun -N 1 -n 1 -c 50 --mem=400G --time=20:00:00 --job-name="virus-attract" --pty R --vanilla
#'
#' Then, when the jobs are done (assuming you're back on your desktop in
#' directory `~/GitHub/Cornell/aeonia/_scripts`):
#'
#' export RDS_DIR="/home2/lan68/13-virus-attract"
#' scp lan68@cbsugreischar.biohpc.cornell.edu:${RDS_DIR}/virus-attract-disps.rds \
#'     ./interm-data/


# source("../03-large-preamble.R")
source("_scripts/03-large-preamble.R")

# Factor for parasitoid wasp responsiveness to aphid densities:
wasp_resp_fct <- factor(1:2, labels = c("weak", "strong"))

# Levels of Pseudomonas densities for large landscapes:
n_pseudo_lvls <- as.integer(c(0, 0:4 * 2000 + 1000))



# large landscape simmer, just focusing on dispersals through time and space
large_dispersal_simmer <- function(n_pseudo,
                                   wasp_resp,
                                   p_load,
                                   fly_p,
                                   virus_attract,
                                   n_sims = 1L) {

    # n_pseudo = 3000L; wasp_resp = "weak"; p_load = 1; fly_p = 0.1
    # virus_attract = 1; n_sims = 1L
    # rm(n_pseudo, wasp_resp, p_load, fly_p, virus_attract, n_sims)
    # rm(landscape, args, sims)

    if (n_pseudo > 0) {
        landscape <- sim_df |>
            filter(n_pseudo == .env$n_pseudo,
                   # uniform, virus starts on uninhabited plant:
                   wt_pp == 1, wt_vp < 1) |>
            getElement("landscape") |>
            getElement(1)
        landscape <- landscape[,,1:n_sims,drop=FALSE]
    } else landscape <- array(c(1L, rep(0L, 99999L)), c(100L, 100L, n_sims))

    args <- list(landscape = landscape,
                 virus_attract = virus_attract,
                 fly_p = fly_p,
                 zeta = ifelse(wasp_resp == "weak", 0.1, 0.9),
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 sd_N = 0,
                 pseudo_repel = 1,
                 Y0 = 250,
                 N0 = 55,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = dim(landscape)[3],
                 summ = "time",
                 out_dispersals = "in",
                 force_disps = TRUE)

    return(do.call(big_plantscape, args))


}


# Takes ~ 15 sec
sims <- large_simmer(n_pseudo = 7000L, wasp_resp = "weak", p_load = 1, virus_attract = 1)
sims$n_infected |> mean()



# with n_pseudo = 7000L, wasp_resp = "weak", p_load = 1, virus_attract = 100
# 8609.44






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
                     summ = "none", out_stages = TRUE,
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


