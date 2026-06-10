
#'
#' Empirical zeta estimates for larger landscape simulations.
#' Do NOT run this locally!
#'
#' Note: Because densities vary little across simulations, I'm only running
#' 12 simulations.
#'
#' Instead run on cluster using an interactive job:
#' cd /home2/lan68/11-large-emp-zeta
#' srun -N 1 -n 1 -c 14 --mem=200G --time=20:00:00 --job-name="lp-zeta" --pty R --vanilla
#'
#' Then, when the jobs are done (assuming you're back on your desktop in
#' directory `~/GitHub/Cornell/aeonia/_scripts`):
#'
#' export RDS_DIR="/home2/lan68/11-large-emp-zeta"
#' scp lan68@cbsugreischar.biohpc.cornell.edu:${RDS_DIR}/large-zeta-sims.csv.gz \
#'     ./interm-data/


source("../03-large-preamble.R")

# Levels of Pseudomonas densities for large landscapes (not including no Pseudo.):
n_pseudo_lvls <- as.integer(0:4 * 2000 + 1000)


large_simmer <- function(landscape,
                         zeta,
                         p_load = 0.5,
                         ...) {

    args <- list(landscape = landscape,
                 sd_N = 0,
                 virus_attract = 1,
                 pseudo_repel = 1,
                 Y0 = 250,
                 N0 = 55,
                 zeta = zeta,
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = 12L,
                 summ = "none")

    args <- list_assign(args, ...)

    sims <- do.call(big_plantscape, args)

    return(sims)

}


one_emp_zeta_sim <- function(np, zeta) {

    if (np > 0) {
        .landscape <- sim_df |>
            filter(n_pseudo == np,
                   # uniform, virus starts on uninhabited plant:
                   wt_pp == 1, wt_vp < 1) |>
            getElement("landscape") |>
            getElement(1) |>
            base::`[`(,,1:12)
    } else {
        .landscape <- array(c(1L, rep(0L, 100L*100L-1L)), c(100L, 100L, 12L))
    }

    sims <- large_simmer(landscape = .landscape, zeta = zeta) |>
        filter(!is.na(x)) |>
        mutate(plant = interaction(y, x),
               aphids = aphids + parasitized) |>
        select(rep, plant, time, virus, aphids, alates, wasps) |>
        split(~ rep) |>
        map(\(x) {
            # Filter for time at which max aphid density occurs.
            # If I don't do this now, resulting data frame will be too large
            max_N_t <- x |>
                group_by(time) |>
                summarize(aphids = max(aphids + alates)) |>
                filter(aphids == max(aphids)) |>
                getElement("time") |> getElement(1)
            x <- x |> filter(time == max_N_t)
            # Now calculate predicted wasps:
            Y <- sum(x$wasps)
            p_plant_names <- which(.landscape[,,x$rep[[1]]] > 1L, arr.ind = TRUE) |>
                # Have to use rev() bc of how x and y dims are arranged
                apply(1, \(x) paste(rev(x), collapse = "."))
            stopifnot(all(p_plant_names %in% levels(x$plant)))
            pseudo <- x$plant %in% p_plant_names
            pred_wasps <- numeric(nrow(x))
            pred_wasps[pseudo] <- 1
            pred_wasps[!pseudo] <- 3.76
            pred_wasps <- Y * pred_wasps / sum(pred_wasps)
            x[["pred_wasps"]] <- pred_wasps
            x[["pseudo"]] <- pseudo
            return(x)
        }) |>
        list_rbind() |>
        mutate(zeta = .env$zeta, n_pseudo = np) |>
        select(zeta, n_pseudo, everything())

    return(sims)

}



# Takes ~10 min
set.seed(1974555786)
emp_zeta_sims <- crossing(np = n_pseudo_lvls,
                          zeta = 2:9 / 10) |>
    pmap(one_emp_zeta_sim,
         .progress = list(clear = FALSE,
                        format = paste("{cli::pb_bar}",
                                       "{cli::pb_percent}",
                                       "[{cli::pb_elapsed}] |",
                                       "ETA: {cli::pb_eta}"))) |>
    list_rbind()


# Compressing this one bc it's pretty large (>300 MB) uncompressed
write_csv(emp_zeta_sims, "large-zeta-sims.csv.gz")
