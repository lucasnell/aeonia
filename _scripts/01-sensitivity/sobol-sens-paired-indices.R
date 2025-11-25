#'
#' Small-scale sensitivity via Sobol indices
#'
#'

source("_scripts/01-sensitivity/00-preamble.R")
source("_scripts/01-sensitivity/sobol-preamble.R")


out_files <- list(summs = "sobol-sims-summs.rds",
                  diff_summs = "sobol-sims-diff-summs.rds") |>
    map(\(x) paste("_scripts/interm-data/", x))



if (!file.exists(out_files$summs)) {


    proc_one_file <- function(j) {
        # j = 1
        # rm(j, sims_j, summs_j)

        # Output from `sobol-sens-paired.R`:
        sims_j <- sprintf("~/_globus/sobol-sims-paired-%i.rds", j) |>
            read_rds()

        # Start for rows in this job's simulations:
        curr_start <- (j - 1L) * (obj_env$n_rows %/% obj_env$n_files) + 1L

        summs_j <- imap(sims_j, \(sim_set, i) {
            out_df <- sim_set |>
                mutate(combo = factor(curr_start + i - 1L,
                                      levels = 1:(obj_env$n_rows)),
                       sims = map(sims, \(s) s[1, names(obj_env$vary_pars)])) |>
                unnest(sims) |>
                select(combo, everything())
            for (yv in c(obj_env$yvars, "p_outbreak", "outbreak_size2", "sd_outbreak_size")) {
                out_df[[yv]] <- 0.0
            }
            for (j in 1:nrow(sim_set)) {
                for (yv in obj_env$yvars) {
                    y <- sim_set$sims[[j]][[yv]]
                    if (yv != "infect_time" && any(is.na(y))) stop(y, " has NA values")
                    out_df[[yv]][[j]] <- mean(y, na.rm = TRUE)
                }
                y <- sim_set$sims[[j]][["outbreak_size"]]
                # now do prob. outbreak happened:
                out_df[["p_outbreak"]][[j]] <- mean(y > 1)
                # and outbreak size when there was one:
                out_df[["outbreak_size2"]][[j]] <- mean(y[y > 1])
                # Lastly, SD(outbreak size):
                out_df[["sd_outbreak_size"]][[j]] <- sd(y)
            }
            return(out_df)

        }) |>
            list_rbind()

        return(summs_j)

    }

    proc_all_files <- function() {

        suppressPackageStartupMessages(library(future.apply))
        with(plan(multisession, workers = options()[["mc.cores"]], gc = TRUE),
             local = TRUE)
        oopts <- options(future.globals.maxSize = 1.5e9)  ## 1.5 GB
        on.exit(options(oopts))

        future_lapply(1:(obj_env$n_files),
                      proc_one_file,
                      future.globals = c("obj_env", "proc_one_file"),
                      future.packages = "tidyverse") |>
            list_rbind()
    }

    obj_env <- list()
    obj_env$n_rows <- nrow(sobol_mat)
    obj_env$n_files <- length(list.files("~/_globus", "sobol-sims-paired-.?.rds"))
    obj_env$vary_pars <- vary_pars
    obj_env$yvars <- yvars

    # Summarize each set of simulations:
    # Takes ~9 min w/ 6 threads
    sobol_summs <- proc_all_files()

    write_rds(sobol_summs, out_files$summs, compress = "gz")

    rm(proc_one_file, obj_env, proc_all_files)

} else {

    sobol_summs <- read_rds(out_files$summs)

}



if (!file.exists(out_files$diff_summs)) {

    # Same thing but looking at differences between with and without Pseudo:
    # Takes a few sec
    diff_sobol_summs <- sobol_summs |>
        group_by(combo, alate_dens, across(all_of(names(vary_pars)))) |>
        summarize(across(all_of(c(yvars, "p_outbreak", "outbreak_size2",
                                  "sd_outbreak_size")),
                         \(x) x[n_pseudo > 0] - x[n_pseudo == 0]),
                  .groups = "drop") |>
        mutate(across(starts_with("outbreak_size"), \(x) round(x, 2)))

    write_rds(diff_sobol_summs, out_files$diff_summs, compress = "gz")

} else {

    diff_sobol_summs <- read_rds(out_files$diff_summs)
}





set.seed(1998643658)
sobol_inds_np3 <- sobol_summs |>
    filter(alate_dens == 1) |>
    (\(ss) {
        Y <- ss[["outbreak_size"]][ss$n_pseudo > 0]
        ind <- sobol_indices(Y = Y, N = N, params = pretty_params(names(vary_pars)),
                             boot = TRUE, R = 2000)
        return(ind)
    })()

set.seed(512925036)
sobol_inds_diff <- diff_sobol_summs |>
    filter(alate_dens == 1) |>
    (\(ss) {
        Y <- ss[["outbreak_size"]]
        ind <- sobol_indices(Y = Y, N = N, params = pretty_params(names(vary_pars)),
                             boot = TRUE, R = 2000)
        return(ind)
    })()


sobol_inds_p1 <- plot(sobol_inds_np3) +
    labs(y = "Sobol' index (outbreak size with<br>three *Pseudomonas*)") +
    theme(axis.text.y = element_markdown(),
          axis.title.x = element_markdown()) +
    coord_flip()
sobol_inds_p2 <- plot(sobol_inds_diff) +
    labs(y = "Sobol' index (effect of *Pseudomonas*<br>on outbreak size)") +
    theme(axis.text.y = element_markdown(),
          axis.title.x = element_markdown()) +
    coord_flip()

sobol_inds_p <- (sobol_inds_p1 | sobol_inds_p2) +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "top")
sobol_inds_p

# save_plot("_plots/sobol-inds.pdf", sobol_inds_p, width = 8, height = 6)






