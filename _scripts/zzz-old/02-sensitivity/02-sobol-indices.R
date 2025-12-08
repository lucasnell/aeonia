#'
#' Small-scale sensitivity via Sobol indices
#'
#'

source("_scripts/00-preamble.R")
source("_scripts/02-sensitivity/00-sobol-preamble.R")


out_files <- list(summs = "sobol-sims-summs.rds",
                  diff_summs = "sobol-sims-diff-summs.rds") |>
    map(\(x) paste0("_scripts/interm-data/", x))




# =============================================================================*
# =============================================================================*
# Sim summaries ----
# =============================================================================*
# =============================================================================*

if (!file.exists(out_files$summs)) {


    proc_one_file <- function(j, prog) {
        # j = 1
        # rm(j, sims_j, summs_j)

        # Output from `sobol-sens-paired.R`:
        sims_j <- sprintf("~/_globus/sobol-sims-paired-%i.rds", j) |>
            read_rds()

        # Start for rows in this job's simulations:
        curr_start <- (j - 1L) * (obj_env$n_rows %/% obj_env$n_files) + 1L

        one_set <- function(sim_set, i) {
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

        }

        summs_j <- imap(sims_j, one_set) |>
            list_rbind()

        prog()

        return(summs_j)

    }

    proc_all_files <- function() {

        suppressPackageStartupMessages(library(future.apply))
        suppressPackageStartupMessages(library(progressr))
        handlers("progress")
        with(plan(multisession, workers = min(4L, options()[["mc.cores"]]),
                  gc = TRUE), local = TRUE)
        oopts <- options(future.globals.maxSize = 1.5e9)  ## 1.5 GB
        on.exit(options(oopts))

        with_progress({
            p <- progressor(nrow(sobol_mat))
            out <- future_lapply(1:(obj_env$n_files),
                                 proc_one_file,
                                 prog = p,
                                 future.globals = c("obj_env", "proc_one_file"),
                                 future.packages = "tidyverse") |>
                list_rbind()
        })
        return(out)
    }

    obj_env <- list()
    obj_env$n_rows <- nrow(sobol_mat)
    obj_env$n_files <- length(list.files("~/_globus", "sobol-sims-paired-.?.rds"))
    obj_env$vary_pars <- vary_pars
    obj_env$yvars <- yvars

    # Summarize each set of simulations:
    # Takes ~13 min w/ 4 threads
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


# =============================================================================*
# =============================================================================*
# Fit Sobol' indices ----
# =============================================================================*
# =============================================================================*

set.seed(1998643658)
sobol_inds_np3 <- sobol_summs |>
    filter(alate_dens == 1) |>
    (\(ss) {
        Y <- ss[["outbreak_size"]][ss$n_pseudo > 0]
        ind <- sobol_indices(Y = Y, N = N, params = names(vary_pars),
                             boot = TRUE, R = 2000)
        return(ind)
    })()

set.seed(512925036)
sobol_inds_diff <- diff_sobol_summs |>
    filter(alate_dens == 1) |>
    (\(ss) {
        Y <- ss[["outbreak_size"]]
        ind <- sobol_indices(Y = Y, N = N, params = names(vary_pars),
                             boot = TRUE, R = 2000)
        return(ind)
    })()


d = tibble(a = 1:3, b = 4:6)

foo <- function(d, x) {
    select(d, {{ x }})
}

foo(d, "b")

# rm(foo, d)


# =============================================================================*
# =============================================================================*
# plot Sobol' indices ----
# =============================================================================*
# =============================================================================*
sobol_plotter <- function(x, y_lab = "Sobol' index", .reorder_vals = NULL) {

    # x = sobol_inds_np3; y_lab = "Sobol' index"; .reorder_vals = "original"
    # rm(x, y_lab, .reorder_vals, dt, gg)

    stopifnot(inherits(x, "sensobol"))

    dt <- x$results |>
        as_tibble() |>
        filter(sensitivity %in% c("Si", "Ti")) |>
        mutate(parameters = pretty_params(parameters) |>
                   factor(levels = pretty_params(names(vary_pars))))

    stopifnot(nrow(dt) > 1)

    if (!is.null(.reorder_vals)) {
        stopifnot(length(.reorder_vals) %in% c(1L, nrow(dt)))
        if ((is.character(.reorder_vals) || is.numeric(.reorder_vals)) &&
            length(.reorder_vals) == nrow(dt)) {
            dt$parameters <- fct_reorder(dt$parameters, .reorder_vals)
        } else if (is.character(.reorder_vals) && length(.reorder_vals) == 1L) {
            if (! .reorder_vals %in% colnames(dt)) {
                stop("\nif .reorder_vals is a single string, it must refer ",
                     "to a column in x$results")
            }
            dt$parameters <- fct_reorder(dt$parameters, dt[[.reorder_vals]])
        } else {
            stop("\n.reorder_vals must be a single string or a numeric / ",
                 "character vector of same length as nrow(x$results)")
        }
    }

    gg <- ggplot(dt, aes(parameters, original, fill = sensitivity)) +
        geom_bar(stat = "identity",
                 position = position_dodge(0.6),
                 color = NA) +
        scale_y_continuous(breaks = scales::pretty_breaks(n = 3)) +
        labs(x = "", y = y_lab) +
        scale_fill_viridis_d(name = "Estimator order:",
                             labels = c("first", "total"),
                             option = "cividis", begin = 0.35, end = 0.9, direction = -1)
    if (any(grepl("high.ci", colnames(dt))) == TRUE) {
        gg <- gg +
            geom_errorbar(aes(ymin = low.ci, ymax = high.ci),
                          position = position_dodge(0.6),
                          width = 0.3, linewidth = 0.75)
    }
    gg +
        theme(legend.text = element_markdown(),
              axis.text.y = element_markdown(),
              axis.title.x = element_markdown()) +
        coord_flip()

}



sobol_inds_p1 <- sobol_plotter(sobol_inds_np3, .reorder_vals = "original",
                               "Sobol' index (outbreak size with<br>three *Pseudomonas* patches)")
sobol_inds_p2 <- sobol_plotter(sobol_inds_diff, .reorder_vals = sobol_inds_np3$results$original,
                               "Sobol' index (effect of *Pseudomonas*<br>on outbreak size)")

sobol_inds_p <- (sobol_inds_p1 | sobol_inds_p2) +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "top")
sobol_inds_p

# save_plot("_plots/sobol-inds.pdf", sobol_inds_p, width = 8, height = 6)






