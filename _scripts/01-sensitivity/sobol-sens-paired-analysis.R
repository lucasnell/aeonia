#'
#' Small-scale sensitivity via Sobol indices
#'
#'

source("_scripts/01-sensitivity/00-preamble.R")
source("_scripts/01-sensitivity/sobol-preamble.R")
source("_scripts/01-sensitivity/sobol-sens-paired-analysis-functions.R")





if (!file.exists("_scripts/interm-data/sobol-sims-summs.rds")) {

    #' Directory containing output from `sobol-sens-paired.R`:
    sobol_dir <- "~/_globus"
    #' Output from `sobol-sens-paired.R`:
    sobol_sims <- paste0(sobol_dir, "/sobol-sims-paired.rds") |>
        read_rds()
    # Summarize each set of simulations:
    # Takes ~1 min  (multithreading doesn't help)
    sobol_summs <- imap(sobol_sims, \(sim_set, i) {
        out_df <- sim_set |>
            mutate(combo = factor(i, levels = 1:length(sobol_sims)),
                   sims = map(sims, \(s) s[1, names(vary_pars)])) |>
            unnest(sims) |>
            select(combo, everything())
        for (yv in c(yvars, "p_outbreak", "outbreak_size2", "sd_outbreak_size")) {
            out_df[[yv]] <- 0.0
        }
        for (j in 1:nrow(sim_set)) {
            for (yv in yvars) {
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

    }, .progress = .prog_args) |>
        list_rbind()

    write_rds(sobol_summs, "_scripts/interm-data/sobol-sims-summs.rds",
              compress = "gz")
    rm(sobol_sims); gc()

} else {

    sobol_summs <- read_rds("_scripts/interm-data/sobol-sims-summs.rds")

}




# Same thing but looking at differences between with and without Pseudo:
diff_sobol_summs <- sobol_summs |>
    group_by(combo, alate_dens, across(all_of(names(vary_pars)))) |>
    summarize(across(all_of(c(yvars, "p_outbreak", "outbreak_size2",
                              "sd_outbreak_size")),
                     \(x) x[n_pseudo > 0] - x[n_pseudo == 0]),
              .groups = "drop") |>
    mutate(across(starts_with("outbreak_size"), \(x) round(x, 2)))










# scatter(diff_sobol_summs |> filter(alate_dens == 1, abs(outbreak_size) > 0.05))



# # diff_sens_scat_p <-
# scatter(diff_sobol_summs, "outbreak_size")
# # scatter(diff_sobol_summs[[1L]], "outbreak_size")
#
# scatter(diff_sobol_summs[[2L]], "p_outbreak")
# # scatter(diff_sobol_summs[[1L]], "p_outbreak")


# sum(diff_sobol_summs[["outbreak_size"]] > 0)
# sum(diff_sobol_summs[["outbreak_size"]] == 0)
# sum(diff_sobol_summs[["outbreak_size"]] < 0)




# scatter(sobol_summs, .filter_vars = "alate_dens", .filter_conds = 1, "outbreak_size")
# scatter(sobol_summs, .filter_vars = "alate_dens", .filter_conds = 1, "sd_outbreak_size")

# sobol_summs |>
#     mutate(n_pseudo = factor(n_pseudo)) |>
#     ggplot(aes(outbreak_size, sd_outbreak_size, color = n_pseudo, fill = n_pseudo)) +
#     geom_point(alpha = 0.1) +
#     # stat_smooth(method = "lm", formula = y ~ poly(x,2,raw = TRUE)) +
#     stat_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
#                 se = TRUE, linewidth = 1) +
#     scale_color_manual(pretty_params("n_pseudo", TRUE),
#                        values = c("goldenrod", "dodgerblue"),
#                        aesthetics = c("color", "fill")) +
#     labs(x = yvar_desc[["outbreak_size"]] |> first_cap(),
#          y = yvar_desc[["sd_outbreak_size"]] |> first_cap()) +
#     facet_wrap(~ alate_dens) +
#     theme(axis.title = element_markdown(),
#           legend.title = element_markdown())




# QUESTION #1 ----
#' For parameter combos that seemed to result in a negative
#' effect of Pseudomonas, does alate ~ density affect outcomes?

# alate_dens_p <- diff_sobol_summs |>
#     group_by(combo) |>
#     summarize(without = outbreak_size[alate_dens == 0],
#               with = outbreak_size[alate_dens == 1]) |>
#     ggplot(aes(with, without)) +
#     geom_point(alpha = 0.1) +
#     geom_abline(slope = 1, intercept = 0, color = "red", linetype = "22") +
#     geom_hline(yintercept = 0, color = "gray70") +
#     geom_vline(xintercept = 0, color = "gray70") +
#     coord_equal(xlim = range(diff_sobol_summs$outbreak_size),
#                 ylim = range(diff_sobol_summs$outbreak_size)) +
#     labs(x = "Effect of *Pseudomonas* on outbreak size - alates ~ density",
#          y = "Effect of *Pseudomonas* on outbreak size - constant alates") +
#     theme(axis.title.x = element_markdown(),
#           axis.title.y = element_markdown())
# alate_dens_p
# save_plot("_plots/alate-dens.pdf", alate_dens_p, width = 6, height = 6)

# Answer: Yes, having the alate ~ density effect makes typically makes
# *Pseudomonas* less beneficial to plants






# LEFT OFF #2 ----
#' What parameter values are associated with Pseudomonas being bad for plants?




# EXTREMES ----
diff_sobol_summs |>
    filter(alate_dens == 1, spat_config %in% c(0, 1)) |>
    arrange(desc(outbreak_size)) |>
    select(combo, all_of(names(vary_pars)), outbreak_size)
diff_sobol_summs |>
    filter(alate_dens == 1, spat_config %in% c(0, 1)) |>
    arrange(outbreak_size) |>
    select(combo, all_of(names(vary_pars)), outbreak_size)



if (!file.exists("_scripts/interm-data/manip-extreme-combos.rds")) {

    # 7112 - Highest outbreak size difference
    # 9264 - Lowest outbreak size difference
    # Takes ~1 min
    set.seed(1783004255)
    combo_sims <- map(c(7112, 9264), one_combo_par_manip_simmer,
                      .progress = .prog_args)

    write_rds(combo_sims, "_scripts/interm-data/manip-extreme-combos.rds",
              compress = "gz")

} else {

    combo_sims <- read_rds("_scripts/interm-data/manip-extreme-combos.rds")

}




combo_p_list <- map(combo_sims, \(x) {
    one_combo_par_manip_plotter(x, .title = paste("combo", x[[1]][["combo"]][[1]]))
})

# combo_p_list[[1]]
# combo_p_list[[2]]

# combo_p_list[[1]] / combo_p_list[[2]] +
#     plot_layout(guides = "collect", axes = "collect")

# one_combo_par_manip_plotter(combo_sims[[1]]) |
#     one_combo_par_manip_plotter(combo_sims[[1]], .yvar = "log_wasps")



# combo_sims[[1]] |>
#     list_rbind() |>
#     filter(is.na(spat_config)) |>
#     filter(is.na(zeta)) |>
#     filter(outbreak_size == min(outbreak_size),
#            log_wasps > 3) |>
#     select(Y0, mean_N:zeta)
#
# combo_sims[[1]] |>
#     list_rbind() |>
#     filter(is.na(spat_config)) |>
#     filter(is.na(zeta)) |>
#     filter(is.na(Y0)) |>
#     ggplot(aes(log_wasps, outbreak_size)) +
#     geom_point(aes(color = n_pseudo)) +
#     scale_color_manual(pretty_params("n_pseudo", TRUE),
#                        values = np_pal)


# LEFT OFF - 17-Nov ----

# Comparing two parameter combinations:
sobol_summs |>
    filter(combo == 7112 | combo == 9264) |>
    filter(alate_dens == 1, n_pseudo > 0) |>
    select(n_pseudo, alate_dens, all_of(names(vary_pars)))



new_sims_args <- sobol_summs |>
    # filter(combo == 7112) |>
    filter(combo == 9264) |>
    filter(alate_dens == 1, n_pseudo > 0) |>
    select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
    mutate(n_sims = 1000, summ = "none") |>
    as.list()

set.seed(97504033)
old_sims <- do.call(one_combo, new_sims_args)
old_sims0 <- do.call(one_combo, list_assign(new_sims_args, n_pseudo = 0L))

new_arg <- list(Y0 = 4)
# set.seed(262255805)
new_sims <- do.call(one_combo, list_assign(new_sims_args, !!!new_arg))
new_sims0 <- do.call(one_combo, list_assign(new_sims_args, !!!c(new_arg, list(n_pseudo = 0L))))


.rep <- 1

sims_maxes <- c("aphids", "alates", "wasps") |>
    set_names() |>
    map(\(spp) {
        if (spp == "wasps") {
            z <- list(old_sims, old_sims0, new_sims, new_sims0) |>
                map(\(x) x |>
                        filter(rep == .rep) |>
                        org_timeseries_vars() |>
                        getElement("wasps")) |>
                do.call(what = c)
        } else {
            z <- c(old_sims[[spp]], old_sims0[[spp]], new_sims[[spp]],
                   new_sims0[[spp]])
            lgl <- !is.na(c(old_sims[["x"]], old_sims0[["x"]], new_sims[["x"]],
                            new_sims0[["x"]]))
            lgl2 <- c(old_sims[["rep"]], old_sims0[["rep"]], new_sims[["rep"]],
                      new_sims0[["rep"]]) == .rep
            if (spp == "aphids") {
                p <- "parasitized"
                z <- z + c(old_sims[[p]], old_sims0[[p]], new_sims[[p]],
                           new_sims0[[p]])
            }
            z <- z[lgl & lgl2]
        }
        return(max(z, na.rm = TRUE))
    })


old_p <- paired_timeseries(old_sims, old_sims0,
                           .rep = .rep,
                           .aphid_max = sims_maxes$aphids,
                           .alate_max = sims_maxes$alates,
                           .wasp_max = sims_maxes$wasps,
                           .tag = serify("",
                                         pretty_params(names(new_arg), TRUE),
                                         sprintf(" = %.2f (original)",
                                                 new_sims_args[[names(new_arg)]])))
new_p <- paired_timeseries(new_sims, new_sims0,
                           .rep = .rep,
                           .aphid_max = sims_maxes$aphids,
                           .alate_max = sims_maxes$alates,
                           .wasp_max = sims_maxes$wasps,
                           .tag = serify("",
                                         pretty_params(names(new_arg), TRUE),
                                         paste(" =", new_arg)))

old_p / new_p

bind_rows(old_sims, old_sims0, new_sims, new_sims0) |>
    filter(is.na(x)) |>
    select(n_pseudo, all_of(names(new_arg)), rep, time, wasps, virus) |>
    group_by(across(all_of(c(names(new_arg), "n_pseudo", "rep")))) |>
    summarize(time = time[wasps == max(wasps)],
              max_wasps = max(wasps),
              outbreak_size = max(virus),
              .groups = "drop_last") |>
    summarize(across(time:last_col(), mean),
              .groups = "drop")





# **Nov-19 - LEFT OFF ----
old_sims_od <- do.call(one_combo, list_assign(new_sims_args, out_dispersals = TRUE,
                                              summ = "all",
                                              # n_pseudo = 0,
                                              n_sims = 1000))
new_sims_od <- do.call(one_combo, list_assign(new_sims_args,
                                              !!!c(new_arg, list(out_dispersals = TRUE,
                                                                 summ = "all",
                                                                 # n_pseudo = 0,
                                                                 n_sims = 1000))))

old_sims_od |> getElement("outbreak_size") |> mean()
new_sims_od |> getElement("outbreak_size") |> mean()


# old_sims_od[["disps"]] |> reduce(`+`) |> (\(x) Matrix::Matrix(x / 1000, sparse = TRUE))()
# new_sims_od[["disps"]] |> reduce(`+`) |> (\(x) Matrix::Matrix(x / 1000, sparse = TRUE))()

# Total alates input to each plant:
# cbind(old = old_sims_od[["disps"]] |> reduce(`+`) |> (\(x) x / 1000)() |> rowSums(),
#       new = new_sims_od[["disps"]] |> reduce(`+`) |> (\(x) x / 1000)() |> rowSums())
old_sims_od[["disps"]] |> reduce(`+`) |> (\(x) x / 1000)() |> rowSums() |> matrix(3, 3)
new_sims_od[["disps"]] |> reduce(`+`) |> (\(x) x / 1000)() |> rowSums() |> matrix(3, 3)




# bind_rows(new_sims, new_sims0) |>
#     filter(is.na(x)) |>
#     group_by(n_pseudo, rep) |>
#     summarize(outbreak_size = max(virus), .groups = "drop")





# 2-par combo sims ----

if (!file.exists("_scripts/interm-data/manip-extreme-combos-2pars.rds")) {

    # Takes ~7.5 min w/ 7 combinations:
    set.seed(1547643870)
    combo_2par_sims <- list(c("Y0", "mean_N"),
                            c("Y0", "sd_N"),
                            c("mean_N", "sd_N"),
                            c("virus_attract", "pseudo_repel"),
                            c("pseudo_repel", "pseudo_surv"),
                            c("zeta", "pseudo_surv"),
                            c("zeta", "sd_N")) |>
        (\(x) set_names(x, map_chr(x, paste, collapse = "+")))() |>
        map(one_combo_2par_manip_simmer, .combo = 7112, progress_ = FALSE,
            .progress = .prog_args)
    write_rds(combo_2par_sims, "_scripts/interm-data/manip-extreme-combos-2pars.rds",
              compress = "gz")

} else {

    combo_2par_sims <- read_rds("_scripts/interm-data/manip-extreme-combos-2pars.rds")

}


cbind(names(combo_2par_sims))
one_combo_2par_manip_plotter(combo_2par_sims[["zeta+sd_N"]], .contour = TRUE)
one_combo_2par_manip_plotter(combo_2par_sims[["pseudo_repel+pseudo_surv"]], .contour = TRUE)


combo_2par_sims[["pseudo_repel+pseudo_surv"]] |>
    filter(n_pseudo == 3) |>
    ggplot(aes(pseudo_surv, outbreak_size, color = pseudo_repel)) +
    # geom_point() +
    geom_line(aes(group = factor(pseudo_repel))) +
    labs(x = pretty_params("pseudo_surv"), y = "Outbreak size") +
    scale_color_viridis_c(serify("", pretty_params("pseudo_repel", TRUE), ""),
                          option = "inferno") +
    theme(axis.title.x = element_markdown(),
          legend.title = element_markdown())
combo_2par_sims[["pseudo_repel+pseudo_surv"]] |>
    filter(n_pseudo == 3) |>
    ggplot(aes(pseudo_repel, outbreak_size, color = pseudo_surv)) +
    # geom_point() +
    geom_line(aes(group = factor(pseudo_surv))) +
    labs(x = pretty_params("pseudo_repel"), y = "Outbreak size") +
    scale_color_viridis_c(serify("", pretty_params("pseudo_surv", TRUE), ""),
                          option = "inferno") +
    theme(axis.title.x = element_markdown(),
          legend.title = element_markdown())


# save_plot("_plots/combo-plot1.pdf", combo_p_list[[1]], width = 8, height = 5)
# save_plot("_plots/combo-plot2.pdf", combo_p_list[[2]], width = 8, height = 5)



max_diff_args <- list(alate_dens = 1L,
                      Y0 = 2.5,
                      mean_N = 25,
                      sd_N = 0,
                      K = 23e3,
                      virus_attract = 1,
                      pseudo_repel = 1,
                      pseudo_surv = 0.85,
                      zeta = 0.0,
                      spat_config = 1L,
                      n_sims = 1e3L,
                      n_pseudo = 3L,
                      summ = "none")
max_diff_sims <- do.call(one_combo, max_diff_args)
max_diff_sims0 <- do.call(one_combo, list_assign(max_diff_args, n_pseudo = 0L))

# paired_timeseries(max_diff_sims, max_diff_sims0, .title = "Original")


max_diff_sims2 <- do.call(one_combo, list_assign(max_diff_args, pseudo_surv = 0.95))
max_diff_sims02 <- do.call(one_combo, list_assign(max_diff_args, pseudo_surv = 0.95, n_pseudo = 0L))

# Make adult alates as susceptible to parasitism as adult apterous:
test_R <- pop_info$R[c(1,2,3,2)]
max_diff_sims3 <- do.call(one_combo, list_assign(max_diff_args, R = test_R))
max_diff_sims03 <- do.call(one_combo, list_assign(max_diff_args, R = test_R,
                                                  n_pseudo = 0L))


bind_rows(calc_metrics(max_diff_sims), calc_metrics(max_diff_sims0),
          calc_metrics(max_diff_sims2), calc_metrics(max_diff_sims02),
          calc_metrics(max_diff_sims3), calc_metrics(max_diff_sims03)) |>
    mutate(n_pseudo = c(3, 0, 3, 0, 3, 0),
           type = c(rep("orig", 2), rep("mod-psi", 2), rep("mod-R", 2))) |>
    select(type, n_pseudo, everything())


#'
#' R has no effect on how Pseudomonas affects outbreaks!
#' Making adult alates as susceptible to parasitism as adult apterous
#' simply makes outbreaks smaller for both with and without Pseudomonas equally.
#' (Confirmed with 10e3 sims.)
#'



# paired_timeseries(max_diff_sims, max_diff_sims0, .tag = "Original",
#                   .aphid_max = 1389, .alate_max = 94, .wasp_max = 39) /
#     paired_timeseries(max_diff_sims2, max_diff_sims02,
#                       .tag = serify("", "&psi;", " = 0.95"),
#                       .aphid_max = 1389, .alate_max = 94, .wasp_max = 39) &
#     theme(panel.grid.major.y = element_line(color = "gray80"))
#
#
# paired_p_alate_hist(max_diff_sims, max_diff_sims0)


# LEFT OFF - 7-Nov ----

# Redo sims with out_stages = TRUE:
stg_max_diff_sims <- max_diff_sims[1,] |>
    select(any_of(names(formals(one_combo)))) |>
    mutate(n_sims = 1e2, out_stages = TRUE, summ = "none") |>
    as.list() |>
    do.call(what = one_combo)
stg_max_diff_sims0 <- max_diff_sims0[1,] |>
    select(any_of(names(formals(one_combo)))) |>
    mutate(n_sims = 1e2, out_stages = TRUE, summ = "none") |>
    as.list() |>
    do.call(what = one_combo)
stg_max_diff_sims2 <- max_diff_sims2[1,] |>
    select(any_of(names(formals(one_combo)))) |>
    mutate(n_sims = 1e2, out_stages = TRUE, summ = "none") |>
    as.list() |>
    do.call(what = one_combo)
stg_max_diff_sims02 <- max_diff_sims02[1,] |>
    select(any_of(names(formals(one_combo)))) |>
    mutate(n_sims = 1e2, out_stages = TRUE, summ = "none") |>
    as.list() |>
    do.call(what = one_combo)


paired_attack_plots(stg_max_diff_sims, stg_max_diff_sims0,
                    .title = serify("", "&psi;", " = 0.85"))
paired_attack_plots(stg_max_diff_sims2, stg_max_diff_sims02,
                    .title = serify("", "&psi;", " = 0.95"))

paired_timeseries(stg_max_diff_sims, stg_max_diff_sims0,
                  .tag = serify("", "&psi;", " = 0.85"),
                  .aphid_max = 1389, .alate_max = 94, .wasp_max = 39) /
    paired_timeseries(stg_max_diff_sims2, stg_max_diff_sims02,
                      .tag = serify("", "&psi;", " = 0.95"),
                      .aphid_max = 1389, .alate_max = 94, .wasp_max = 39) &
    theme(panel.grid.major.y = element_line(color = "gray80"))

paired_attack_plots(stg_max_diff_sims, stg_max_diff_sims0,
                  .tag = serify("", "&psi;", " = 0.85"),
                  .y_min = 0.43) /
    paired_attack_plots(stg_max_diff_sims2, stg_max_diff_sims02,
                      .tag = serify("", "&psi;", " = 0.95"),
                      .y_min = 0.43) &
    theme(panel.grid.major.y = element_line(color = "gray80"))



stg_max_diff_sims |>
    filter(rep == 1) |>
    select(zeta, n_pseudo, rep:last_col()) |>
    org_timeseries_vars()




 max_diff_sims_ob <- max_diff_sims |>
    select(rep:last_col()) |>
    filter(!is.na(wasps)) |>
    group_by(rep) |>
    summarize(outbreak_size = max(virus)) |>
    getElement("outbreak_size") |>
    mean()
max_diff_sims0_ob <- max_diff_sims0 |>
    select(rep:last_col()) |>
    filter(!is.na(wasps)) |>
    group_by(rep) |>
    summarize(outbreak_size = max(virus)) |>
    getElement("outbreak_size") |>
    mean()

mean(max_diff_sims_ob); mean(max_diff_sims0_ob)
print_diff_mean_w_boot_ci(max_diff_sims_ob, max_diff_sims0_ob)





# largest so far: 5.29   (5.15 - 5.42)
#' alate_dens = 1L,
#' Y0 = 2.5,
#' mean_N = 25,
#' sd_N = 0,
#' K = 23e3,
#' virus_attract = 1,
#' pseudo_repel = 1,
#' pseudo_surv = 0.85,
#' zeta = 0.0,
#' spat_config = 1L


paired_target_sims(7112L, .n_sims = 1000, .summ = "all",
                   pseudo_repel = 2,
                   zeta = 0.2) |>
    split(~ zeta +
              pseudo_repel +
              n_pseudo, drop = TRUE) |>
    map(\(x) calc_metrics(x) |>
            mutate(zeta = x$zeta[[1]],
                   pseudo_repel = x$pseudo_repel[[1]],
                   n_pseudo = x$n_pseudo[[1]])) |>
    list_rbind() |>
    arrange(zeta,
            pseudo_repel,
            n_pseudo) |>
    select(zeta,
           pseudo_repel,
           n_pseudo, everything())


# With just changing zeta = 0.2

# # A tibble: 4 × 6
#     zeta n_pseudo outbreak_size log_aphids log_alates log_wasps
#    <dbl>    <int>         <dbl>      <dbl>      <dbl>     <dbl>
# 1 0.0322        0          4.24       5.04       1.81      3.01
# 2 0.0322        3          7.27       5.22       2.06      3.11
# 3 0.2           0          4.25       5.03       1.80      3.01
# 4 0.2           3          5.92       5.09       1.90      3.00


sims <- paired_target_sims(7112L, zeta = 0.5)




# set.seed(1032439295)
sim_df <- sobol_summs |>
    filter(combo == i) |>
    filter(n_pseudo > 0, alate_dens == 1) |>
    select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
    as.list() |>
    c(list(n_sims = 4, summ = "none")) |>
    do.call(what = one_combo) |>
    select(n_pseudo, rep:wasps)
sim0_df <- sobol_summs |>
    filter(combo == i) |>
    filter(n_pseudo == 0, alate_dens == 1) |>
    select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
    as.list() |>
    c(list(n_sims = 4, summ = "none")) |>
    do.call(what = one_combo) |>
    select(n_pseudo, rep:wasps)

sim_df |> calc_metrics()
sim0_df |> calc_metrics()

paired_timeseries(sim_df, sim0_df, .title = sprintf("combo %i", i))
paired_timeseries(sim_df, sim0_df, TRUE, 4, .title = sprintf("combo %i", i))



diff_sobol_summs |>
    filter(alate_dens == 1, spat_config == 1) |>
    ggplot(aes(zeta, Y0)) +
    geom_point(aes(color = outbreak_size)) +
    scale_color_viridis_c()

diff_sobol_summs |>
    filter(alate_dens == 1, spat_config == 1) |>
    ggplot(aes(zeta, outbreak_size)) +
    geom_point(alpha = 0.1) +
    # stat_smooth(method = "gam",
    #             formula = y ~ s(x, bs = "cs"),
    #             se = TRUE, linewidth = 1) +
    labs(x = pretty_params("zeta"), y = yvar_desc[["outbreak_size"]]) +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown())




# LEFT OFF #3 ----
# Why does the effect of Pseudomonas on P(outbreak) not coincide closely with
# its effect on outbreak size?


# names(vary_pars)
# [1] "Y0"            "mean_N"        "sd_N"          "K"             "virus_attract"
# [6] "pseudo_repel"  "pseudo_surv"   "zeta"          "spat_config"

dd <- diff_sobol_summs[[2L]] |>
    filter(spat_config == 0) |>
    mutate(spat_config = factor(spat_config),
           resid = residuals(lm(outbreak_size ~ p_outbreak)))

dd |>
    ggplot(aes(log(outbreak_size), p_outbreak)) +
    geom_point() +
    # facet_wrap(~ spat_config) +
    scale_color_viridis_c()

mod <- lm(resid ~ K * Y0, dd)
mod |> summary()


dd |>
    # group_by(spat_config) |>
    # mutate(resid = residuals(lm(outbreak_size ~ p_outbreak))) |>
    # ungroup() |>
    mutate(resid = residuals(mod)) |>
    ggplot(aes(resid, Y0)) +
    geom_point()




# sens_diff_scat_p <- scatter(diff_sobol_summs, .title = "")
# save_plot("_plots/sens-diff-scatter.pdf", sens_diff_scat_p,
#           width = 8, height = 5)
# save_plot("_plots/sens-diff-scatter.png", sens_diff_scat_p, dpi = 150,
#           width = 8, height = 5)


# sens_scat_p <- scatter(sobol_summs[[2L]], "outbreak_size", .title = "")
# save_plot("_plots/sens-scatter.pdf", sens_scat_p, width = 8, height = 5)
#
# sens_scat0_p <- scatter(sobol_summs[[1L]], "outbreak_size", .title = "")
# save_plot("_plots/sens-scatter-b0.pdf", sens_scat0_p, width = 8, height = 5)
#
# # scatter2(sobol_summs[[1L]])
# sens_scat2_p <- scatter2(sobol_summs[[2L]], .title = "")
#
# save_plot("_plots/sens-scatter2.pdf", sens_scat2_p, width = 8, height = 5)



spat_config_p <- sobol_summs[[2]] |>
    mutate(spat_config = factor(spat_config), n_pseudo = factor(n_pseudo)) |>
    ggplot(aes(spat_config, outbreak_size, color = n_pseudo)) +
    geom_violin() +
    geom_hline(data = sobol_summs[[2]] |>
                   mutate(spat_config = factor(spat_config), n_pseudo = factor(n_pseudo)) |>
                   filter(spat_config == 1) |>
                   group_by(spat_config, n_pseudo) |>
                   summarize(outbreak_size = mean(outbreak_size), .groups = "drop"),
               aes(yintercept = outbreak_size, color = n_pseudo),
               linetype = "22") +
    stat_summary(fun = "mean", size = 4, geom = "point") +
    labs(x = pretty_params("spat_config") |> str_to_sentence(),
         y = yvar_desc[["outbreak_size"]] |> str_to_sentence()) +
    scale_color_manual(pretty_params("n_pseudo", TRUE),
                       values = c("goldenrod", "dodgerblue")) +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown(),
          legend.title = element_markdown(),
          plot.title = element_markdown())

# save_plot("_plots/spat-config.pdf", spat_config_p, width = 6, height = 4)

diff_spat_config_p <- diff_sobol_summs[[2]] |>
    mutate(spat_config = factor(spat_config)) |>
    ggplot(aes(spat_config, outbreak_size)) +
    geom_violin() +
    geom_hline(data = diff_sobol_summs[[2]] |>
                   mutate(spat_config = factor(spat_config)) |>
                   filter(spat_config == 1) |>
                   summarize(outbreak_size = mean(outbreak_size), .groups = "drop"),
               aes(yintercept = outbreak_size),
               linetype = "22") +
    stat_summary(fun = "mean", size = 4, geom = "point") +
    labs(x = pretty_params("spat_config") |> str_to_sentence(),
         y = paste("Effect of *Pseudomonas* on", yvar_desc[["outbreak_size"]])) +
    theme(axis.title.y = element_markdown(),
          axis.title.x = element_markdown(),
          legend.title = element_markdown(),
          plot.title = element_markdown())

# save_plot("_plots/spat-config-diff.pdf", diff_spat_config_p, width = 6, height = 4)




# multi_scatter(sobol_summs[[2L]], np = 0L, .N = 2^8, axis.text = element_markdown(size = 6),
#               axis.title = element_markdown(family = "serif", size = 8))
# multi_scatter(sobol_summs[[2L]], np = 3L, .N = 2^8, axis.text = element_markdown(size = 6),
#               axis.title = element_markdown(family = "serif", size = 8))




set.seed(1998643658)
sobol_inds_np3 <- map(sobol_summs, \(ss) {
    Y <- ss[["outbreak_size"]][ss$n_pseudo > 0]
    ind <- sobol_indices(Y = Y, N = N, params = pretty_params(names(vary_pars)),
                         boot = TRUE, R = 2000)
    return(ind)
})

set.seed(512925036)
sobol_inds_diff <- map(diff_sobol_sims, \(ss) {
    Y <- ss[["outbreak_size"]]
    ind <- sobol_indices(Y = Y, N = N, params = pretty_params(names(vary_pars)),
                         boot = TRUE, R = 2000)
    return(ind)
})


sobol_inds_p1 <- plot(sobol_inds_np3[[2]]) +
    labs(y = "Sobol' index (outbreak size with<br>three *Pseudomonas*)") +
    theme(axis.text.y = element_markdown(),
          axis.title.x = element_markdown()) +
    coord_flip()
sobol_inds_p2 <- plot(sobol_inds_diff[[2]]) +
    labs(y = "Sobol' index (effect of *Pseudomonas*<br>on outbreak size)") +
    theme(axis.text.y = element_markdown(),
          axis.title.x = element_markdown()) +
    coord_flip()

sobol_inds_p <- (sobol_inds_p1 | sobol_inds_p2) +
    plot_layout(guides = "collect", axes = "collect") &
    theme(legend.position = "top")

save_plot("_plots/sobol-inds.pdf", sobol_inds_p, width = 8, height = 6)






