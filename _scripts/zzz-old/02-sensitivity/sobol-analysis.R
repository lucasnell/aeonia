#'
#' Small-scale sensitivity via Sobol indices
#'
#'

source("_scripts/00-preamble.R")
source("_scripts/02-sensitivity/00-sobol-preamble.R")
source("_scripts/02-sensitivity/sobol-analysis-functions.R")


sobol_summs <- read_rds("_scripts/interm-data/sobol-sims-summs.rds")

diff_sobol_summs <- read_rds("_scripts/interm-data/sobol-sims-diff-summs.rds")








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




# =============================================================================*
# alate ~ density effect? ----
# =============================================================================*

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






# =============================================================================*
# Extreme param combos ----
# =============================================================================*

#' What parameter values are associated with Pseudomonas being especially bad
#' or good for plants?


diff_sobol_summs |>
    filter(alate_dens == 1, spat_config %in% c(0, 1)) |>
    arrange(desc(outbreak_size)) |>
    select(combo, all_of(names(vary_pars)), outbreak_size)
diff_sobol_summs |>
    filter(alate_dens == 1, spat_config %in% c(0, 1)) |>
    arrange(outbreak_size) |>
    select(combo, all_of(names(vary_pars)), outbreak_size)



extreme_combos <- diff_sobol_summs |>
    filter(alate_dens == 1, spat_config %in% c(0, 1)) |>
    filter(outbreak_size %in% range(outbreak_size)) |>
    arrange(desc(outbreak_size)) |>
    getElement("combo") |>
    set_names(c("high", "low"))


# Comparing two most extreme parameter combinations:
diff_sobol_summs |>
    filter(combo %in% extreme_combos) |>
    filter(alate_dens == 1) |>
    select(combo, outbreak_size, all_of(names(vary_pars)))

sobol_summs |>
    filter(combo %in% extreme_combos) |>
    filter(alate_dens == 1) |>
    select(combo, n_pseudo, outbreak_size, all_of(names(vary_pars))) |>
    rename(OBS = outbreak_size, virus_attr = virus_attract)



set.seed(1381664252)
extreme_sims <- sobol_summs |>
    filter(combo %in% extreme_combos) |>
    filter(alate_dens == 1) |>
    select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
    mutate(spat_config = 1,  # << for easier interpretation (doesn't affect outcome)
           n_sims = 10, summ = "none", out_stages = TRUE) |>
    (\(x) {
        x$sims  <- pmap(x, one_combo) |>
            map(\(x) select(x, rep:last_col()))
        new_pars <- c("n_pseudo", "zeta")
        for (i in 1:nrow(x)) {
            for (z in new_pars) {
                x$sims[[i]][[z]] <- x[[z]][[i]]
            }
            x$sims[[i]] <- select(x$sims[[i]], all_of(new_pars), everything())
        }
        return(x)
    })()



paired_timeseries(extreme_sims$sims[[1]], extreme_sims$sims[[2]], .tag = "High",
                  .aphid_max = 1610, .alate_max = 301, .wasp_max = 50) /
    paired_timeseries(extreme_sims$sims[[3]], extreme_sims$sims[[4]],
                      .tag = "Low",
                      .aphid_max = 1610, .alate_max = 301, .wasp_max = 50)


paired_attack_plots(extreme_sims$sims[[1]], extreme_sims$sims[[2]], .tag = "High",
                    .y_min = 0.4) /
    paired_attack_plots(extreme_sims$sims[[3]], extreme_sims$sims[[4]], .tag = "Low",
                        .y_min = 0.4)





infect_time <- extreme_sims$sims[[3]] |>
    filter(rep == 1, x == 2, y == 2) |>
    filter(virus == 1) |>
    getElement("time") |>
    getElement(1)


# Weird changes in attack survival are due to influx of adult alates when
# other aphid densities are low:
extreme_sims$sims[[3]] |>
    filter(rep == 1) |>
    filter(time > infect_time - 5, time < infect_time + 12) |>
    mutate(plant = interaction(x, y, drop = TRUE)) |>
    split(~ time, drop = TRUE) |>
    map(add_Y_A,
        zeta = extreme_sims$sims[[3]][["zeta"]][[1]],
        a = pop_info[["a"]],
        h = pop_info[["h"]],
        k = pop_info[["k"]],
        R = pop_info[["R"]],
        keep_stages = TRUE) |>
    list_rbind() |>
    filter(plant == "2.2") |>
    select(time, aphids_juv:wasps, A) |>
    pivot_longer(aphids_juv:A) |>
    ggplot(aes(time, value)) +
    geom_line(aes(color = name)) +
    geom_point(aes(color = name)) +
    geom_vline(xintercept = infect_time, linetype = "22") +
    geom_vline(xintercept = c(36, 38), linetype = "22", color = "gray70") +
    scale_color_viridis_d(option = "turbo", guide = "none") +
    facet_wrap(~ name, scales = "free")



# manip extremes ----

if (!file.exists("_scripts/interm-data/manip-extreme-combos.rds")) {
    # Takes ~1 min
    set.seed(1783004255)
    combo_sims <- map(extreme_combos, one_combo_par_manip_simmer,
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

# save_plot("_plots/combo-plot1.pdf", combo_p_list[[1]], width = 8, height = 5)
# save_plot("_plots/combo-plot2.pdf", combo_p_list[[2]], width = 8, height = 5)


# 2-par manip extremes ----

if (!file.exists("_scripts/interm-data/manip-extreme-combos-2pars.rds")) {

    # Takes ~15 min w/ 2 sets of 7 combinations:
    set.seed(1547643870)
    combo_2par_sims <- extreme_combos |>
        map(\(Z) {
            list(c("Y0", "mean_N"),
                 c("Y0", "sd_N"),
                 c("mean_N", "sd_N"),
                 c("virus_attract", "pseudo_repel"),
                 c("pseudo_repel", "pseudo_surv"),
                 c("zeta", "pseudo_surv"),
                 c("zeta", "sd_N")) |>
                (\(x) set_names(x, map_chr(x, paste, collapse = "+")))() |>
                map(one_combo_2par_manip_simmer, .combo = Z,
                    progress_ = FALSE, .spat_config = 1L,
                    .progress = .prog_args)
        })

    write_rds(combo_2par_sims, "_scripts/interm-data/manip-extreme-combos-2pars.rds",
              compress = "gz")

} else {

    combo_2par_sims <- read_rds("_scripts/interm-data/manip-extreme-combos-2pars.rds")

}

# make sure these are both within range [-5, 5]:
combo_2par_sims |>
    map(\(x) {
        map(x,
            \(xx) {
                xx |>
                    group_by(across(2:3))  |>
                    summarize(outbreak_size = outbreak_size[n_pseudo != "0"] -
                                  outbreak_size[n_pseudo == "0"],
                              .groups = "drop") |>
                    mutate(outbreak_size = round(outbreak_size, 3)) |>
                    getElement("outbreak_size") |>
                    range()
            }) |>
            do.call(what = rbind)
    }) |>
    do.call(what = rbind) |>
    (\(xx) return(c(min(xx[,1]), max(xx[,2]))))()


# Combinations that result in greater outbreak size with Pseudomonas:
combo_2par_sims |>
    map(\(.df) {
        .df |>
            list_rbind() |>
            mutate(grp = map(1:(n() %/% 2L), \(i) return(rep(i, 2))) |> list_c()) |>
            group_by(grp, .drop = FALSE) |>
            summarize(outbreak_size = outbreak_size[n_pseudo != 0] -
                          outbreak_size[n_pseudo == 0],
                      across(any_of(names(vary_pars)), mean)) |>
            mutate(outbreak_size = round(outbreak_size, 3)) |>
            filter(outbreak_size > 0) |>
            arrange(desc(outbreak_size))
    })


cbind(names(combo_2par_sims[[1]]))
#      [,1]
# [1,] "Y0+mean_N"
# [2,] "Y0+sd_N"
# [3,] "mean_N+sd_N"
# [4,] "virus_attract+pseudo_repel"
# [5,] "pseudo_repel+pseudo_surv"
# [6,] "zeta+pseudo_surv"
# [7,] "zeta+sd_N"


which_2par <- "Y0+mean_N"

(one_combo_2par_manip_plotter(combo_2par_sims[[1]][[which_2par]],
                              .contour = TRUE, .tag = "Highest")) /
    (one_combo_2par_manip_plotter(combo_2par_sims[[2]][[which_2par]],
                                  .contour = TRUE, .tag = "Lowest")) +
    plot_layout(guides = "collect")






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




# For below, I'm manually manipulating parameters to go from lowest effect of
# Pseudomonas on outbreak size, to highest (or close)
# **Manipulate pars ----

# # A tibble: 2 × 10
#   combo   *obs    Y0 mean_N  sd_N      K   † va    ‡ pr    ¶ ps   zeta
#   <fct>  <dbl> <dbl>  <dbl> <dbl>  <dbl>  <dbl>   <dbl>   <dbl>  <dbl>
# 1 7064    2.93  2.15   16.5  4.43 18358.   1.38    5.99   0.852 0.0159
# 2 34255  -3.98  1.65   87.7 24.6  15634.   2.97    9.03   0.892 0.816
#
# * = outbreak_size
# † = virus_attract
# ‡ = pseudo_repel
# ¶ = pseudo_surv


old_sims_args <- sobol_summs |>
    # filter(combo == extreme_combos[["high"]]) |>
    filter(combo == extreme_combos[["low"]]) |>
    filter(alate_dens == 1, n_pseudo > 0) |>
    select(n_pseudo, alate_dens, all_of(names(vary_pars))) |>
    mutate(n_sims = 1000, summ = "none", spat_config = -1L) |>
    as.list()

set.seed(97504033)
old_sims <- do.call(one_combo, old_sims_args)
old_sims0 <- do.call(one_combo, list_assign(old_sims_args, n_pseudo = 0L))

# Values for candidate arguments that could flip outcome:
new_args <- list(zeta = 0,
                 mean_N = 15,
                 # Y0 = 2,
                 virus_attract = 1,
                 pseudo_repel = 1,
                 sd_N = 0)


# Just takes a few sec
set.seed(1826998848)
multi_pert_sims <- 0:length(new_args) |>
    map(\(k) {
        # k = 2L
        # rm(k, m, out, i, arg_list_i, sims)
        if (k == 0) {
            out <- tibble(k = 0L,
                          params = "",
                          n_pseudo = c(0L, old_sims_args[["n_pseudo"]]),
                          outbreak_size = 0)
            out$outbreak_size <- list(old_sims0, old_sims) |>
                map_dbl(\(x) {
                    x |>

                        filter(is.na(x)) |>
                        group_by(rep) |>
                        summarize(outbreak_size = max(virus)) |>
                        getElement("outbreak_size") |>
                        mean()
                })
            return(out)
        }
        m <- combn(length(new_args), k)
        out <- crossing(k = k,
                        params = map(1:ncol(m), \(i) names(new_args)[m[,i]]),
                        n_pseudo = c(0L, old_sims_args[["n_pseudo"]]),
                        outbreak_size = 0)
        for (i in 1:nrow(out)) {
            arg_list_i <- list_assign(old_sims_args, !!!new_args[out$params[[i]]])
            arg_list_i[["n_pseudo"]] <- out$n_pseudo[[i]]
            arg_list_i[["summ"]] <- "all"
            sims <- do.call(one_combo, arg_list_i)
            out$outbreak_size[[i]] <- mean(sims$outbreak_size)
        }
        out$params <- map_chr(out$params, \(x) paste(x, collapse = " & "))
        return(out)
    }) |>
    list_rbind() |>
    mutate(n_pseudo = factor(n_pseudo))


old_obs <- list()
old_obs[["0"]] <- multi_pert_sims |>
    filter(n_pseudo == 0, k == 0) |>
    getElement("outbreak_size")
old_obs[["3"]] <- multi_pert_sims |>
    filter(n_pseudo == 3, k == 0) |>
    getElement("outbreak_size")
old_obs[["diff"]] <- old_obs[["3"]] - old_obs[["0"]]



multi_pert_plotter <- function(d, ob_lim = NULL, dob_lim = NULL) {
    dd <- d |>
        group_by(k, params) |>
        summarize(outbreak_size = outbreak_size[n_pseudo == 3] -
                      outbreak_size[n_pseudo == 0],
                  .groups = "drop")
    p1 <- d |>
        ggplot(aes(outbreak_size)) +
        geom_vline(xintercept = old_obs[["0"]], color = np_pal[["0"]],
                   linetype = "22") +
        geom_vline(xintercept = old_obs[["3"]], color = np_pal[["3"]],
                   linetype = "22") +
        geom_segment(aes(xend = origin, y = params, color = n_pseudo),
                     position = position_dodge(width = 0.3)) +
        geom_point(aes(y = params, color = n_pseudo),
                   position = position_dodge(width = 0.3)) +
        coord_cartesian(xlim = ob_lim) +
        labs(x = "Outbreak size") +
        scale_color_manual(pretty_params("n_pseudo", TRUE), values = np_pal)
    p2 <- dd |>
        ggplot(aes(outbreak_size)) +
        geom_vline(xintercept = 0, color = "gray70", linewidth = 1) +
        geom_vline(xintercept = old_obs[["diff"]], color = "gray70", linetype = "22") +
        geom_segment(aes(xend = old_obs[["diff"]], y = params)) +
        geom_point(aes(y = params)) +
        coord_cartesian(xlim = dob_lim) +
        labs(x = "Effect of *Pseudomonas*<br>on outbreak size")

    (p1 | p2) +
        plot_layout(guides = "collect", axes = "collect") &
        theme(axis.text.y = element_markdown(),
              axis.title.y = element_blank(),
              axis.title.x = element_markdown(),
              legend.title = element_markdown())
}



multi_pert_p <- multi_pert_sims |>
    filter(k == length(new_args) | k == (length(new_args)-1)) |>
    mutate(params = map2_chr(k, params,
                             \(k, params) {
                                 if (k == length(new_args)) return("all")
                                 p <- str_split(params, " & ")[[1]]
                                 no <- names(new_args)[!names(new_args) %in% p]
                                 serify("no ", pretty_params(no, TRUE), "")
                             }) |>
               factor() |>
               fct_relevel("all", after = Inf)) |>
    mutate(origin = ifelse(n_pseudo == 0, old_obs[["0"]], old_obs[["3"]])) |>
    multi_pert_plotter(c(1,9), c(-3.6, 3.9))

multi_pert_p2 <- multi_pert_sims |>
    filter(k == 0 | k == 1) |>
    mutate(params = map2_chr(k, params,
                             \(k, params) {
                                 if (k == 0) return("none")
                                 if (k == 1) return(serify("", pretty_params(params, TRUE),
                                                           ""))
                                 stop("k must be 0 or 1 for this fxn")
                             }) |>
               factor() |>
               fct_relevel("none", after = Inf)) |>
    mutate(origin = ifelse(n_pseudo == 0, old_obs[["0"]], old_obs[["3"]])) |>
    multi_pert_plotter(c(1,9), c(-3.6, 3.9))

multi_pert_p + multi_pert_p2 +
    plot_layout(widths = c(1, 1, 2))


multi_pert_sims |>
    # filter(k %in% c(1, 2)) |>
    # filter((k == 1 & !grepl("spat_config", params)) |
    #            (k == 2 & grepl("spat_config", params))) |>
    mutate(params = str_split(params, " & ") |>
               map_chr(\(x) paste(serify("", pretty_params(x, TRUE), ""),
                                  collapse = " & "))) |>
    mutate(dummy = paste(k, params, sep = "_"),
           params = fct_reorder(params, dummy, .fun = \(x) x[1], .desc = TRUE)) |>
    mutate(origin = ifelse(n_pseudo == 0, old_obs[["0"]], old_obs[["3"]])) |>
    multi_pert_plotter(c(1,9), c(-3.6, 3.9))


cbind(new = new_args, old = old_sims_args[names(new_args)])


# set.seed(262255805)
new_sims <- do.call(one_combo, list_assign(old_sims_args, !!!new_args))
new_sims0 <- do.call(one_combo, list_assign(old_sims_args, !!!c(new_args, list(n_pseudo = 0L))))


bind_rows(bind_rows(old_sims, old_sims0) |> mutate(params = "old"),
          bind_rows(new_sims, new_sims0) |> mutate(params = "new")) |>
    mutate(params = factor(params, levels = c("old", "new"))) |>
    filter(is.na(x)) |>
    select(n_pseudo, params, rep, time, virus:wasps) |>
    group_by(params, n_pseudo, rep) |>
    summarize(wasps_t = time[wasps == max(wasps)],
              max_wasps = max(wasps),
              outbreak_size = max(virus),
              .groups = "drop_last") |>
    summarize(across(wasps_t:last_col(), mean),
              .groups = "drop_last") |>
    (\(x) {print(ungroup(x)); return(x)})() |>
    summarize(across(wasps_t:last_col(), \(x) x[n_pseudo > 0] - x[n_pseudo == 0]))



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
                           # .tag = serify("",
                           #               pretty_params(names(new_arg), TRUE),
                           #               sprintf(" = %.2f",
                           #                       old_sims_args[names(new_arg)])) |>
                           #     paste(collapse = " | ") |>
                           #     paste("(original)"))
                           .tag = "original")
new_p <- paired_timeseries(new_sims, new_sims0,
                           .rep = .rep,
                           .aphid_max = sims_maxes$aphids,
                           .alate_max = sims_maxes$alates,
                           .wasp_max = sims_maxes$wasps,
                           # .tag = serify("",
                           #               pretty_params(names(new_arg), TRUE),
                           #               paste(" =", new_arg)) |>
                           #     paste(collapse = " | "))
                           .tag = "new")

old_p / new_p





# Nov-19 - LEFT OFF ----
old_sims_od <- do.call(one_combo, list_assign(old_sims_args, out_dispersals = TRUE,
                                              summ = "all",
                                              # n_pseudo = 0,
                                              n_sims = 1000))
new_sims_od <- do.call(one_combo, list_assign(old_sims_args,
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







# 2-par sims extras ----



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


paired_target_sims(extreme_combos[["high"]], .n_sims = 1000, .summ = "all",
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


sims <- paired_target_sims(extreme_combos[["high"]], zeta = 0.5)




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






