
#'
#' Find par values for larger landscape simulations
#' 04-large-plantscapes-par-vals.sh must be run on the cluster first, then
#' its output sent to ./_scripts/interm-data/
#'


source("_scripts/03-large-preamble.R")



# ============================================================================*
# Small outbreaks ----
# ============================================================================*


test_sims <- list.files("_scripts/interm-data", "large-plantscapes-par-vals-.?.?.rds",
                        full.names = TRUE) |>
    map(read_rds) |>
    list_rbind() |>
    mutate(scenario = interaction(Y0, N0, zeta, drop = TRUE)) |>
    select(scenario, p_load, everything()) |>
    pivot_longer(outbreak_size:p_emerge, names_to = "outcome") |>
    filter(!is.na(value)) |>
    mutate(outcome = factor(outcome,
                            levels = c("outbreak_size", "p_emerge"),
                            labels = c("outbreak size", "prob. emerge")),
           n_pseudo = factor(n_pseudo)) |>
    arrange(scenario, p_load, outcome, n_pseudo) |>
    # Because this doesn't vary:
    select(-N0)


np_pal <- c("#999999", "#0046D2") |>
    set_names(levels(test_sims$n_pseudo))



# Based on below, I'll choose p_load = 0.075
test_sims |>
    # group_by(p_load, scenario) |>
    # mutate(within = max(value[outcome == "prob. emerge"]) < 1 &
    #            min(value[outcome == "prob. emerge"]) > 0) |>
    # ungroup() |>
    # filter(within) |>
    # filter(zeta == min(zeta)) |>
    # filter(zeta == 0.1, p_load == 0.15) |>
    # filter(p_load == 0.1) |>
    filter(outcome == "outbreak size") |>
    # filter(outcome == "prob. emerge") |>
    group_by(p_load, scenario) |>
    summarize(diff = value[n_pseudo != "0"] / value[n_pseudo == "0"], .groups = "drop") |>
    # filter(diff > 0) |>
    arrange(desc(abs(diff)))



# For outbreak size and ascending sort:
#    p_load scenario   diff
#     <dbl> <fct>     <dbl>
#  1   0.15 275.55.1 0.0613
#  2   0.1  225.55.1 0.0617
#  3   0.15 300.55.1 0.0629
#  4   0.15 350.55.1 0.0640
#  5   0.1  175.55.1 0.0643
#  6   0.1  200.55.1 0.0656
#  7   0.15 375.55.1 0.0660
#  8   0.1  250.55.1 0.0665
#  9   0.15 250.55.1 0.0680
# 10   0.1  275.55.1 0.0687
#
# For outbreak size and ascending sort:
#    p_load scenario     diff
#     <dbl> <fct>       <dbl>
#  1   0.1  400.55.0.1   4.88
#  2   0.1  375.55.0.1   4.32
#  3   0.1  350.55.0.1   4.07
#  4   0.15 400.55.0.1   4.06
#  5   0.15 375.55.0.1   3.93
#  6   0.1  325.55.0.1   3.75
#  7   0.15 400.55.0.15  3.40
#  8   0.15 325.55.0.1   3.22
#  9   0.1  400.55.0.15  3.22
# 10   0.15 350.55.0.1   3.17
#
#
#
# For prob. emerge and ascending sort:
#    p_load scenario    diff
#     <dbl> <fct>      <dbl>
#  1   0.05 400.55.1   0.322
#  2   0.05 375.55.1   0.404
#  3   0.05 375.55.0.9 0.444
#  4   0.05 350.55.1   0.489
#  5   0.05 300.55.1   0.541
#  6   0.05 400.55.0.9 0.560
#  7   0.05 275.55.1   0.573
#  8   0.05 350.55.0.9 0.580
#  9   0.05 325.55.0.9 0.594
# 10   0.05 225.55.1   0.6
#
# For prob. emerge and descending sort:
#    p_load scenario     diff
#     <dbl> <fct>       <dbl>
#  1   0.05 400.55.0.1   1.17
#  2   0.05 400.55.0.15  1.14
#  3   0.05 350.55.0.15  1.11
#  4   0.05 325.55.0.1   1.08
#  5   0.05 375.55.0.15  1.07
#  6   0.05 375.55.0.1   1.05
#  7   0.05 300.55.0.15  1.05
#  8   0.05 350.55.0.3   1.04
#  9   0.05 325.55.0.15  1.03
# 10   0.05 275.55.0.15  1.03



test_sims |>
    filter(p_load == 0.05, Y0 == 400, zeta %in% c(0.1, 0.9)) |>
    select(outcome, zeta, n_pseudo, value) |>
    arrange(outcome, zeta, n_pseudo)

#   outcome        zeta n_pseudo value
#   <fct>         <dbl> <fct>    <dbl>
# 1 outbreak size   0.1 0         6.15
# 2 outbreak size   0.1 7000     14.7
# 3 outbreak size   0.9 0         6.43
# 4 outbreak size   0.9 7000      2.49
# 5 prob. emerge    0.1 0         0.84
# 6 prob. emerge    0.1 7000      0.98
# 7 prob. emerge    0.9 0         0.91
# 8 prob. emerge    0.9 7000      0.51




# _ conclusions ----
# For small outbreaks:
# - Keep N0 at 55
# - p_load = 0.05
# - Y0 = 400
# - zeta = 0.1 or 0.9




# ============================================================================*
# Big outbreaks ----
# ============================================================================*

#
# This section is testing for simulating large outbreaks (where p_emerge = 1
# and outbreak_size approaches # plants)
#


# If running below on the cluster (takes much less time):
#
# cd /home2/lan68/04-large-plantscapes-par-vals/
# srun -N 1 -n 1 -c 50 --mem=100G --time=1-20:00:00 --job-name="aeonia-test" --pty R --vanilla
#
# source("../03-large-preamble.R")


landscape1 <- sim_df |>
    filter(n_pseudo == 7000, wt_vp == 1e-6, wt_pp == 1) |>
    getElement("landscape") |> getElement(1)
landscape0 <- array(c(1L, rep(0L, 99999L)), c(100L, 100L, 1L))



one_test <- function(p_load, Y0, N0, zeta = c(0.1, 0.9)) {

    large_simmer <- function(landscape, Y0, N0, zeta, p_load) {

        args <- list(landscape = landscape,
                     sd_N = 0,
                     virus_attract = 5,
                     pseudo_repel = 1,
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

    np <- sum(landscape1[,,1] > 1L)

    one_pair <- function(zeta) {
        n_inf1 <- large_simmer(landscape = landscape1,
                               Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
            getElement("n_infected")
        n_inf0 <- large_simmer(landscape = landscape0,
                               Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
            getElement("n_infected")

        tibble(# p_load = .env$p_load, Y0 = .env$Y0, N0 = .env$N0,
               zeta = .env$zeta,
               n_pseudo = c(np, 0L),
               outbreak_size = c(mean(n_inf1[n_inf1 > 1]), mean(n_inf0[n_inf0 > 1])),
               p_emerge = c(mean(n_inf1 > 1), mean(n_inf0 > 1)))
    }

    map(zeta, one_pair) |>
        list_rbind()

}


# > (sim <- one_test(p_load = 0.5, Y0 = 400, N0 = 55))
# # A tibble: 4 × 4
#    zeta n_pseudo outbreak_size p_emerge
#   <dbl>    <int>         <dbl>    <dbl>
# 1   0.1     7000         5567.        1
# 2   0.1        0         1835.        1
# 3   0.9     7000          397.        1
# 4   0.9        0         1914.        1

# > (sim <- one_test(p_load = 0.5, Y0 = 300, N0 = 55))
# # A tibble: 4 × 4
#    zeta n_pseudo outbreak_size p_emerge
#   <dbl>    <int>         <dbl>    <dbl>
# 1   0.1     7000         6585.        1
# 2   0.1        0         6612.        1
# 3   0.9     7000          662.        1
# 4   0.9        0         3202.        1

# > (sim <- one_test(p_load = 0.5, Y0 = 400, N0 = 65))
# # A tibble: 4 × 4
#    zeta n_pseudo outbreak_size p_emerge
#   <dbl>    <int>         <dbl>    <dbl>
# 1   0.1     7000         5086.        1
# 2   0.1        0         2312.        1
# 3   0.9     7000          531.        1
# 4   0.9        0         2375.        1

# > (sim <- one_test(p_load = 1, Y0 = 400, N0 = 55))
# # A tibble: 4 × 4
#    zeta n_pseudo outbreak_size p_emerge
#   <dbl>    <int>         <dbl>    <dbl>
# 1   0.1     7000         8824.        1
# 2   0.1        0         3873.        1
# 3   0.9     7000         1271.        1
# 4   0.9        0         4016.        1

# > (sim <- one_test(p_load = 1, Y0 = 200, N0 = 55))
# # A tibble: 4 × 4
#    zeta n_pseudo outbreak_size p_emerge
#   <dbl>    <int>         <dbl>    <dbl>
# 1   0.1     7000         9854.        1
# 2   0.1        0         8556.        1
# 3   0.9     7000         3771.        1
# 4   0.9        0         8438.        1

# > (sim <- one_test(p_load = 1, Y0 = 100, N0 = 55))
# # A tibble: 4 × 4
#    zeta n_pseudo outbreak_size p_emerge
#   <dbl>    <int>         <dbl>    <dbl>
# 1   0.1     7000         9987.        1
# 2   0.1        0         9841.        1
# 3   0.9     7000         7053.        1
# 4   0.9        0         9853.        1


# _ conclusions ----
# For big outbreaks:
# - Keep N0 at 55
# - Always use p_load = 1
# - When zeta is low: Y0 = 400
# - When zeta is high: Y0 = 200
