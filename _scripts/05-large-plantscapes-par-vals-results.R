
#'
#' Find par values for larger landscape simulations
#' 04-large-plantscapes-par-vals.sh must be run on the cluster first, then
#' its output sent to ./_scripts/interm-data/
#'


source("_scripts/03-large-preamble.R")

test_sims <- read_rds("_scripts/interm-data/large-plantscapes-par-vals.rds") |>
    # determine which two scenario: Pseudomonas promotes or inhibits viruses, resp.:
    mutate(scenario = factor(zeta < 0.2 & N0 < 20, levels = c(TRUE, FALSE),
                             labels = c("promotes", "inhibits"))) |>
    select(scenario, everything()) |>
    mutate(outbreak_size = log10(outbreak_size)) |>
    pivot_longer(outbreak_size:p_emerge, names_to = "outcome") |>
    filter(!is.na(value)) |>
    mutate(outcome = factor(outcome,
                            levels = c("outbreak_size", "p_emerge"),
                            labels = c("log10(outbreak size)", "prob. emerge")),
           n_pseudo = factor(n_pseudo))


np_pal <- c(`0` = "#999999", `5000` = "#0046D2")



# Based on below, I'll choose p_load = 0.075
test_sims |>
    ggplot(aes(p_load, value, color = n_pseudo, shape = scenario)) +
    geom_hline(yintercept = 0) +
    geom_point(position = position_jitterdodge(jitter.width = 0.005, dodge.width = 0.01)) +
    facet_wrap(~ outcome, ncol = 1, scales = "free") +
    scale_color_manual(values = np_pal) +
    scale_shape_manual(values = c(1, 2))



p_load <- 0.075

#'
#' From below graph and table, I should use
#' high zeta (0.12) when *Pseudomonas* promotes outbreaks, and
#' medium zeta (0.95) when *Pseudomonas* inhibits outbreaks.
#' Because I simulated 5000 *Pseudomonas* patches instead of 7000 (the density that
#' typically produces the strongest *Pseudomonas* effect),
#' I will use zeta = 0.14 when *Pseudomonas* promotes outbreaks and
#' zeta = 1.0 when *Pseudomonas* inhibits outbreaks.
#'


test_sims |>
    filter(p_load == .env$p_load) |>
    select(-p_load) |>
    filter(Y0 == 100) |>
    filter(N0 %in% c(10, 110)) |>
    group_by(scenario) |>
    mutate(zeta_lvl = case_when(zeta == min(zeta) ~ 1L,
                                zeta == max(zeta) ~ 3L,
                                .default = 2L)) |>
    ungroup() |>
    ggplot(aes(zeta_lvl, value, color = n_pseudo)) +
    # geom_hline(yintercept = 0) +
    geom_point(aes(shape = scenario)) +
    geom_line(aes(linetype = scenario)) +
    scale_x_continuous(breaks = 1:3, labels = c("low", "mid", "high")) +
    facet_grid(outcome ~ ., scales = "free") +
    scale_color_manual(values = np_pal)


test_sims |>
    filter(p_load == .env$p_load) |>
    filter(Y0 == 100) |>
    filter(N0 %in% c(10, 110)) |>
    select(-p_load, -Y0, -N0) |>
    group_by(scenario) |>
    mutate(zeta_lvl = case_when(zeta == min(zeta) ~ 1L,
                            zeta == max(zeta) ~ 3L,
                            .default = 2L)) |>
    ungroup() |>
    filter((scenario == "promotes" & zeta_lvl == 3L) |
               (scenario == "inhibits" & zeta_lvl == 2L)) |>
    arrange(outcome, scenario, n_pseudo) |>
    select(outcome, scenario, n_pseudo, value)



