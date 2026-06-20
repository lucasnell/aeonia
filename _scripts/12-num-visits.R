
library(gtools)
source("_scripts/00-preamble.R")

.overwrite <- FALSE

# # Not used:
# #
# # Probability of at least one infected plant visited:
# p_inf_visit <- function(p_infected, n_visits) {
#     1 - (1 - p_infected)^n_visits
# }
# # Probability of at least one UNinfected plant visited:
# p_uninf_visit <- function(p_infected, n_visits) {
#     1 - p_infected^n_visits
# }


# Expected number of events where vector visits infected then uninfected plant.
# Vectorized across `p_infected` values.
exp_inf2uninf <- function(p_infected, n_visits, virus_attract = 1) {

    # p_infected = c(0.1, 0.5, 0.9); n_visits = 5; virus_attract = 1
    # rm(p_infected, n_visits, virus_attract)

    stopifnot(is.numeric(p_infected) && length(p_infected) >= 1)
    stopifnot(is.numeric(n_visits) && length(n_visits) == 1)
    stopifnot(is.numeric(virus_attract) && length(virus_attract) == 1)

    # Possible permutations of `n_visits+1` visits of uninfected (`0`) and
    # infected (`1`) plants
    # Using `+1` to include the starting plant
    perms <- permutations(n = 2, r = n_visits + 1L, v = 0:1,
                          repeats.allowed = TRUE)
    # For each permutation, how many times did the alate go from infected
    # to uninfected?
    n_events <- sapply(1:nrow(perms), \(i) sum(diff(perms[i,]) == -1))

    # Scaling each for the probability that they would occur (based on
    # `p_infected` and `virus_attract`)
    out <- numeric(length(p_infected))
    for (j in 1:length(p_infected)) {
        p_inf <- p_infected[j]
        if (virus_attract != 1) {# see (Donnelly et al. 2019)
            p_inf <- virus_attract * p_inf / ((1 - p_inf) + virus_attract * p_inf)
        }
        for (i in 1:nrow(perms)) {
            p_i <- (1 - p_inf)^sum(perms[i,] == 0) * p_inf^sum(perms[i,] == 1)
            out[j] <- out[j] + n_events[i] * p_i
        }
    }

    return(out)
}

n_visits_df <- crossing(n_visits = 1:6, virus_attract = c(1, 5, 100)) |>
    mutate(n_events = map2(n_visits, virus_attract, \(nv, va) {
        p_i <- 0:100 / 100
        ev <- exp_inf2uninf(p_i, nv, virus_attract = va)
        return(tibble(p_infected = p_i, n_events = ev))
    })) |>
    unnest(n_events)



n_visits_p <- n_visits_df |>
    mutate(n_visits = factor(n_visits),
           virus_attract = factor(virus_attract, levels = sort(unique(virus_attract)),
                                  labels = sort(unique(virus_attract)) |>
                                      sprintf(fmt = "&nu; = %.0f") |>
                                      serify(prefix = "", suffix = ""))) |>
    ggplot(aes(p_infected, n_events, color = n_visits)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(linewidth = 1) +
    scale_color_viridis_d("New<br>plants<br>visited", begin = 0.1, end = 0.9) +
    facet_wrap(~ virus_attract) +
    scale_x_continuous(breaks = 0:2 / 2) +
    labs(x = "Proportion plants infected",
         y = "Expected infected &rarr; uninfected<br>vector visits per flight") +
    theme(strip.text.x = element_markdown(size = 14))

# n_visits_p


n_visits_df |>
    mutate(n_visits = factor(n_visits, levels = sort(unique(n_visits)),
                             labels = sort(unique(n_visits)) |>
                                 sprintf(fmt = "%i new plants<br>visited")),
           virus_attract = factor(virus_attract)) |>
    ggplot(aes(p_infected, n_events, color = virus_attract)) +
    geom_hline(yintercept = 0, color = "gray70") +
    geom_line(linewidth = 1) +
    scale_color_viridis_d(pretty_params("virus_attract", FALSE, TRUE, TRUE) |>
                              # having 2 of these is intentional:
                              str_replace(" ", "<br>") |>
                              str_replace(" ", "<br>"),
                          begin = 0.1, end = 0.9) +
    facet_wrap(~ n_visits, scales = "free_y") +
    scale_x_continuous(breaks = 0:2 / 2) +
    labs(x = "Proportion plants infected",
         y = "Expected infected &rarr; uninfected<br>vector visits per flight")



if (.overwrite) {
    save_plot("_plots/infected-to-uninfected.pdf", n_visits_p, width = 6.5, height = 3)
}




