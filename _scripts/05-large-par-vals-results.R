
#'
#' Find par values for larger landscape simulations
#' 04-large-par-vals.sh must be run on the cluster first, then
#' its output sent to ./_scripts/interm-data/
#'


source("_scripts/03-large-preamble.R")
library(patchwork)
library(ggtext)
library(viridisLite)


# ============================================================================*
# Small outbreaks ----
# ============================================================================*


test_sims <- list.files("_scripts/interm-data", "large-par-vals-.?.?.rds",
                        full.names = TRUE) |>
    map(read_rds) |>
    list_rbind() |>
    # Because this doesn't vary:
    select(-N0) |>
    pivot_longer(outbreak_size:n_infected, names_to = "outcome") |>
    # filter(!is.na(value)) |>
    mutate(outcome = factor(outcome, levels = c("outbreak_size", "p_emerge",
                                                "n_infected"),
                            labels = c("outbreak size", "prob. emerge",
                                       "infected plants")),
           n_pseudo = factor(n_pseudo)) |>
    arrange(p_load, Y0, outcome, n_pseudo)


np_pal <- plasma(10, begin = 0.1, end = 0.9, direction = -1L)[c(1, 7)] |>
    set_names(levels(test_sims$n_pseudo))

# Potential value to use for Y0:
Y0_line <- 250


test_plots <- levels(test_sims$outcome) |>
    set_names() |>
    map(\(yvar) {
        sort(unique(test_sims$p_load)) |>
            map(\(pl) {
                # rm(pl, dd, yvar, p1, p2, y_line)
                dd <- test_sims |>
                    filter(p_load == pl, outcome == yvar) |>
                    mutate(zeta_fct =
                               factor(zeta, labels = sprintf(
                                   "&zeta; = %.1f", sort(unique(test_sims$zeta)))))
                if (str_detect(yvar, "prob")) {
                    y_line <- 0
                    y_lim <- NULL
                } else if (str_detect(yvar, "size")) {
                    y_line <- 2
                    y_lim <- c(-0.3, NA)
                } else {
                    y_line <- 1
                    y_lim <- c(-0.3, NA)
                }
                p1 <- dd |>
                    group_by(zeta, zeta_fct, Y0) |>
                    summarize(value = value[n_pseudo != "0"] - value[n_pseudo == "0"],
                              .groups = "drop") |>
                    filter(!is.na(value)) |>
                    mutate(value = ifelse(zeta > 0.5, - value, value)) |>
                    ggplot(aes(Y0, value, linetype = zeta_fct)) +
                    geom_hline(yintercept = 0, color = "gray70") +
                    geom_vline(xintercept = Y0_line, color = "gray70") +
                    geom_line(aes(linewidth = zeta_fct)) +
                    labs(x = "Y<sub>0</sub>",
                         y = paste("Standardized difference<br>in", yvar),
                         title = paste("&delta;<sub>p</sub> = &delta;<sub>a</sub> =", pl)) +
                    coord_cartesian(ylim = y_lim) +
                    scale_linetype_manual(NULL, values = c("solid", "22")) +
                    scale_linewidth_manual(NULL, values = c(0.75, 1))
                p2 <- dd |>
                    filter(!is.na(value)) |>
                    ggplot(aes(Y0, value, color = n_pseudo, linetype = zeta_fct)) +
                    geom_hline(yintercept = y_line, color = "gray70") +
                    geom_vline(xintercept = Y0_line, color = "gray70") +
                    geom_line(aes(linewidth = zeta_fct)) +
                    scale_x_continuous(breaks = c(100, 300, 500)) +
                    labs(x = "Y<sub>0</sub>", y = str_to_sentence(yvar)) +
                    # facet_wrap(~ zeta_fct, nrow = 1) +
                    scale_color_manual("*Pseudo.*<br>plants", values = np_pal) +
                    scale_linetype_manual(NULL, values = c("solid", "22")) +
                    scale_linewidth_manual(NULL, values = c(0.75, 1))
                return(list(p1, p2))
            }) |>
            do.call(what = c) |>
            wrap_plots(nrow = 2, guides = "collect", axis_titles = "collect",
                       byrow = FALSE) +
            plot_annotation(title = yvar |>
                                str_replace("prob.", "probability of") |>
                                str_replace("emerge", "emergence") |>
                                str_to_sentence(),
                            theme = theme(plot.title = element_markdown(face = "bold")))
    })




test_plots$`prob. emerge`
# test_plots$`outbreak size`
test_plots$`infected plants`

