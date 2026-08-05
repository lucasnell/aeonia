
#'
#' Empirical zeta from large landscape simulations.
#'

source("_scripts/00-preamble.R")

.overwrite <- FALSE







# In Ives et al. (1999), they found that parasitoids spent ~3.76 more time
# foraging at plants where they encountered an aphid.
# Bestrong simulates our model with varying zeta, then compares the observed
# wasp abundances to those predicted when parasitoids never encounter an aphid
# on a Pseudomonas-inhabited plant but always do on plants without Pseudomonas.

if (!file.exists(interm_files$large_zeta_sims)) stop("Run 02-large-hpc/04-large-emp-zeta.R first!")

# Takes a little while to read
zeta_sims <- read_csv(interm_files$large_zeta_sims, col_types = "diicidddddld")

add_n_pseudo_fct <- function(d, label_both = FALSE) {
    fmt <- ifelse(label_both, "%.0f%% *Pseudo.*", "%.0f%%")
    n_pseudo_labs <- n_pseudo_lvls |>
        (\(x) x / 10e3 * 100)() |>
        (\(x) sprintf(fmt, x))()
    d |>
        mutate(n_pseudo_fct = factor(n_pseudo, levels = n_pseudo_lvls,
                                     labels = n_pseudo_labs))
}

# Takes ~ 2 min
set.seed(1540192361)
zeta_sim_summs <- zeta_sims |>
    mutate(rel = wasps / pred_wasps) |>
    group_by(zeta, n_pseudo, pseudo) |>
    summarize(ci = list(booter(rel)[c("Lower", "Upper")] |>
                            as.list() |> as_tibble()),
              min = min(rel), max = max(rel),
              rel = mean(rel), .groups = "drop") |>
    unnest(ci) |>
    mutate(pseudo = factor(pseudo, levels = c(FALSE, TRUE),
                           labels = c("Plants without<br>*Pseudomonas*",
                                      "Plants with<br>*Pseudomonas*"))) |>
    add_n_pseudo_fct(label_both = TRUE)


# approximate using linear interpolation where lines equal 1:
zeta_one_ests <- zeta_sim_summs |>
    filter(n_pseudo > 0) |>
    group_by(n_pseudo_fct, n_pseudo, pseudo) |>
    summarize(zeta = (\(z, r) {
        # Linear interpolation:
        f <- approxfun(z, r)
        # Find where it equals 1:
        approx_z <- uniroot(function(x) f(x) - 1, range(z))$root
        return(approx_z)
    })(zeta, rel), .groups = "drop_last") |>
    summarize(zeta_d = abs(diff(range(zeta))),
              zeta = mean(zeta), .groups = "drop")

# They are the same across whether Pseudomonas is on plant:
zeta_one_ests$zeta_d |> max()
# [1] 1.332268e-15

zeta_one_ests |> select(n_pseudo, zeta)
#   n_pseudo  zeta
#      <int> <dbl>
# 1     1000 0.749
# 2     3000 0.698
# 3     5000 0.625
# 4     7000 0.509
# 5     9000 0.296


zeta_p <- zeta_sim_summs |>
    filter(n_pseudo > 0) |>
    ggplot(aes(zeta, rel, color = pseudo)) +
    geom_hline(yintercept = 1, color = "gray70", linewidth = 1) +
    geom_vline(data = zeta_one_ests, aes(xintercept = zeta),
               color = "black", linetype = "22") +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = pseudo),
                # geom_ribbon(aes(ymin = min, ymax = max, fill = pseudo),
                color = NA, alpha = 0.25) +
    geom_line(linewidth = 0.75) +
    facet_wrap(~ n_pseudo_fct, nrow = 2, axes = "all", axis.labels = "margins") +
    labs(x = pretty_params("zeta", cap1 = TRUE),
         y = "Observed / predicted parasitoid density") +
    scale_color_manual(NULL, values = np_pal |> set_names(nm = NULL),
                       aesthetics = c("color", "fill")) +
    theme(legend.position = c(5/6, 1/4),
          legend.justification = c(0.5, 0.5),
          legend.key.spacing.y = unit(1, "lines"))

# zeta_p

if (.overwrite) {
    save_plot("_plots/large-empirical-zeta.pdf", zeta_p, width = 6, height = 5)
}



# Now plot for zeta ~ n_infected plants
set.seed(19813000)
zeta_ninf_sim_summs <- zeta_sims |>
    filter(n_pseudo > 0) |>
    group_by(zeta, n_pseudo, rep) |>
    summarize(n_infected = sum(max_virus), .groups = "drop_last") |>
    summarize(ci = list(booter(n_infected)[c("Lower", "Upper")] |>
                            as.list() |> as_tibble()),
              n_infected = mean(n_infected), .groups = "drop") |>
    unnest(ci)

set.seed(694769807)
no_pseudo_df <- zeta_sims |>
    filter(n_pseudo == 0) |>
    group_by(zeta, rep) |>
    summarize(n_infected = sum(max_virus), .groups = "drop") |>
    summarize(ci = list(booter(n_infected)[c("Lower", "Upper")] |>
                            as.list() |> as_tibble()),
              n_infected = mean(n_infected), .groups = "drop") |>
    unnest(ci)





zeta_ninf_p <- zeta_ninf_sim_summs |>
    filter(n_pseudo > 0) |>
    add_n_pseudo_fct() |>
    ggplot(aes(zeta, n_infected, color = n_pseudo_fct)) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = n_pseudo_fct), alpha = 0.25, color = NA) +
    geom_line() +
    geom_hline(yintercept = no_pseudo_df$n_infected) +
    geom_hline(yintercept = as.numeric(no_pseudo_df[,c("Lower", "Upper")]),
               linetype = "22") +
    geom_point(data = zeta_one_ests |>
                   # Use linear interpolation to estimate value of n_infected for each zeta:
                   mutate(n_infected = map2_dbl(n_pseudo, zeta, \(np, z) {
                       d <- zeta_ninf_sim_summs |> filter(n_pseudo == np)
                       f <- approxfun(d$zeta, d$n_infected)
                       return(f(z))
                   })) |>
                   add_n_pseudo_fct(),
               size = 3) +
    scale_color_manual("*Pseudo.*<br>plants",
                       values = full_np_pal[paste(n_pseudo_lvls)] |>
                           set_names(sprintf("%.0f%%", n_pseudo_lvls / 10e3 * 100)),
                       aesthetics = c("color", "fill")) +
    labs(x = pretty_params("zeta", cap1 = TRUE),
         y = paste0("Peak infected plants",
                    "<br><span style='font-size:8pt;'>",
                    "(mean across 100 sims)</span>"))


# zeta_ninf_p



if (.overwrite) {
    save_plot("_plots/large-zeta-n_infected.pdf", zeta_ninf_p, width = 6, height = 4)
}
