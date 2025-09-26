
library(tidyverse)
library(gameofclones)
library(aeonia)



#' This only gets run if my local .Rprofile has been run and if it's an
#' interactive session:
if (interactive() && exists("LAN_USER")) {
    setHook(packageEvent("grDevices", "onLoad"),
            function(...) grDevices::quartz.options(width = 4, height = 4,
                                                    pointsize = 10))
    options("device" = "quartz")
    grDevices::graphics.off()
}


# Carrying capacity with no alates or natural enemies:
CC <- function(K, B, m){
    L  <- rbind(c(pop_info$surv_j, pop_info$fecund),
                c(pop_info$recruit, pop_info$surv_a))
    K * (max(abs(eigen(L)[["values"]])) * (1 - B) * (1 - m) - 1)
}


spp_pal <- viridisLite::plasma(101)[c(10, 50, 90)] |>
    set_names(c("aphids", "alates", "enemies"))


test_insect_pops(max_t = 10,
                 s = pop_info$s * 0,
                 B = 0,
                 A0 = 100,
                 W0 = 0,
                 P0 = 1) |>
    mutate(aphids = aphids + alates) |>
    select(-alates) |>
    pivot_longer(aphids:enemies, names_to = "spp", values_to = "N") |>
    mutate(spp = factor(spp, levels = c("aphids", "enemies"))) |>
    # filter(morph == "aphids") |>
    ggplot(aes(time, N, color = spp)) +
    # geom_hline(yintercept = CC(12500, 0, 0.1), color = "gray80", linewidth = 1) +
    geom_line(linewidth = 1) +
    # scale_color_viridis_d(begin = 0.2, end = 0.9, drop = FALSE) +
    scale_color_manual(NULL, values = spp_pal) +
    ylab("Abundance") +
    # coord_cartesian(ylim = c(0, 1600)) +
    theme()






insect_sims <- map(1:10, \(i){
    test_insect_pops(max_t = 100,
                     B = 0,
                     h = 5,
                     alate_0 = -Inf,
                     alate_1 = 0,
                     demog_error = TRUE,
                     disaster_p = 0.02,
                     disaster_s = 0.1,
                     A0 = 10,
                     W0 = 0,
                     P0 = 0) |>
        mutate(aphids = aphids + alates) |>
        select(-alates) |>
        pivot_longer(aphids:enemies, names_to = "species", values_to = "N") |>
        mutate(species = factor(species, levels = c("aphids", "enemies")),
               rep = i)
}) |>
    list_rbind() |>
    mutate(rep = factor(rep))
insect_sims2 <- test_insect_pops(max_t = 100,
                                 B = 0,
                                 h = 5,
                                 alate_0 = -Inf,
                                 alate_1 = 0,
                                 demog_error = FALSE,
                                 A0 = 10,
                                 W0 = 0,
                                 P0 = 0) |>
    mutate(aphids = aphids + alates) |>
    select(-alates) |>
    pivot_longer(aphids:enemies, names_to = "species", values_to = "N") |>
    mutate(species = factor(species, levels = c("aphids", "enemies")))

# insect_time_pts <- c(375, 434, 500, 577)
# insect_sims |>
#     filter(time %in% insect_time_pts)

insect_sims |>
    # filter(species != "alates") |>
    ggplot(aes(time, N, color = species)) +
    geom_hline(yintercept = CC(12500, 0, 0.1), color = "gray80", linewidth = 1) +
    geom_line(aes(group = interaction(rep, species)), linewidth = 1, alpha = 0.25) +
    geom_line(data = insect_sims2 |> filter(species == "aphids"), linewidth = 1,
              color = "red") +
    # geom_vline(xintercept = insect_time_pts, linetype = "22",
    #            color = "gray70", linewidth = 1) +
    scale_color_viridis_d(begin = 0.2, end = 0.9) +
    scale_y_continuous("Abundance") +
    theme_minimal()

















# ===========================================================================*
# Predator version ----
# ===========================================================================*


run_pred <- function() {
    A <- numeric(max_t+1)
    P <- numeric(max_t+1)
    A[1] <- 100
    P[1] <- 10
    for (t in 1:max_t) {
        # predator per-capita consumption per prey:
        # pcc = pred_r
        pcc = A[t]^(pred_k-1) * pred_r / (1 + pred_r * pred_h * A[t]^pred_k)
        # A[t+1] = A[t] * exp(r - alpha * A[t] - P[t]^pred_k)
        A[t+1] = A[t] * exp(r - alpha * A[t] - pcc * P[t])
        # P[t+1] = pred_c * exp(pred_0 + pred_1 * log10(A[t]))
        P[t+1] = P[t] * exp(pred_cc * pcc * A[t] - pred_m)
    }
    p <- tibble(time = 0:max_t, aphids = A, predator = P) |>
        pivot_longer(-time, names_to = "species", values_to = "density") |>
        ggplot(aes(time, density, color = species)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(linewidth = 0.75) +
        scale_color_manual(values = c("dodgerblue", "firebrick")) +
        theme_classic()
    return(p)
}



r = 0.2560
alpha = r / 1000

pred_c = 0.6
pred_0 = -1.6181
pred_1 = 0.4374
max_t = 1000

pred_cc = 0.5
pred_m = 0.1


# pred_h = 0
# pred_k = 2
# pred_r = 8e-6

# pred_h = 0
# pred_k = 1
# pred_r = 4e-3

pred_h = 0.1
pred_k = 1
pred_r = 5e-3


run_pred()



pred_h = 0.02
pred_k = 1
pred_r = 2e-3
max_t = 1000
run_pred()





# ===========================================================================*
# Parasitoid version ----
# ===========================================================================*


run_paras <- function(aphids0 = 500, wasps0 = 1,
                      a = 2.32, k = 0.35, h = 0.008, r = 0.2560,
                      alpha = 1 / 2400, s_y = 0.69,
                      s_p = 0.9745, r_m = 1/10,
                      max_t = 1000, mummies = TRUE) {
    .spp_pal <- c(aphids = "dodgerblue", wasps = "firebrick",
                  parasitized = "goldenrod")
    aphids <- numeric(max_t+1)
    paras <- numeric(max_t+1)
    wasps <- numeric(max_t+1)
    aphids[1] <- aphids0
    wasps[1] <- wasps0
    if (!mummies) {
        for (t in 1:max_t) {
            attack_surv <- (1 + a * wasps[t] / (k * (h * aphids[t] + 1)))^(-k)
            aphids[t+1] <- aphids[t] * exp(r * (1 - alpha * aphids[t])) * attack_surv
            wasps[t+1] <- wasps[t] * s_y + aphids[t] *
                exp(r * (1 - alpha * aphids[t])) * (1 - attack_surv)
        }
        d <- tibble(time = 0:max_t, aphids = aphids, wasps = wasps)
    } else {
        for (t in 1:max_t) {
            attack_surv <- (1 + a * wasps[t] / (k * (h * aphids[t] + 1)))^(-k)
            aphids[t+1] <- aphids[t] * exp(r * (1 - alpha * aphids[t])) * attack_surv
            paras[t+1] <- paras[t] * s_p * (1 - r_m) + aphids[t] *
                exp(r * (1 - alpha * aphids[t])) * (1 - attack_surv)
            wasps[t+1] <- wasps[t] * s_y + paras[t] * s_p * r_m
        }
        d <- tibble(time = 0:max_t, aphids = aphids, parasitized = paras, wasps = wasps)
    }

    p <- d |>
        pivot_longer(-time, names_to = "species", values_to = "density") |>
        ggplot(aes(time, density, color = species)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(linewidth = 0.75) +
        scale_color_manual(values = .spp_pal)

    return(p)
}



run_paras(10, 1, mummies = FALSE)
run_paras(10, 1)


library(gameofclones)

line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,10), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")

sim_experiments(line_s, 1, 1000, wasp_density_0 = 1, wasp_delay = 0, extinct_N = 0) |>
    (\(x) {
        x$aphids |>
            mutate(species = ifelse(type %in% c("mummy", "parasitized"),
                                "parasitized", "aphids")) |>
            group_by(time, species) |>
            summarize(density = sum(N), .groups = "drop") |>
            bind_rows(x$wasps |> select(time, wasps) |>
                          mutate(species = "wasps") |> rename(density = wasps))
    })() |>
    ggplot(aes(time, density, color = species)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(linewidth = 0.75) +
        scale_color_manual(values = c(aphids = "dodgerblue", wasps = "firebrick",
                                      parasitized = "goldenrod"))





