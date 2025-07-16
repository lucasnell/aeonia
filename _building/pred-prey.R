
library(tidyverse)
library(gameofclones)
library(aeonia)


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


run_paras <- function(aphids0 = 500, paras0 = 1,
                      a = 2.32, k = 0.35, h = 0.008, r = 0.2560,
                      alpha = 1 / 2400, s_y = 0.69, max_t = 1000) {
    aphids <- numeric(max_t+1)
    paras <- numeric(max_t+1)
    aphids[1] <- aphids0
    paras[1] <- paras0
    for (t in 1:max_t) {
        attack_surv <- (1 + a * paras[t] / (k * (h * aphids[t] + 1)))^(-k)
        aphids[t+1] = aphids[t] * exp(r * (1 - alpha * aphids[t])) * attack_surv
        paras[t+1] = paras[t] * s_y + aphids[t] * exp(r * (1 - alpha * aphids[t])) * (1 - attack_surv)

    }
    p <- tibble(time = 0:max_t, aphids = aphids, parasitoid = paras) |>
        pivot_longer(-time, names_to = "species", values_to = "density") |>
        ggplot(aes(time, density, color = species)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(linewidth = 0.75) +
        scale_color_manual(values = c("dodgerblue", "firebrick")) +
        theme_classic()
    return(p)
}



run_paras(10, 1, h = 0.008, a = 0.01)
