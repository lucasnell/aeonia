
library(tidyverse)
library(gameofclones)
library(aphidsync)  # calc_lambda
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


line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")

L0 <- line_s$leslie[,,1]
sad0 <- gameofclones:::sad_leslie(L0) |>
    (\(s) {
        out <- cbind(numeric(5))
        for (i in 1:4) out[i] <- sum(s[(2*i-1):(2*i)])
        out[5] <- sum(s[9:length(s)])
        return(out)
    })()


populations$surv_juv
populations$surv_adult$high[1:20]
populations$repro$high[1:30]
populations$s_y

L <- matrix(0, 5, 5)
L[1,5] <- mean(L0[1, 10:29])
# L[cbind(2:5, 1:4)] <- median(L0[row(L0) - col(L0) == 1][1:8])

# Use proportion that stay adults to match Leslie from gameofclones
lambda <- calc_lambda(line_s$leslie[,,1])
op <- nloptr::directL(\(x) {
    L1 <- L
    # L1[5,5] <- x[1]
    # L1[1,5] <- x[2]
    L1[cbind(2:5, 1:4)] <- x[1:4]
    L1[5,5] <- x[5]
    # ld <- abs(calc_lambda(L1) - lambda)
    sad1 <- gameofclones:::sad_leslie(L1)
    if (!identical(dim(sad1), dim(sad0))) return(100)
    sadd <- sum(abs(sad0 - sad1))
    return(sadd)
},
              lower = rep(0, 5), upper = rep(1, 5),
              control = list(stopval = 0, maxeval = 10e3,
                             xtol_rel = .Machine$double.eps^0.5))
op

# L[5,5] <- op[["par"]][[1]]
# L[1,5] <- op[["par"]][[2]]

L[cbind(2:5, 1:4)] <- op[["par"]][1:4]
L[5,5] <- op[["par"]][5]

Matrix::expm(L0)



cbind(calc_lambda(L), calc_lambda(L0))
cbind(gameofclones:::sad_leslie(L), sad0)




run_paras <- function(aphids0 = 500, paras0 = 1) {
    aphids <- numeric(max_t+1)
    paras <- numeric(max_t+1)
    X <- numeric(5)
    aphids[1] <- aphids0
    paras[1] <- paras0
    for (t in 1:max_t) {
        attack_surv <- (1 + a * paras[t] / (k * (h * aphids[t] + 1)))^(-k)
        aphids[t+1] = aphids[t] * exp(r - alpha * aphids[t]) * attack_surv
        paras[t+1] = paras[t] * exp(- pred_m) + aphids[t] * exp(r - alpha * aphids[t]) * (1 - attack_surv)

    }
    p <- tibble(time = 0:max_t, aphids = A, parasitoid = P) |>
        pivot_longer(-time, names_to = "species", values_to = "density") |>
        ggplot(aes(time, density, color = species)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(linewidth = 0.75) +
        scale_color_manual(values = c("dodgerblue", "firebrick")) +
        theme_classic()
    return(p)
}

a = 2.32
k = 0.35
h = 0.008

r = 0.2560
alpha = r / 1000
max_t = 1000



run_paras()



