library(tidyverse)
library(gameofclones)
library(aphidsync)  # calc_lambda
library(aeonia)



line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high",
                      p_instar_smooth = 0)

L0 <- line_s$leslie[,,1]


X0 <- gameofclones:::sad_leslie(L0)

# Juveniles staying juveniles:
# Note: indices in `X0[1:7]` and `X0[1:8]` below are intentional!
sj <- sum(L0[row(L0) - col(L0) == 1][1:7] * (X0[1:7] / sum(X0[1:8])))

# Recruitment: juveniles transitioning to adults:
r <- L0[row(L0) - col(L0) == 1][8] * (X0[8] / sum(X0[1:8]))

# Adults staying adults:
sa <- sum(L0[row(L0) - col(L0) == 1][9:28] * (X0[9:28] / sum(X0[9:28])))

# Fecundity:
f <- sum(L0[1,9:29] * (X0[9:29] / sum(X0[9:29])))

# Create simplified transition matrix:
L <- rbind(c(sj, f),
           c(r, sa))

# Are growth rates the same?
cbind(calc_lambda(L0), calc_lambda(L))

# Are stable age distributions the same?
cbind(orig = c(sum(X0[1:8]), sum(X0[9:29])),
      new = as.numeric(gameofclones:::sad_leslie(L)))


# Get max abundance given leslie matrix and K:
.K = 12500
.K * (max(abs(eigen(L)[["values"]])) - 1)
.K * (max(abs(eigen(L0)[["values"]])) - 1)


.P <- numeric(2)
.P[1] <- (gameofclones::wasp_attack[["rel_attack"]][1:4] * sapply(1:4, \(i) sum(X0[(2*i-1):(2*i)])) / sum(X0[1:8])) |> sum()
.P[2] <- gameofclones::wasp_attack[["rel_attack"]][5]



run_paras <- function(L, aphids0 = 10, paras0 = 1,
                      a = 2.32, k = 0.35, h = 0.008, P = .P,
                      K = 12500, s_y = 0.69, max_t = 1000,
                      .plot = TRUE) {

    aphids <- numeric(max_t+1)
    paras <- numeric(max_t+1)
    aphids[1] <- aphids0
    paras[1] <- paras0
    X <- gameofclones:::sad_leslie(L) * aphids0
    P <- unname(cbind(P))
    for (t in 1:max_t) {
        A <- (1 + P * a * paras[t] / (k * (h * aphids[t] + 1)))^(-k)
        S <- 1 / (1 + aphids[t] / K);
        paras[t+1] = paras[t] * s_y + sum(S * (1 - A) * (L %*% X))
        X <- S * A * (L %*% X)
        aphids[t+1] <- sum(X)
    }
    d <- tibble(time = 0:max_t, aphids = aphids, parasitoid = paras) |>
        pivot_longer(-time, names_to = "species", values_to = "density")
    if (!.plot) return(d)
    p <- d |>
        ggplot(aes(time, density, color = species)) +
        geom_hline(yintercept = 0, color = "gray70") +
        geom_line(linewidth = 0.75) +
        scale_color_manual(values = c("dodgerblue", "firebrick")) +
        theme_classic()
    return(p)
}




run_paras(L, a = 2.32, k = 0.35, h = 0.008, P = .P)

run_paras(L, a = 1e-3, h = 0.002, k = 1, P = c(1,1))
