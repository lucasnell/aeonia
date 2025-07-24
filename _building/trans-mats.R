library(tidyverse)
library(gameofclones)
library(aphidsync)  # calc_lambda
library(aeonia)



line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,10), rep(0, 5)),
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


.p <- numeric(2)
.p[1] <- (gameofclones::wasp_attack[["rel_attack"]][1:4] * sapply(1:4, \(i) sum(X0[(2*i-1):(2*i)])) / sum(X0[1:8])) |> sum()
.p[2] <- gameofclones::wasp_attack[["rel_attack"]][5]



L <- rbind(c(pop_info$surv_j, pop_info$fecund),
           c(pop_info$recruit, pop_info$surv_a))

run_paras <- function(L, aphids0 = 10, paras0 = 1,
                      a = 2.32, k = 0.35, h = 0.008, p = c(1,1),
                      K = 12500, s_y = 0.69, s_p =  0.9745, s_m = 0.9, max_t = 1000,
                      K_p_mult = 0.6369427,
                      .plot = TRUE) {

    K_p <- K * K_p_mult
    if (length(p) == 1) p <- rep(p, 2)

    aphids <- numeric(max_t+1)
    paras <- numeric(max_t+1)
    aphids[1] <- aphids0
    paras[1] <- paras0
    X <- gameofclones:::sad_leslie(L) * aphids0 # unparasitized aphids
    P <- 0  # parasitized, alive aphids
    M <- 0  # parasitized, dead aphids
    p <- unname(cbind(p))
    for (t in 1:max_t) {
        A <- (1 + p * a * paras[t] / (k * (h * sum(X) + 1)))^(-k)
        S <- 1 / (1 + aphids[t] / K)
        # >>>>>>>>>>>>>>>>
        # # not including P and M:
        paras[t+1] = paras[t] * s_y + sum((1 - A) * (L %*% X))
        X <- S * A * (L %*% X)
        aphids[t+1] <- sum(X)
        # <<<<<<<<<<<<<<<<
        # S_p <- 1 / (1 + aphids[t] / K_p)
        # new_adult_paras <- 0.3125 * M
        # new_M <- 0.0928 * P
        # new_P <- sum((1 - A) * (L %*% X))
        # paras[t+1] = paras[t] * s_y + 0.5 * new_adult_paras
        # M <- s_m * (M - new_adult_paras + new_M)
        # if (M < 0) M <- 0
        # P <- S_p * (s_p * P - new_M + new_P)
        # if (P < 0) P <- 0
        # X <- S * A * (L %*% X)
        # aphids[t+1] <- sum(X) + P
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


run_paras2 <- function(aphids0 = 10, paras0 = 1,
                       a = wasp_attack$a, k = wasp_attack$k, h = wasp_attack$h,
                       s = 0.69,
                       K = 12500, max_t = 1000) {
    test_insect_pops(max_t = max_t, A0 = aphids0, W0 = 0, P0 = paras0, B = 0,
                     a = a, h = h, k = k, s = s, alate_0 = -Inf, alate_1 = 0) |>
        select(-alates) |>
        rename(parasitoid = enemies) |>
        pivot_longer(-time, names_to = "species", values_to = "density")
}



# .h = 3;
# .h = wasp_attack$h * 1; .k = wasp_attack$k * 1
.a = wasp_attack$a * 1; .h = wasp_attack$h * 625; .k = wasp_attack$k * 1
# run_paras(L, a = .a, h = .h, k = .k, max_t = 1e3, p = 1, .plot = FALSE) |>
run_paras2(a = .a, h = .h, k = .k, max_t = 1e3) |>
    filter(time > 120) |>
    mutate(time = time - min(time)) |>
    ggplot(aes(time, log10(0.5 + density), color = species)) +
    geom_hline(yintercept = log10(0.5), color = "gray70") +
    geom_line(linewidth = 0.75) +
    scale_color_manual(values = c("dodgerblue", "firebrick")) +
    scale_x_continuous(breaks = 0:5 * 200) +
    theme_classic() +
    ggtitle(sprintf("a = %.2f, k = %.2f, h = %.2g", .a, .k, .h))


sims <- sim_experiments(line_s, 1, 1000, alate_b0 = -Inf, alate_b1 = 0,
                        wasp_density_0 = 1, wasp_delay = 0, extinct_N = 0)


left_join(sims$aphids |>
              filter(!is.na(type)) |>
              group_by(time) |>
              summarize(aphids = sum(N)),
          sims$wasps |>
              group_by(time) |>
              summarize(wasps = sum(wasps)),
          by = join_by(time)) |>
    pivot_longer(aphids:wasps, names_to = "species", values_to = "density") |>
    ggplot(aes(time, density, color = species)) +
    geom_line() +
    scale_color_viridis_d(begin = 0.2) +
    facet_wrap( ~ species, scales = "free_y", ncol = 1) +
    scale_x_continuous(breaks = 0:5 * 200) +
    theme_classic()

