## code to prepare `pop_info` dataset goes here

library(gameofclones)


# ------------------------*
# Background
# ------------------------*

# Susceptible aphid line from `gameofclones`:
line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,10), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high")
# Leslie matrix (for apterous, non-parasitized):
L0 <- clonal_line("susceptible",
                  density_0 = cbind(c(0,0,0,0,10), rep(0, 5)),
                  surv_juv_apterous = "high",
                  surv_adult_apterous = "high",
                  repro_apterous = "high",
                  p_instar_smooth = 0) |>
    getElement("leslie") |>
    base::`[`(,,1)
# And calculate stable age distribution for it:
X0 <- gameofclones:::sad_leslie(L0)

# Used below to calculate growth rate from Leslie matrix:
calc_lambda <- function(L) {
    eigen(L, only.values = TRUE)[["values"]] |>
        abs() |>
        max()
}


# ------------------------*
# Calculate parameters
# ------------------------*

# Juveniles staying juveniles:
# Note: indices in `X0[1:7]` and `X0[1:8]` below are intentional!
sj <- sum(L0[row(L0) - col(L0) == 1][1:7] * (X0[1:7] / sum(X0[1:8])))

# Recruitment: juveniles transitioning to adults:
r <- L0[row(L0) - col(L0) == 1][8] * (X0[8] / sum(X0[1:8]))

# Adults staying adults:
sa <- sum(L0[row(L0) - col(L0) == 1][9:28] * (X0[9:28] / sum(X0[9:28])))

# Fecundity:
f <- sum(L0[1,9:29] * (X0[9:29] / sum(X0[9:29])))


# ------------------------*
# Create simplified Leslie and compare to old
# ------------------------*

# Create simplified Leslie:
L <- rbind(c(sj, f),
           c(r, sa))

# Are growth rates quite close?
cbind(orig = calc_lambda(L0), new = calc_lambda(L))
#          orig      new
# [1,] 1.331508 1.331524

# Are stable age distributions quite close?
cbind(orig = c(sum(X0[1:8]), sum(X0[9:29])),
      new = as.numeric(gameofclones:::sad_leslie(L)))
#            orig       new
# [1,] 0.91119723 0.9111939
# [2,] 0.08880277 0.0888061



# ------------------------*
# Define other population info
# ------------------------*

# These are just the defaults from this package that were vetted to work
# well in the experiments:
K <- eval(formals(gameofclones::sim_experiments)[["K"]])
K_p_mult <- eval(formals(gameofclones::sim_experiments)[["K_y_mult"]])
s_p <- line_s$leslie[,,3][2,1]
pred_surv <- 1 - eval(formals(gameofclones::sim_experiments)[["pred_rate"]])
a0 <- eval(formals(gameofclones::sim_experiments)[["alate_b0"]])
a1 <- eval(formals(gameofclones::sim_experiments)[["alate_b1"]])
a <- wasp_attack$a
h <- wasp_attack$h
k <- wasp_attack$k
s_y <- populations$s_y
# Relative attack rates by stage:
R <- (X0 * rep(wasp_attack$rel_attack, c(rep(2, 4), 21))) |>
    # Group by juvenile vs adult:
    (\(x) c(sum(x[1:8]), sum(x[-1:-8])))() |>
    # repeat for winged forms, where adult winged are not attacked at all:
    (\(x) c(x, x[1], 0))()




# ------------------------*
# Optimizations for trans_pm and trans_ma
# ------------------------*

# Simplified R version of host-parasitoid model (no alates):
run_paras_pm1 <- function(aphids0, wasps0, max_t = 100,
                          trans_ma = 1,
                          trans_pm = 1) {
    R <- R[1:2]  # <:-- Because we're not simulating alates
    K_p <- K * K_p_mult
    aphids <- numeric(max_t+1)
    wasps <- numeric(max_t+1)
    aphids[1] <- aphids0
    wasps[1] <- wasps0
    X <- gameofclones:::sad_leslie(L) * aphids0 # unparasitized aphids
    P <- 0  # parasitized, alive aphids
    M <- 0  # mummies
    for (t in 1:max_t) {
        A <- (1 + R * a * wasps[t] / (k * (h * sum(X) + 1)))^(-k)
        S <- 1 / (1 + aphids[t] / K)
        S_p <- 1 / (1 + aphids[t] / K_p)
        new_adult_wasps <- trans_ma * M
        new_M <- trans_pm * P
        new_P <- sum((1 - A) * (L %*% X))
        wasps[t+1] = wasps[t] * s_y + 0.5 * new_adult_wasps
        M <- pred_surv * (M + new_M) - new_adult_wasps
        if (M < 0) stop("M < 0")
        P <- pred_surv * S_p * (s_p * P + new_P) - new_M
        if (P < 0) stop("P < 0")
        X <- pred_surv * S * A * (L %*% X)
        aphids[t+1] <- sum(X) + P
    }
    d <- tibble(time = 0:max_t, aphids = aphids, wasps = wasps) |>
        pivot_longer(-time, names_to = "species", values_to = "density")
    return(d)
}
# Run full model from `gameofclones`:
run_paras_got <- function(aphids0, wasps0, max_t = 100) {

    .line_s <- line_s
    .line_s$density_0 <- cbind(X0, 0) * aphids0
    sim_experiments(.line_s, 1, max_t = max_t, alate_b0 = -Inf, alate_b1 = 0,
                    wasp_density_0 = wasps0, wasp_delay = 0, extinct_N = 0) |>
        (\(sims) {
            left_join(sims$aphids |>
                          filter(type != "mummy", type != "parasitized") |>
                          group_by(time) |>
                          summarize(aphids = sum(N)),
                      sims$wasps |>
                          group_by(time) |>
                          summarize(wasps = sum(wasps)),
                      by = join_by(time))
        })() |>
        pivot_longer(aphids:wasps, names_to = "species", values_to = "density")
}
combos <- crossing(a0 = 10^(0:2), w0 = 10^(-3:-1)) |>
    mutate(w0 = a0 * w0)
sims_full <- pmap(combos, \(a0, w0) {
    run_paras_got(aphids0 = a0, wasps0 = w0) |>
        mutate(combo = paste(a0, w0, sep = "_"))
}) |>
    list_rbind() |>
    mutate(combo = factor(combo),
           sim_type = "full")
sse <- function(sims1, sims2) {
    stopifnot(nrow(sims1) == nrow(sims2))
    errs <- sims1[["density"]][sims1$species == "aphids"] -
        sims2[["density"]][sims2$species == "aphids"]
    return(sum(errs^2))
}
fit_pm1 <- function(pars) {
    pars <- inv_logit(pars)
    sims <- pmap(combos, \(a0, w0) {
        run_paras_pm1(aphids0 = a0, wasps0 = w0,
                      trans_ma = pars[1],
                      trans_pm = pars[2])
    }) |>
        list_rbind()
    .sse <- sse(sims, sims_full)
    return(.sse)
}

set.seed(1974834691)
op_pm1 <- optim(logit(c(0.25, 0.1)), fit_pm1)
# inv_logit(op_pm1$par)
# [1] 0.05206727 0.17706238
trans_ma <- inv_logit(op_pm1$par[1])
trans_pm <- inv_logit(op_pm1$par[2])





# ------------------------*
# Create dataset
# ------------------------*
pop_info <- list(surv_j = sj, surv_a = sa, recruit = r, fecund = f,
                 alate_0 = a0, alate_1 = a1, K = K, pred_surv = pred_surv,
                 a = a, h = h, k = k, s_y = s_y,
                 s_p = s_p, K_p_mult = K_p_mult, R = R,
                 trans_ma = trans_ma, trans_pm = trans_pm)

usethis::use_data(pop_info, overwrite = TRUE)
