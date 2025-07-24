## code to prepare `pop_info` dataset goes here

library(gameofclones)


# ------------------------*
# Background
# ------------------------*

# Creates Leslie matrix for previously studied pea aphid clonal line:
line_s <- clonal_line("susceptible",
                      density_0 = cbind(c(0,0,0,0,32), rep(0, 5)),
                      surv_juv_apterous = "high",
                      surv_adult_apterous = "high",
                      repro_apterous = "high",
                      p_instar_smooth = 0)
# Extract that Leslie matrix to base simpler one off of:
L0 <- line_s$leslie[,,1]
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

# Are stable age distributions quite close?
cbind(orig = c(sum(X0[1:8]), sum(X0[9:29])),
      new = as.numeric(gameofclones:::sad_leslie(L)))


# ------------------------*
# Define other population info
# ------------------------*

# These are just the defaults from this package that were vetted to work
# well in the experiments:
.K <- eval(formals(gameofclones::sim_experiments)[["K"]])
a0 <- eval(formals(gameofclones::sim_experiments)[["alate_b0"]])
a1 <- eval(formals(gameofclones::sim_experiments)[["alate_b1"]])
a <- wasp_attack$a
h <- wasp_attack$h
k <- wasp_attack$k
s <- populations$s_y


# ------------------------*
# Create dataset
# ------------------------*
pop_info <- list(surv_j = sj, surv_a = sa, recruit = r, fecund = f,
                 alate_0 = a0, alate_1 = a1, K = .K,
                 a = a, h = h, k = k, s = s)

usethis::use_data(pop_info, overwrite = TRUE)
