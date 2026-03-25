# library(bench)
#
#
# Rcpp::cppFunction(code =
# 'SEXP make_ptr(const arma::uword& r, const arma::uword& c) {
#     XPtr<arma::Mat<arma::uhword>> xp(new arma::Mat<arma::uhword>(
#                                      r, c, arma::fill::none), true);
#     arma::Mat<arma::uhword>& x(*xp);
#     typedef arma::uword uint32;
#     for (uint32 i = 0; i < x.n_rows; i++) {
#         for (uint32 j = 0; j < x.n_cols; j++) {
#             arma::uhword& xij(x.at(i,j));
#             xij = (arma::uhword)0U;
#             if (R::runif(0,1) > 0.5) xij += (arma::uhword)1;
#             if (R::runif(0,1) > 0.5) xij += (arma::uhword)2;
#         }
#     }
#     // double x_max = arma::max(arma::max(x));
#     // double x_min = arma::min(arma::min(x));
#     // Rcout << "min = " << x_min << ", ";
#     // Rcout << "max = " << x_max << std::endl;
#     return xp;
# }
# ', depends = "RcppArmadillo")
#
#
#
#
# Rcpp::cppFunction(code =
# 'arma::Mat<arma::uhword> getset_bits(SEXP x_ptr) {
#     typedef arma::uhword uint16;
#     typedef arma::uword uint32;
#     XPtr<arma::Mat<uint16>> xp(x_ptr);
#     const arma::Mat<uint16>& x(*xp);
#
#     arma::Mat<uint16> z(arma::size(x), arma::fill::none);
#     uint16 k0 = 0;
#     uint16 k1 = 1;
#     uint16 mask1 = 1;
#     bool b;
#     for (uint32 i = 0; i < x.n_rows; i++) {
#         for (uint32 j = 0; j < x.n_cols; j++) {
#             const uint16& xij(x.at(i,j));
#             uint16& zij(z.at(i,j));
#             zij = 0;
#             // b = (xij & ( 1 << k0 )) >> k0;
#             b = xij & mask1;
#             if (b) zij |= (mask1 << k0);
#             // b = (xij & ( 1 << k1 )) >> k1;
#             b = (xij >> mask1) & mask1;
#             if (b) zij |= ((uint32)1 << k1);
#         }
#     }
#     return z;
# }
# ', depends = "RcppArmadillo")
# Rcpp::cppFunction(code =
# 'arma::Mat<arma::uhword> getset_bool(SEXP x_ptr) {
#     typedef arma::uhword uint16;
#     typedef arma::uword uint32;
#     XPtr<arma::Mat<uint16>> xp(x_ptr);
#     const arma::Mat<uint16>& x(*xp);
#
#     arma::Mat<uint16> z(arma::size(x), arma::fill::none);
#     bool b;
#     for (uint32 i = 0; i < x.n_rows; i++) {
#         for (uint32 j = 0; j < x.n_cols; j++) {
#             const uint16& xij(x.at(i,j));
#             uint16& zij(z.at(i,j));
#             zij = 0;
#             b = xij == 1U || xij == 3U;
#             if (b) zij += 1U;
#             // b = xij == 2U || xij == 3U;
#             b = xij > 1U;
#             if (b) zij += 2U;
#         }
#     }
#     return z;
# }
# ', depends = "RcppArmadillo")
# # Rcpp::cppFunction(code =
# # 'arma::umat getset_2bools(const LogicalMatrix& x, const LogicalMatrix& y) {
# #     arma::umat z(x.nrow(), x.ncol(), arma::fill::none);
# #     typedef arma::uword uint32;
# #     for (size_t i = 0; i < x.nrow(); i++) {
# #         for (size_t j = 0; j < x.ncol(); j++) {
# #             uint32& zij(z.at(i,j));
# #             if (x[i,j]) zij += 1U;
# #             if (y[i,j]) zij += 2U;
# #         }
# #     }
# #     return z;
# # }
# # ', depends = "RcppArmadillo")
#
#
# n <- 100L
# ptr <- make_ptr(100, 100)
#
#
#
# mark(bits = getset_bits(ptr),
#      bool = getset_bool(ptr),
#      iterations = 2000, check = FALSE, memory = FALSE)





#'
#' Find par values for larger landscape simulations
#' 04-large-plantscapes-par-vals.sh must be run on the cluster first, then
#' its output sent to ./_scripts/interm-data/
#'


source("_scripts/03-large-preamble.R")




test_sims <- read_rds("_scripts/interm-data/large-plantscapes-par-vals.rds") |>
    # # determine which two scenario: Pseudomonas promotes or inhibits viruses, resp.:
    # mutate(scenario = factor(zeta < 0.2 & N0 < 20, levels = c(TRUE, FALSE),
    #                          labels = c("promotes", "inhibits"))) |>
    mutate(scenario = interaction(Y0, N0, zeta, drop = TRUE)) |>
    select(scenario, p_load, everything()) |>
    mutate(outbreak_size = outbreak_size) |>
    pivot_longer(outbreak_size:p_emerge, names_to = "outcome") |>
    filter(!is.na(value)) |>
    mutate(outcome = factor(outcome,
                            levels = c("outbreak_size", "p_emerge"),
                            labels = c("outbreak size", "prob. emerge")),
           n_pseudo = factor(n_pseudo)) |>
    arrange(scenario, p_load, outcome, n_pseudo)


np_pal <- c("#999999", "#0046D2") |>
    set_names(levels(test_sims$n_pseudo))



# Based on below, I'll choose p_load = 0.075
test_sims |>
    # group_by(p_load, scenario) |>
    # mutate(within = max(value[outcome == "prob. emerge"]) < 1 &
    #            min(value[outcome == "prob. emerge"]) > 0) |>
    # ungroup() |>
    # filter(within) |>
    # filter(zeta == min(zeta)) |>
    # filter(zeta == 0.1, p_load == 0.15) |>
    # filter(p_load == 0.1) |>
    filter(outcome == "outbreak size") |>
    group_by(p_load, scenario) |>
    summarize(diff = value[n_pseudo != "0"] / value[n_pseudo == "0"], .groups = "drop") |>
    # filter(diff > 0) |>
    arrange((abs(diff)))


scenarios = c(low = "200.60.1", high = "200.45.0.125")
p_load = 0.1




test_sims |>
    # filter(p_load == .env$p_load) |>
    # filter(p_load == min(p_load)) |>
    filter(p_load == max(p_load)) |>
    filter(scenario %in% scenarios) |>
    arrange(scenario, outcome, n_pseudo)


#     ggplot(aes(Y0, value, color = factor(N0), shape = n_pseudo)) +
#     geom_hline(yintercept = 0) +
#     geom_point(position = position_jitterdodge(jitter.width = 0.0, dodge.width = 3)) +
#     # geom_point() +
#     facet_wrap(~ outcome, ncol = 1, scales = "free") +
#     # scale_color_manual(values = np_pal)
#     scale_color_viridis_d(end = 0.9) +
#     scale_shape_manual(values = c(`0` = 1, `7000` = 8))


test_sims |>
    # filter(scenario == "100.60.1", p_load == 0.075) |>
    filter(scenario == "105.60.1", p_load == 0.2) |>
    ggplot(aes(p_load, value, color = factor(N0), shape = n_pseudo)) +
    geom_hline(yintercept = 0) +
    geom_point(position = position_jitterdodge(jitter.width = 0.0, dodge.width = 3)) +
    # geom_point() +
    facet_wrap(~ outcome, ncol = 1, scales = "free") +
    # scale_color_manual(values = np_pal)
    scale_color_viridis_d(end = 0.9) +
    scale_shape_manual(values = c(`0` = 1, `7000` = 8))



large_simmer <- function(landscape, Y0, N0, zeta, p_load) {

    args <- list(landscape = landscape,
                 sd_N = 0,
                 virus_attract = 1,
                 pseudo_repel = 1,
                 Y0 = Y0,
                 N0 = N0,
                 zeta = zeta,
                 p_load_alate = p_load,
                 p_load_plant = p_load,
                 K = 12.5e3,
                 pseudo_surv = 0.85,
                 n_sims = 100L,
                 summ = "all")

    return(do.call(big_plantscape, args))

}



landscape1 <- sim_df |>
    filter(n_pseudo == 7000, wt_vp == 1e-6, wt_pp == 1,
           # These do not affect landscape:
           type == "low", sd_N == 0, virus_attract == 1,
           pseudo_repel == 1) |>
    getElement("landscape") |> getElement(1)
landscape0 <- array(c(1L, rep(0L, 99999L)), c(100L, 100L, 1L))



one_test <- function(p_load, Y0, N0, zeta) {

    n_inf1 <- large_simmer(landscape = landscape1,
                           Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
        getElement("n_infected")
    n_inf0 <- large_simmer(landscape = landscape0,
                           Y0 = Y0, N0 = N0, zeta = zeta, p_load = p_load) |>
        getElement("n_infected")

    tibble(n_pseudo = c(7000L, 0L),
           p_load = .env$p_load, Y0 = .env$Y0, N0 = .env$N0, zeta = .env$zeta,
           outbreak_size = c(mean(n_inf1[n_inf1 > 1]), mean(n_inf0[n_inf0 > 1])),
           p_emerge = c(mean(n_inf1 > 1), mean(n_inf0 > 1)))

}

# Both are with virus_attract = 5
#
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000   0.05   300    35   0.1         23.7      0.99
# 2        0   0.05   300    35   0.1          5.40     0.84
#
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000   0.05   220    60     1          3.64     0.74
# 2        0   0.05   220    60     1         23.6      0.99


# > (sim <- one_test(p_load = 0.5, Y0 = 220, N0 = 60, zeta = 1))
# # A tibble: 2 × 7
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000    0.5   220    60     1         1254.        1
# 2        0    0.5   220    60     1         5296.        1
# > (sim <- one_test(p_load = 0.5, Y0 = 220, N0 = 90, zeta = 1))
# # A tibble: 2 × 7
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000    0.5   220    90     1         2105.        1
# 2        0    0.5   220    90     1         6690.        1
# > (sim <- one_test(p_load = 0.5, Y0 = 220, N0 = 100, zeta = 1))
# # A tibble: 2 × 7
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000    0.5   220   100     1         2202.        1
# 2        0    0.5   220   100     1         7080.        1
# > (sim <- one_test(p_load = 0.5, Y0 = 220, N0 = 150, zeta = 1))
# # A tibble: 2 × 7
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000    0.5   220   150     1         3476.        1
# 2        0    0.5   220   150     1         8888.        1
# > (sim <- one_test(p_load = 0.5, Y0 = 280, N0 = 150, zeta = 1))
# # A tibble: 2 × 7
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000    0.5   280   150     1         2492.        1
# 2        0    0.5   280   150     1         7917.        1
# > (sim <- one_test(p_load = 0.5, Y0 = 200, N0 = 150, zeta = 1))
# # A tibble: 2 × 7
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000    0.5   200   150     1         3652.        1
# 2        0    0.5   200   150     1         8986.        1
# > (sim <- one_test(p_load = 0.5, Y0 = 150, N0 = 150, zeta = 1))
# # A tibble: 2 × 7
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000    0.5   150   150     1         4786.        1
# 2        0    0.5   150   150     1         9518.        1

# > (sim <- one_test(p_load = 0.5, Y0 = 150, N0 = 20, zeta = 0.1))
# # A tibble: 2 × 7
#   n_pseudo p_load    Y0    N0  zeta outbreak_size p_emerge
#      <int>  <dbl> <dbl> <dbl> <dbl>         <dbl>    <dbl>
# 1     7000    0.5   150    20   0.1         9415.        1
# 2        0    0.5   150    20   0.1         3763.        1




# srun -N 1 -n 1 -c 50 --mem=100G --time=1-20:00:00 --job-name="aeonia-test" --pty bash -l
# Takes ~ 10 min (1 min on cluster w 50 threads)
(sim <- one_test(p_load = 0.2, Y0 = 210, N0 = 60, zeta = 1))
# 1 / 0.127 = 7.874016

for (n0 in c(5, 20, 45, 90, 180)) for (y0 in c(5, 25, 50, 75, 100, 200, 300, 400)) {
    cat("Y0 =", y0, ", N0 =", n0, "\n")
    print(one_test(p_load = 0.2, Y0 = y0, N0 = n0, zeta = 0.1))
    cat("\n\n")
}






p_load <- 0.075

#'
#' From below graph and table, I should use
#' high zeta (0.12) when *Pseudomonas* promotes outbreaks, and
#' medium zeta (0.95) when *Pseudomonas* inhibits outbreaks.
#' Because I simulated 5000 *Pseudomonas* patches instead of 7000 (the density that
#' typically produces the strongest *Pseudomonas* effect),
#' I will use zeta = 0.14 when *Pseudomonas* promotes outbreaks and
#' zeta = 1.0 when *Pseudomonas* inhibits outbreaks.
#'


test_sims |>
    filter(p_load == .env$p_load) |>
    select(-p_load) |>
    filter(Y0 == 100) |>
    filter(N0 %in% c(10, 110)) |>
    group_by(scenario) |>
    mutate(zeta_lvl = case_when(zeta == min(zeta) ~ 1L,
                                zeta == max(zeta) ~ 3L,
                                .default = 2L)) |>
    ungroup() |>
    ggplot(aes(zeta_lvl, value, color = n_pseudo)) +
    # geom_hline(yintercept = 0) +
    geom_point(aes(shape = scenario)) +
    geom_line(aes(linetype = scenario)) +
    scale_x_continuous(breaks = 1:3, labels = c("low", "mid", "high")) +
    facet_grid(outcome ~ ., scales = "free") +
    scale_color_manual(values = np_pal)


test_sims |>
    filter(p_load == .env$p_load) |>
    filter(Y0 == 100) |>
    filter(N0 %in% c(10, 110)) |>
    select(-p_load, -Y0, -N0) |>
    group_by(scenario) |>
    mutate(zeta_lvl = case_when(zeta == min(zeta) ~ 1L,
                            zeta == max(zeta) ~ 3L,
                            .default = 2L)) |>
    ungroup() |>
    filter((scenario == "promotes" & zeta_lvl == 3L) |
               (scenario == "inhibits" & zeta_lvl == 2L)) |>
    arrange(outcome, scenario, n_pseudo) |>
    select(outcome, scenario, n_pseudo, value)



