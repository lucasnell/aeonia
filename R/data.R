#' Basic pea aphid population information
#'
#'
#' @format ## `pop_info`
#' A named list with the following items:
#' \describe{
#'   \item{`surv_j`}{Aphid juvenile survival.}
#'   \item{`surv_a`}{Aphid adult survival.}
#'   \item{`recruit`}{Aphid recruitment.}
#'   \item{`fecund`}{Aphid fecundity.}
#'   \item{`alate_infl`}{Inflection point for sigmoid relationship between aphid
#'                       density and alate offspring proportion.
#'                       The proportion of winged offspring
#'                       from apterous aphids is `1 / {1 + 10^((alate_infl - z) * alate_slope)}`
#'                       where `z` is the total number of aphids on that plant.}
#'   \item{`alate_slope`}{Slope for sigmoid relationship between aphid
#'                        density and alate offspring proportion.
#'                        See `alate_infl` above for the equation.}
#'   \item{`K`}{Unparasitized aphid density dependence.}
#'   \item{`pred_surv`}{Aphid and mummy survival from generalist predators.}
#'   \item{`a`}{Parasitoid attack rate.}
#'   \item{`h`}{Parasitoid handling time.}
#'   \item{`k`}{Parasitoid aggregation parameter.}
#'   \item{`s_y`}{Adult parasitoid daily survival.}
#'   \item{`s_p`}{Parasitized aphid daily survival.}
#'   \item{`K_p_mult`}{Multiplier to get parasitized aphid density dependence (`K_p = K_p_mult * K`).}
#'   \item{`R`}{Relative parasitoid attack rates by aphid stage.}
#'   \item{`trans_ma`}{Proportion of mummies that transition to adult parasitoids each day.}
#'   \item{`trans_pm`}{Proportion of parasitized aphids that transition to mummies each day.}
#' }
#'
#' All are single numeric values, except for `R` that is a vector of length 4.
#'
#' @source <https://www.github.com/lucasnell/gameofclones>
"pop_info"



