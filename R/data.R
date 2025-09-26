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
#'   \item{`alate_0`}{Alate offspring proportion when
#'                    aphid density is low. The proportion of winged offspring
#'                    from apterous aphids is `inv_logit(alate_0 + alate_1 * z)`
#'                    where `z` is the total number of aphids on that plant.}
#'   \item{`alate_1`}{How strongly aphid density affects alate production.
#'                    See `alate_0` above for the equation.}
#'   \item{K}{Unparasitized aphid density dependence.}
#'   \item{m}{Aphid and mummy mortality due to generalist predators.}
#'   \item{a}{Parasitoid attack rate.}
#'   \item{h}{Parasitoid handling time.}
#'   \item{k}{Parasitoid aggregation parameter.}
#'   \item{s_y}{Adult parasitoid daily survival.}
#'   \item{s_p}{Parasitized aphid daily survival.}
#'   \item{K_p}{Parasitized aphid density dependence.}
#'   \item{R}{Relative parasitoid attack rates by aphid stage.}
#'   \item{theta_m}{Proportion of mummies that transition to adult parasitoids each day.}
#'   \item{theta_p}{Proportion of parasitized aphids that transition to mummies each day.}
#' }
#'
#' All are single numeric values, except for `R` that is a vector of length 4.
#'
#' @source <https://www.github.com/lucasnell/gameofclones>
"pop_info"



