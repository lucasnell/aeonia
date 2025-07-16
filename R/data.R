#' Basic pea aphid population information
#'
#'
#' @format ## `pop_info`
#' A named list with the following items:
#' \describe{
#'   \item{surv_j}{Single numeric indicating aphid juvenile survival.}
#'   \item{surv_a}{Single numeric indicating aphid adult survival.}
#'   \item{recruit}{Single numeric indicating aphid recruitment.}
#'   \item{fecund}{Single numeric indicating aphid fecundity.}
#'   \item{K}{Single numeric indicating pea aphid density dependence.}
#'   \item{alate_0}{Single numeric affecting alate offspring proportion when
#'                  aphid density is low. The proportion of winged offspring
#'                  from apterous aphids is `inv_logit(alate_0 + alate_1 * z)`
#'                  where `z` is the total number of aphids on that plant.}
#'   \item{alate_1}{Single numeric affecting how strongly aphid density affects
#'                  alate production. See `alate_0` above for the equation.}
#' }
#' @source <https://www.github.com/lucasnell/gameofclones>
"pop_info"



