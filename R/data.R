


#' Development times for aphids and wasps.
#'
#' See `vignette("internals", "aeonia")` for the exact values here.
#'
#'
#' @format A list with two items:
#' \describe{
#'   \item{instar_days}{Number of days per instar, for low (\code{lowT}) and
#'         high (\code{highT}) temperatures (20º and 27º C).}
#'   \item{mum_days}{1x2 matrix with the number of days per stage of an aphid
#'                   being parasitized: living and dead ("mummy"), respectively.}
#' }
#'
#' @source \url{http://doi.wiley.com/10.1890/13-1933.1}
#'
"an_dev_times"


#' Population rates and starting values for aphids and wasps.
#'
#' See `vignette("internals", "aeonia")` for more information and
#' for the exact values used.
#'
#'
#' @format A list with ten items:
#' \describe{
#'   \item{surv_juv}{
#'       List of length two, each item a length-of-one numeric vector containing
#'       juvenile daily survival for aphid lines with low and high survivals
#'       (\code{low} and \code{high}).
#'   }
#'   \item{surv_adult}{
#'       List of length two, each item a 1x200 matrix containing
#'       adults daily survival for aphid lines with low and high survivals
#'       (\code{low} and \code{high}).
#'       Although each matrix has 200 columns, most of them are filled with zeros.
#'   }
#'   \item{repro}{
#'       List of length two, each item a 1x200 matrix containing
#'       daily reproduction for aphid lines with low and high reproduction
#'       (\code{low} and \code{high}).
#'       Although each matrix has 200 columns, most of them are filled with zeros.
#'   }
#'   \item{K}{
#'     A single number representing aphid density dependence.
#'   }
#'   \item{K_p_mult}{
#'     A single number representing the number multiplied by `K` to get
#'     parasitized aphid density dependence.
#'   }
#'   \item{s_y}{
#'     A single number representing parasitoid adult daily survival.
#'   }
#'   \item{sex_ratio}{
#'     A single number representing proportion of female wasps.
#'   }
#' }
#'
#' @source \url{http://doi.wiley.com/10.1890/13-1933.1}
#'
"an_populations"


#' Wasp attack rate parameters.
#'
#' See `vignette("internals", "aeonia")` for more information and
#' for the exact values used.
#'
#'
#' @format A list with five items:
#' \describe{
#'   \item{a}{
#'     parasitoid attack rate
#'   }
#'   \item{k}{
#'     aggregation parameter of the negative binomial distribution
#'   }
#'   \item{h}{
#'     parasitoid attack rate handling time
#'   }
#'   \item{rel_attack}{
#'     A vector of length 5, giving  relative attack rates on the
#'     different instars.
#'   }
#'   \item{attack_surv}{
#'     A numeric vector of length two, of survivals of singly attacked and
#'     multiply attacked resistant aphids.
#'     \emph{Note:} This is not from either paper, but from unpublished code by
#'     Anthony Ives.
#'   }
#' }
#'
#' @source \url{http://doi.wiley.com/10.1890/13-1933.1}
#' @source \url{http://www.journals.uchicago.edu/doi/10.1086/303269}
#'
"an_wasp_attack"


#' Parameters associated with environmental effects and stochasticity.
#'
#' See `vignette("internals", "aeonia")` for more information and
#' for the exact values used.
#'
#'
#' @format A list of length 9:
#' \describe{
#'   \item{disp_start}{
#'     List of length 2, each item of which contains a 1-length numeric vector
#'     indicating the day at which aphids begin dispersing for 20ºC  (\code{lowT})
#'     and 27ºC (\code{highT}).
#'     It's assumed that only adults are dispersing, but the day at which
#'     this occurs depends on how quickly the aphids are developing.
#'   }
#'   \item{sigma_x}{
#'     Numeric vector of length 1, indicating environmental std dev for aphids.
#'   }
#'   \item{sigma_y}{
#'     Numeric vector of length 1, indicating environmental std dev for wasps.
#'   }
#'   \item{rho}{
#'     Numeric vector of length 1, indicating environmental correlation among instars.
#'   }
#'   \item{radius}{
#'     Numeric vector of length 1, indicating radius of alate dispersal.
#'     Value taken from \url{https://doi.org/10.3929/ethz-a-000909462}.
#'   }
#' }
#'
#' @source \url{http://doi.wiley.com/10.1890/13-1933.1}
#'
"an_environ"

