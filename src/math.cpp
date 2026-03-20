
#include "math.hpp"


using namespace Rcpp;

/*
 =====================================================================================
 =====================================================================================
 Logit and inverse logit
 =====================================================================================
 =====================================================================================
 */


//' Logit and inverse logit functions.
//'
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector logit(NumericVector p) {
    NumericVector out(p.size());
    for (uint32 i = 0; i < p.size(); i++) {
        logit__(p[i], out[i]);
    }
    return out;
}
//' @describeIn logit
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector inv_logit(NumericVector a){
    NumericVector out(a.size());
    for (uint32 i = 0; i < a.size(); i++) {
        inv_logit__(a[i], out[i]);
    }
    return out;
}




/*
 =====================================================================================
 =====================================================================================
 Stable age distribution
 =====================================================================================
 =====================================================================================
 */

//' Compute stable age distribution from transition matrix.
//'
//' The stable age distribution is the proportion of different stages that
//' is required for the population to grow exponentially.
//'
//' @param L Numeric transition matrix. Must be square.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector sad_leslie(arma::mat L) {
    if (!L.is_square()) stop("L must be square");
    arma::vec out;
    sad_leslie__(L, out);
    return wrap(out);
}

/*
 Create Leslie matrix from aphid info

 @param instar_days Integer vector of the number of stages (days) per aphid instar.
 @param surv_juv Single numeric of daily juvenile survival.
 @param surv_adult Numeric vector of aphid adult survivals by stage.
 @param repro Numeric vector of aphid reproductive rates by stage.

 */
//[[Rcpp::export]]
NumericMatrix leslie_matrix(IntegerVector instar_days, const double& surv_juv,
                            NumericVector surv_adult, NumericVector repro) {
    arma::mat LL;
    leslie_matrix__(as<arma::uvec>(instar_days), surv_juv,
                    as<arma::vec>(surv_adult),
                    as<arma::vec>(repro),
                    LL);
    return wrap(LL);
}
