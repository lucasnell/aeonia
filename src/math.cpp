#include <RcppArmadillo.h>
#include <cmath>

#include "math.hpp"
#include "aeonia_types.hpp"

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
