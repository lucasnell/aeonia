#include <RcppArmadillo.h>
#include <cmath>

#include "math.hpp"
#include "aeonia-types.hpp"

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

