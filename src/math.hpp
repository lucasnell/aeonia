# ifndef __AEONIA_MATH_H
# define __AEONIA_MATH_H


#include <RcppArmadillo.h>
#include <cmath>
#include <random>

#include <pcg/pcg_random.hpp>   // pcg prng
#include "aeonia_types.hpp"
#include "pcg.hpp"              // runif_ab fxn


using namespace Rcpp;








/*
 =====================================================================================
 =====================================================================================
 Logit and inverse logit
 =====================================================================================
 =====================================================================================
 */

inline void logit__(const double& p, double& out) {
    out = std::log(p / (1-p));
    return;
}
inline void inv_logit__(const double& a, double& out) {
    out = 1 / (1 + std::exp(-a));
    return;
}

inline double logit__(const double& p) {
    double a = std::log(p / (1-p));
    return a;
}
inline double inv_logit__(const double& a) {
    double p = 1 / (1 + std::exp(-a));
    return p;
}





/*
 =====================================================================================
 =====================================================================================
 Leslie matrix
 =====================================================================================
 =====================================================================================
 */

// This computes the "stable age distribution" from the transition matrix, which is
// the proportion of different classes that is required for the population to grow
// exponentially
inline void sad_leslie__(const arma::mat& L, arma::vec& out) {

    arma::cx_vec r_cx;
    arma::cx_mat SAD;

    arma::eig_gen(r_cx, SAD, L);

    arma::vec r = arma::abs(r_cx);

    double rmax = arma::max(r);

    arma::cx_mat SADdist = SAD.cols(arma::find(r == rmax));
    arma::cx_double all_SAD = arma::accu(SADdist);
    SADdist /= all_SAD;
    SADdist.resize(SADdist.n_elem, 1);

    out = arma::real(SADdist);

    return;
}






#endif

