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

inline double logit__(double& p) {
    a = std::log(p / (1-p));
    return a;
}
inline double inv_logit__(double& a) {
    double p = 1 / (1 + std::exp(-a));
    return p;
}








#endif

