#ifndef __AEONIA_APHID_POPULATION_H
#define __AEONIA_APHID_POPULATION_H

/*
 This contains code for aphid and predator population dynamics in a single patch.
 */

#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "aeonia_types.hpp"     // integer types
#include "convert-dims.hpp"     // XY
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;





class Populations {

    typedef std::binomial_distribution<uint32> Binom;
    typedef std::binomial_distribution<uint32>::param_type BinomParams;

    Binom distr;

    // constants:
    double r;  // aphid population growth rate
    double a;  // aphid density dependence
    double B;  // effect of Pseudomonas bacteria on aphid growth
    double pred_0;  // intercept for predators ~ log10(aphid density)
    double pred_1;  // intercept for predators ~ log10(aphid density)
    double pred_c;  // "proportionality constant" for predators ~ log10(aphid density)
    double pred_k;  // power for effect of predators on aphids
    double alate_0; // intercept for Pr(alates) ~ log(aphid density)
    double alate_1; // intercept for Pr(alates) ~ log(aphid density)
    double fly_p;   // probability of an alate leaving patch each day

public:

    Populations(const double& r_,
                const double& a_,
                const double& B_,
                const double& pred_0_,
                const double& pred_1_,
                const double& pred_c_,
                const double& pred_k_,
                const double& alate_0_,
                const double& alate_1_,
                const double& fly_p_,
                const double& A0,
                const double& W0,
                const double& P0)
        : distr(1, 0.5),
          r(r_),
          a(a_),
          B(B_),
          pred_0(pred_0_),
          pred_1(pred_1_),
          pred_c(pred_c_),
          pred_k(pred_k_),
          alate_0(alate_0_),
          alate_1(alate_1_),
          fly_p(fly_p_),
          A(A0),
          W(W0),
          P(P0) {}

    // iterate and output number of alates moving from this patch:
    uint32 iterate(pcg32& eng) {

        double z = A + W;
        double alate_p = inv_logit__(alate_0 + alate_1 * z);
        // proportional abundance change for both winged and non-winged aphids
        // (i.e., (X_t+1 - X_t) / X_t, where X is W and A):
        double pac = std::exp(r - a * z - std::pow(P, pred_k) - B) - 1;
        double dA = A * pac;
        double dW = W * pac;
        if (pac > 0) {
            dW += (dA * alate_p);
            dA *= (1 - alate_p);
        }
        A += dA;
        W += dW;
        P = pred_c * std::exp(pred_0 + pred_1 * std::log10(z));

        uint32 n_alates = 0;
        uint32 distr_n = std::floor(W);
        if (distr_n > 0) {
            distr.param(BinomParams(distr_n, fly_p));
            n_alates = distr(eng);
            W -= static_cast<double>(n_alates);
        }

        return n_alates;
    }


    double A;  // non-winged aphid population density
    double W;  // winged aphid population density
    double P;  // predator population density

};




#endif
