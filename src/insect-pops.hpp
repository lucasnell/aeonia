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





class InsectPops {

    typedef std::binomial_distribution<uint32> Binom;
    typedef std::binomial_distribution<uint32>::param_type BinomParams;

    Binom distr;

    // constants:
    double r;       // aphid population growth rate
    double K;       // aphid carrying capacity (in absence of predators / Pseudomonas)
    double B;       // effect of Pseudomonas bacteria on aphid growth
    double pred_a;  // predator attack rate
    double pred_h;  // predator handling time
    double pred_c;  // consumption efficiency for predators
    double pred_m;  // predator mortality
    double alate_0; // intercept for Pr(alates) ~ log(aphid density)
    double alate_1; // intercept for Pr(alates) ~ log(aphid density)
    double fly_p;   // probability of an alate leaving patch each day

public:

    InsectPops(const double& r_,
               const double& K_,
               const double& B_,
               const double& pred_a_,
               const double& pred_h_,
               const double& pred_c_,
               const double& pred_m_,
               const double& alate_0_,
               const double& alate_1_,
               const double& fly_p_,
               const double& A0,
               const double& W0,
               const double& P0)
        : distr(1, 0.5),
          r(r_),
          K(K_),
          B(B_),
          pred_c(pred_c_),
          pred_a(pred_a_),
          pred_h(pred_h_),
          pred_m(pred_m_),
          alate_0(alate_0_),
          alate_1(alate_1_),
          fly_p(fly_p_),
          A(A0),
          W(W0),
          P(P0) {};

    InsectPops(const InsectPops& other)
    : distr(other.distr), r(other.r), K(other.K), B(other.B),
      pred_c(other.pred_c), pred_a(other.pred_a), pred_h(other.pred_h),
      pred_m(other.pred_m), alate_0(other.alate_0), alate_1(other.alate_1),
      fly_p(other.fly_p), A(other.A), W(other.W), P(other.P) {};


    // iterate and output number of alates moving from this population:
    uint32 iterate(pcg32& eng) {

        double z = A + W;
        // Proportion of new aphids (from apterous females) that are alates
        double alate_p = inv_logit__(alate_0 + alate_1 * z);
        // predator per-capita consumption per prey:
        double consumpt = pred_a / (1 + pred_a * pred_h * z);
        // proportional abundance change for both winged and non-winged aphids
        // (i.e., (X_t+1 - X_t) / X_t, where X is W and A):
        double pac = std::exp(r * (1 - z / K) - consumpt * P - B) - 1;
        double dA = A * pac;
        double dW = W * pac;
        if (pac > 0) {
            /*
             Next two lines are because alates produce only apterous aphids,
             so the only new alates come from apterous mothers.
             */
            dA += dW;
            dW = dA * alate_p; // << should be `=`, NOT `+=`
            // To balance number of new aphids:
            dA *= (1 - alate_p);
        }
        A += dA;
        W += dW;
        P *= std::exp(pred_c * consumpt * z - pred_m)

        uint32 n_alates = 0;
        uint32 binom_n = std::floor(W);
        if (binom_n > 0) {
            distr.param(BinomParams(binom_n, fly_p));
            n_alates = distr(eng);
            W -= static_cast<double>(n_alates);
        }

        return n_alates;
    }

    /*
     Set B, which is useful for using a single InsectPops to population
     everything, then going back and changing some values of `B` based on
     the landscape:
     */
    void set_B(const double& new_B) {
        B = new_B;
        return;
    }


    double A;  // non-winged aphid population density
    double W;  // winged aphid population density
    double P;  // predator population density

};




#endif
