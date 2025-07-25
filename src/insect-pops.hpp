#ifndef __AEONIA_APHID_POPULATION_H
#define __AEONIA_APHID_POPULATION_H

/*
 This contains code for aphid and natural enemy population dynamics in a single patch.
 */

#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "aeonia_types.hpp"     // integer types
#include "convert-dims.hpp"     // XY
#include "math.hpp"              // inv_logit__ and sad_leslie__ fxns
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;





class InsectPops {

    typedef std::normal_distribution<double> Norm;
    typedef std::binomial_distribution<uint32> Binom;
    typedef std::binomial_distribution<uint32>::param_type BinomParams;

    Norm norm = Norm(0, 1);
    Binom distr;

    // constants:
    arma::mat L;    // transition matrix
    double K;       // aphid density dependence
    double B;       // effect of Pseudomonas bacteria on aphid growth
    double disaster_p; // probability of disaster
    double disaster_s; // disaster survival
    double extinct_N;  // extinction threshold
    bool demog_error;// whether to include demographic stochasticity
    double sigma_x; // environmental stochasticity
    double a;       // natural enemy attack rate
    double h;       // natural enemy handling time
    double k;       // natural enemy aggregation parameter
    double s;       // natural enemy daily survival
    double alate_0; // intercept for Pr(alates) ~ log(aphid density)
    double alate_1; // intercept for Pr(alates) ~ log(aphid density)
    double fly_p;   // probability of an alate leaving patch each day
    double wasp_d_p;// proportion of adult parasitoids added to dispersal pool each day


    // Non-dispersal of `iterate`, updates `A`, `W`, and `P` without
    // producing alates:
    void no_disp_iterate(pcg32& eng) {

        // max adults there could be (used for if stochasticity is added):
        arma::vec max_adults = {arma::accu(X.head(2)), arma::accu(X.tail(2))};

        // Total aphids:
        double z = arma::accu(X);
        // Density dependence:
        double S = 1 / (1 + z / K);
        // survival from natural enemies:
        double attack_surv;
        if (P > 0) attack_surv = std::pow(1 + a * P / (k * (h * z + 1)), -k);
        else attack_surv = 1;

        // Proportion of new aphids (from apterous females) that are alates
        double alate_p = inv_logit__(alate_0 + alate_1 * z);
        // Now adjust Leslie matrix for alate proportion:
        L(0,1) = L(0,3) * (1 - alate_p); // non-winged offspring from winged adults
        L(2,1) = L(0,3) * alate_p; // winged offspring from non-winged adults
        // Note: because all offspring from winged adults are non-winged,
        // L(0,3) is always the "full-strength" fecundity

        arma::vec LX = L * X;

        if (P > 0) {
            P *= s; // daily survival of existing
            P += arma::accu((1 - B) * S * (1 - attack_surv) * LX); // new individuals
        }

        X = (1 - B) * S * attack_surv * LX;


        // Variance for all process error:
        if (sigma_x > 0 || demog_error) {
            double sigma2 = 0;
            if (demog_error) sigma2 += std::min(0.5, 1 / (1 + arma::accu(X)));
            if (sigma_x > 0) sigma2 += (sigma_x * sigma_x);
            double sigma = std::sqrt(sigma2);
            for (double& x : X) x *= std::exp(norm(eng) * sigma);
            // Don't allow adults to exceed all stages from previous time point:
            if (X(1) > max_adults(0)) X(1) = max_adults(0);
            if (X(3) > max_adults(1)) X(3) = max_adults(1);
        }

        // Sample for disaster:
        if (disaster_p > 0 && runif_01(eng) < disaster_p) {
            X *= disaster_s;
        }

        return;

    }



public:

    InsectPops(const double& surv_j,
               const double& surv_a,
               const double& recruit,
               const double& fecund,
               const double& K_,
               const double& B_,
               const double& disaster_p_,
               const double& disaster_s_,
               const double& extinct_N_,
               const bool& demog_error_,
               const double& sigma_x_,
               const double& a_,
               const double& h_,
               const double& k_,
               const double& s_,
               const double& alate_0_,
               const double& alate_1_,
               const double& fly_p_,
               const double& wasp_d_p_,
               const double& A0,
               const double& W0,
               const double& P0)
        : distr(1, 0.5),
          L(4, 4, arma::fill::zeros),
          K(K_),
          B(B_),
          disaster_p(disaster_p_),
          disaster_s(disaster_s_),
          extinct_N(extinct_N_),
          demog_error(demog_error_),
          sigma_x(sigma_x_),
          a(a_),
          h(h_),
          k(k_),
          s(s_),
          alate_0(alate_0_),
          alate_1(alate_1_),
          fly_p(fly_p_),
          wasp_d_p(wasp_d_p_),
          X(4, arma::fill::zeros),
          P(P0) {

        L(0,0) = surv_j;
        L(2,2) = surv_j;
        L(1,1) = surv_a;
        L(3,3) = surv_a;
        /*
         Note: changing these next two lines doesn't change the beginning
         densities further below because I'm separately specifying winged
         vs non-winged aphid abundance.
         */
        L(0,1) = fecund * 0.5;
        L(2,1) = fecund * 0.5;
        L(0,3) = fecund;
        L(1,0) = recruit;
        L(3,2) = recruit;

        set_aphids(A0, W0);

    };

    InsectPops(const InsectPops& other)
        : distr(other.distr), L(other.L), K(other.K), B(other.B),
          disaster_p(other.disaster_p), disaster_s(other.disaster_s),
          extinct_N(other.extinct_N),
          demog_error(other.demog_error), sigma_x(other.sigma_x),
          a(other.a), h(other.h), k(other.k), s(other.s),
          alate_0(other.alate_0), alate_1(other.alate_1), fly_p(other.fly_p),
          wasp_d_p(other.wasp_d_p),
          X(other.X), P(other.P) {};


    // iterate and set number of alates and parasitoids moving from this population:
    void iterate(uint32& n_alates, double& n_wasp_d, pcg32& eng) {

        no_disp_iterate(eng);

        if (wasp_d_p > 0) {
            n_wasp_d = P * wasp_d_p;
            P -= n_wasp_d;
        } else n_wasp_d = 0;

        n_alates = 0;
        double& adult_winged(X(3));
        uint32 binom_n = std::floor(adult_winged);
        if (binom_n > 0 && fly_p > 0) {
            distr.param(BinomParams(binom_n, fly_p));
            n_alates = distr(eng);
            adult_winged -= static_cast<double>(n_alates);
        }

        return;
    }
    // Overloaded for not using output object (used in `test_insect_pops`)
    void iterate(pcg32& eng) {
        no_disp_iterate(eng);
        return;
    }

    /*
     Set B, which is useful for using a single InsectPops to populate a vector
     of them, then going back and changing some values of `B` based on
     the landscape:
     */
    void set_B(const double& new_B) {
        B = new_B;
        return;
    }
    /*
     Set K, which is useful for using a single InsectPops to populate a vector
     of them, then going back and changing some values of `K` if you want
     the landscape to have variation:
     */
    void set_K(const double& new_K) {
        K = new_K;
        return;
    }


    // Get non-winged aphid population density:
    double A() const {
        return arma::accu(X.head(2));
    }
    // Get winged aphid population density:
    double W() const {
        return arma::accu(X.tail(2));
    }
    double& winged_adults() {
        return X.back();
    }

    // Fill non-winged and winged aphids using stable age distributions,
    // for each separately:
    void set_aphids(const double& A0, const double& W0) {
        if (A0 <= 0 && W0 <= 0) return;
        arma::vec dens;
        sad_leslie__(L, dens);
        if (A0 > 0) {
            X.head(2) = A0 * dens.head(2) / arma::accu(dens.head(2));
        }
        if (W0 > 0) {
            X.tail(2) = W0 * dens.tail(2) / arma::accu(dens.tail(2));
        }
        return;
    }


    arma::vec X;    // aphids by stage
    double P;       // natural enemy population density

};




#endif
