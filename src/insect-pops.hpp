#ifndef __AEONIA_APHID_POPULATION_H
#define __AEONIA_APHID_POPULATION_H

/*
 This contains code for aphid and parasitoid population dynamics in a single patch.
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
    double K_p;     // parasitized aphid density dependence
    double s_p;     // parasitized aphid daily survival
    arma::vec R;    // relative parasitoid attack rates by aphid stage
    double trans_ma; // proportion of mummies that transition to adults each day
    double trans_pm; // proportion of parasitized aphids that transition to mummies each day
    double pred_surv;       // aphid survival from generalist predators and host plant death
    double pseudo_surv;       // survival from Pseudomonas infection
    double disaster_p; // probability of disaster
    double disaster_s; // disaster survival
    double extinct_N;  // extinction threshold
    bool demog_error;// whether to include demographic stochasticity
    double sigma_x; // environmental stochasticity
    double a;       // parasitoid attack rate
    double h;       // parasitoid handling time
    double k;       // parasitoid aggregation parameter
    double s_y;     // parasitoid adult daily survival
    double alate_0; // intercept for Pr(alates) ~ log(aphid density)
    double alate_1; // slope for Pr(alates) ~ log(aphid density)
    double fly_p;   // probability of an alate leaving patch each day
    double wasp_disp_m0;// proportion of adult parasitoids added to dispersal pool each day when no aphids present
    double wasp_disp_m1;// effect of aphid density on parasitoid emigration


    // Non-dispersal of `iterate`, updates `N`, `W`, and `Y` without
    // producing alates:
    void no_disp_iterate(pcg32& eng) {

        // max adults there could be (used for if stochasticity is added):
        arma::vec max_adults = {arma::accu(X.head(2)), arma::accu(X.tail(2))};

        // Total aphids:
        double x = arma::accu(X);
        double z = x + P;
        // Density dependence:
        double S = 1 / (1 + z / K);
        double S_p = 1 / (1 + z / K_p);
        // survival from parasitoids:
        arma::vec A;
        if (Y > 0) {
            A = arma::pow(1 + R * a * Y / (k * (h * x + 1)), -k);
        } else A = arma::ones(R.n_elem);

        // Proportion of new aphids (from apterous females) that are alates
        double alate_p = inv_logit__(alate_0 + alate_1 * z);
        // Now adjust Leslie matrix for alate proportion:
        L(0,1) = L(0,3) * (1 - alate_p); // non-winged offspring from winged adults
        L(2,1) = L(0,3) * alate_p; // winged offspring from non-winged adults
        // Note: because all offspring from winged adults are non-winged,
        // L(0,3) is always the "full-strength" fecundity

        arma::vec LX = L * X;

        // new parasitized aphids:
        double new_P = arma::as_scalar((1 - A).t() * LX);
        // new mummies:
        double new_M = trans_pm * P;
        // new adult parasitoids (both male and female!):
        double new_Y = trans_ma * M;

        P = pred_surv * pseudo_surv * S_p * (s_p * P + new_P) - new_M;
        M = pred_surv * (M + new_M) - new_Y;
        Y = s_y * Y + 0.5 * new_Y;

        X = (pred_surv * pseudo_surv * S * A) % LX;

        bool extinct_X = (arma::accu(X) + P) < extinct_N;
        if (extinct_X) {
            X.zeros();
            P = 0;
        }
        if (M < extinct_N) M = 0;
        if (Y < extinct_N) Y = 0;


        // Variance for all process error:
        if (!extinct_X && (sigma_x > 0 || demog_error)) {
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
        if (!extinct_X && (disaster_p > 0 && runif_01(eng) < disaster_p)) {
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
               const double& K_p_mult,
               const double& s_p_,
               const arma::vec& R_,
               const double& trans_ma_,
               const double& trans_pm_,
               const double& pred_surv_,
               const double& pseudo_surv_,
               const double& disaster_p_,
               const double& disaster_s_,
               const double& extinct_N_,
               const bool& demog_error_,
               const double& sigma_x_,
               const double& a_,
               const double& h_,
               const double& k_,
               const double& s_y_,
               const double& alate_0_,
               const double& alate_1_,
               const double& fly_p_,
               const double& wasp_disp_m0_,
               const double& wasp_disp_m1_,
               const double& N0,
               const double& W0,
               const double& Y0)
        : distr(1, 0.5),
          L(4, 4, arma::fill::zeros),
          K(K_),
          K_p(K * K_p_mult),
          s_p(s_p_),
          R(R_),
          trans_ma(trans_ma_),
          trans_pm(trans_pm_),
          pred_surv(pred_surv_),
          pseudo_surv(pseudo_surv_),
          disaster_p(disaster_p_),
          disaster_s(disaster_s_),
          extinct_N(extinct_N_),
          demog_error(demog_error_),
          sigma_x(sigma_x_),
          a(a_),
          h(h_),
          k(k_),
          s_y(s_y_),
          alate_0(alate_0_),
          alate_1(alate_1_),
          fly_p(fly_p_),
          wasp_disp_m0(wasp_disp_m0_),
          wasp_disp_m1(wasp_disp_m1_),
          X(4, arma::fill::zeros),
          P(0),
          M(0),
          Y(Y0) {

        // Checks for parameter values that could cause negative numbers:
        if (trans_ma >= pred_surv) stop("trans_ma cannot be >= pred_surv");
        if (trans_pm >= pred_surv * pseudo_surv * 0.65 * s_p) {
            std::string err_msg = "trans_pm cannot be >= pred_surv * pseudo_surv * 0.65 ";
            err_msg += "* s_p (0.65 is about as low as S_p(z(t)) goes)";
            stop(err_msg.c_str());
        }

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

        set_aphids(N0, W0);

    };

    InsectPops(const InsectPops& other)
        : distr(other.distr), L(other.L), K(other.K), K_p(other.K_p),
          s_p(other.s_p), R(other.R), trans_ma(other.trans_ma),
          trans_pm(other.trans_pm), pred_surv(other.pred_surv), pseudo_surv(other.pseudo_surv),
          disaster_p(other.disaster_p), disaster_s(other.disaster_s),
          extinct_N(other.extinct_N),
          demog_error(other.demog_error), sigma_x(other.sigma_x),
          a(other.a), h(other.h), k(other.k), s_y(other.s_y),
          alate_0(other.alate_0), alate_1(other.alate_1), fly_p(other.fly_p),
          wasp_disp_m0(other.wasp_disp_m0), wasp_disp_m1(other.wasp_disp_m1),
          X(other.X), P(other.P), M(other.M), Y(other.Y) {};


    // iterate and set number of alates and add parasitoids moving from this
    // population to the parasitoid dispersal pool:
    void iterate(uint32& n_alates, double& wasp_disp_pool, pcg32& eng) {

        no_disp_iterate(eng);

        if (wasp_disp_m0 > 0) {
            double p_out;
            if (wasp_disp_m1 != 0) {
                double lz = std::log(arma::accu(X));
                p_out = wasp_disp_m0 * std::exp(-wasp_disp_m1 * lz);
                if (p_out > 1) p_out = 1; // this can happen when accu(X) < 1
            } else {
                p_out = wasp_disp_m0;
            }
            wasp_disp_pool += Y * p_out;
            Y *= (1 - p_out);
        };

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
     Set pseudo_surv, which is useful for using a single InsectPops to populate a vector
     of them, then going back and changing some values of `pseudo_surv` based on
     the landscape:
     */
    void set_pseudo_surv(const double& new_pseudo_surv) {
        if (trans_pm >= pred_surv * pseudo_surv * 0.65 * s_p) {
            std::string err_msg = "trans_pm cannot be >= pred_surv * pseudo_surv * 0.65 ";
            err_msg += "* s_p (0.65 is about as low as S_p(z(t)) goes)";
            stop(err_msg.c_str());
        }
        pseudo_surv = new_pseudo_surv;
        return;
    }
    /*
     Set K, which is useful for using a single InsectPops to populate a vector
     of them, then going back and changing some values of `K` if you want
     the landscape to have variation:
     */
    void set_K(const double& new_K) {
        K_p /= K;
        K = new_K;
        K_p *= new_K;
        return;
    }


    // Get non-winged aphid population density:
    double aphids() const {
        return arma::accu(X.head(2));
    }
    // Get winged aphid population density:
    double alates() const {
        return arma::accu(X.tail(2));
    }
    double& alate_adults() {
        return X.back();
    }

    // Fill non-winged and winged aphids using stable age distributions,
    // for each separately:
    void set_aphids(const double& N0, const double& W0) {
        if (N0 <= 0 && W0 <= 0) return;
        arma::vec dens;
        sad_leslie__(L, dens);
        if (N0 > 0) {
            X.head(2) = N0 * dens.head(2) / arma::accu(dens.head(2));
        }
        if (W0 > 0) {
            X.tail(2) = W0 * dens.tail(2) / arma::accu(dens.tail(2));
        }
        return;
    }


    arma::vec X;    // aphids by stage
    double P;       // parasitized aphid density
    double M;       // mummy density
    double Y;       // parasitoid population density

};




#endif
