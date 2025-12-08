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





/*
 * ============================================================================
 * Aphid class
 * ============================================================================
 */

class AphidPops {

    typedef std::normal_distribution<double> Norm;
    typedef std::binomial_distribution<uint32> Binom;
    typedef std::binomial_distribution<uint32>::param_type BinomParams;

    // For demographic stochasticity:
    Norm norm = Norm(0, 1);
    // For sampling numbers of alate leaving plant:
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
    double pseudo_surv;     // survival from Pseudomonas infection
    double extinct_N;  // extinction threshold
    bool demog_error;// whether to include demographic stochasticity
    double sigma_x; // environmental stochasticity
    double a;       // parasitoid attack rate
    double h;       // parasitoid handling time
    double k;       // parasitoid aggregation parameter
    double alate_infl; // inflection point for Pr(alates) ~ aphid density
    double alate_slope; // slope for Pr(alates) ~ aphid density
    double alate_max;   // max Pr(alates), should be 1 unless for testing
    double fly_p;   // probability of an alate leaving patch each day


    /*
     Non-dispersal of `iterate`, updates `N`, `W` without
     producing alates.
     Returns number of new adult (both male and female) parasitoids
     */
    double no_disp_iterate(const double& Yi, pcg32& eng) {

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
        if (Yi > 0) {
            A = arma::pow(1 + R * a * Yi / (k * (h * x + 1)), -k);
        } else A = arma::ones(R.n_elem);

        // Proportion of new aphids (from apterous females) that are alates
        double alate_p = alate_max / (1 + std::pow(10, ((alate_infl - z) * alate_slope)));

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

        X = (pred_surv * pseudo_surv * S * A) % LX;

        bool extinct_X = (arma::accu(X) + P) < extinct_N;
        if (extinct_X) {
            X.zeros();
            P = 0;
        }
        if (M < extinct_N) M = 0;


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

        return new_Y;

    }



public:

    AphidPops(const double& surv_j,
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
              const double& extinct_N_,
              const bool& demog_error_,
              const double& sigma_x_,
              const double& a_,
              const double& h_,
              const double& k_,
              const double& alate_infl_,
              const double& alate_slope_,
              const double& alate_max_,
              const double& fly_p_,
              const double& N0,
              const double& W0)
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
          extinct_N(extinct_N_),
          demog_error(demog_error_),
          sigma_x(sigma_x_),
          a(a_),
          h(h_),
          k(k_),
          alate_infl(alate_infl_),
          alate_slope(alate_slope_),
          alate_max(alate_max_),
          fly_p(fly_p_),
          X(4, arma::fill::zeros),
          P(0),
          M(0) {

        // Checks for parameter values that could cause negative numbers:
        if (trans_ma >= pred_surv) throw std::runtime_error("trans_ma cannot be >= pred_surv");
        if (trans_pm >= pred_surv * pseudo_surv * 0.65 * s_p) {
            std::string err_msg = "trans_pm cannot be >= pred_surv * pseudo_surv * 0.65 ";
            err_msg += "* s_p (0.65 is about as low as S_p(z(t)) goes)";
            throw std::runtime_error(err_msg.c_str());
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

    AphidPops(const AphidPops& other)
        : distr(other.distr), L(other.L), K(other.K), K_p(other.K_p),
          s_p(other.s_p), R(other.R), trans_ma(other.trans_ma),
          trans_pm(other.trans_pm), pred_surv(other.pred_surv), pseudo_surv(other.pseudo_surv),
          extinct_N(other.extinct_N),
          demog_error(other.demog_error), sigma_x(other.sigma_x),
          a(other.a), h(other.h), k(other.k),
          alate_infl(other.alate_infl), alate_slope(other.alate_slope),
          alate_max(other.alate_max), fly_p(other.fly_p),
          X(other.X), P(other.P), M(other.M) {};


    // iterate and set number of alates and add mummies transitioning to adult
    // parasitoids from this patch to the adult parasitoid population.
    // Note that `Yi` is the number of parasitoids on this patch only.
    void iterate(const double& Yi, uint32& n_alates, double& new_Ys, pcg32& eng) {

        new_Ys += no_disp_iterate(Yi, eng);

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
    void iterate(const double& Yi, double& new_Ys, pcg32& eng) {
        new_Ys += no_disp_iterate(Yi, eng);
        return;
    }

    /*
     Set pseudo_surv, which is useful for using a single AphidPops to populate a vector
     of them, then going back and changing some values of `pseudo_surv` based on
     the landscape:
     */
    void set_pseudo_surv(const double& new_pseudo_surv) {
        if (trans_pm >= pred_surv * pseudo_surv * 0.65 * s_p) {
            std::string err_msg = "trans_pm cannot be >= pred_surv * pseudo_surv * 0.65 ";
            err_msg += "* s_p (0.65 is about as low as S_p(z(t)) goes)";
            throw std::runtime_error(err_msg.c_str());
        }
        pseudo_surv = new_pseudo_surv;
        return;
    }
    /*
     Set K, which is useful for using a single AphidPops to populate a vector
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
    // Total alive aphids:
    double z() const {
        double out = P + arma::accu(X);
        return out;
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
};





/*
 * ============================================================================
 * Adult wasp class
 * ============================================================================
 */


class AdultWaspPop {

    double s_y;         // parasitoid adult daily survival
    double zeta;        // constant between 0 and 1 that affects the extent to
                        // which parasitoids respond to aphid density
    double extinct_N;   // extinction threshold
    arma::mat z_mat;    // Total aphid abundances by patch

public:

    AdultWaspPop(const double& s_y_,
                 const double& zeta_,
                 const double& extinct_N_,
                 const double& Y0)
    : s_y(s_y_),
      zeta(zeta_),
      extinct_N(extinct_N_),
      z_mat(),
      Y(Y0),
      Yi_mat() {};

    AdultWaspPop(const AdultWaspPop& other)
        : s_y(other.s_y),
          zeta(other.zeta),
          extinct_N(other.extinct_N),
          z_mat(other.z_mat),
          Y(other.Y),
          Yi_mat(other.Yi_mat) {};

    /*
     Fill vector of parasitoid abundances by patch
     Note: have to define this as template because otherwise it'd need to
           be after OnePlant definition.
           This works bc fill_Yi is used inside `PlantScape` class
           in file `plantscape.hpp` that includes `one-plant.hpp` before
           the class is defined.
     */
    template <typename P>
    void fill_Yi(const std::vector<std::vector<P>>& plants,
                 const arma::mat& wasp_attract) {
        uint32 n_x = plants.size();
        if (n_x == 0) throw std::runtime_error("ERROR: EMPTY PLANT");
        if (n_x != wasp_attract.n_rows)
            throw std::runtime_error("ERROR: `plants` DIMS DON'T MATCH `wasp_attract`");
        uint32 n_y = plants[0].size();
        if (n_y == 0) throw std::runtime_error("ERROR: EMPTY PLANT ROW");
        if (n_y != wasp_attract.n_cols)
            throw std::runtime_error("ERROR: `plants` DIMS DON'T MATCH `wasp_attract`");

        if (z_mat.n_rows != n_x || z_mat.n_cols != n_y) z_mat.set_size(n_x, n_y);
        if (Yi_mat.n_rows != n_x || Yi_mat.n_cols != n_y) Yi_mat.set_size(n_x, n_y);

        double z_tot = 0;
        for (uint32 x = 0; x < n_x; x++) {
            if (plants[x].size() != n_y)
                throw std::runtime_error("ERROR: INCONSISTENT `plants` VECTOR");
            for (uint32 y = 0; y < n_y; y++) {
                const AphidPops& aphids_xy(plants[x][y].aphids);
                z_mat(x,y) = aphids_xy.z();
                z_tot += z_mat(x,y);
            }
        }

        // Now go back through and calculate Y for each patch:
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                z_mat(x,y) /= z_tot;
                // Note: `wasp_attract` sums to 1 (verified inside `PlantScape`
                // constructor), and by default is the same for all patches
                Yi_mat(x,y) = Y * ((1-zeta) * wasp_attract(x,y) + zeta * z_mat(x,y));
            }
        }

        return;
    }

    void iterate(const double& new_Y) {
        Y = s_y * Y + 0.5 * new_Y;
        if (Y < extinct_N) Y = 0;
        return;
    }

    double Y;                   // Total adult, female parasitoid abundance
    arma::mat Yi_mat;            // Adult, female parasitoid abundance by patch

};




/*
 * ============================================================================
 * Wrapper class to include both aphids and wasps
 * ============================================================================
 */

struct InsectPops {

    AphidPops aphids;
    AdultWaspPop wasps;

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
               const double& extinct_N_,
               const bool& demog_error_,
               const double& sigma_x_,
               const double& a_,
               const double& h_,
               const double& k_,
               const double& alate_infl_,
               const double& alate_slope_,
               const double& alate_max_,
               const double& fly_p_,
               const double& N0,
               const double& W0,
               const double& s_y_,
               const double& zeta_,
               const double& Y0)
        : aphids(surv_j, surv_a, recruit, fecund, K_, K_p_mult, s_p_, R_,
                 trans_ma_, trans_pm_, pred_surv_, pseudo_surv_,
                 extinct_N_, demog_error_, sigma_x_, a_, h_, k_,
                 alate_infl_, alate_slope_, alate_max_, fly_p_, N0, W0),
          wasps(s_y_, zeta_, extinct_N_, Y0) {};


    // Iterate aphid and wasp populations. For use in `test_insect_pops`
    void iterate(pcg32& eng) {
        double new_Ys = 0;
        aphids.iterate(wasps.Y, new_Ys, eng);
        wasps.iterate(new_Ys);
        return;
    }

};


#endif
