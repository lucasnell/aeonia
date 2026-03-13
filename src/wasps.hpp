# ifndef __AEONIA_WASP_POPULATION_H
# define __AEONIA_WASP_POPULATION_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <cmath>                // log, exp
#include <random>               // normal distribution
#include <cstdint>              // integer types
#include <algorithm>            // find
#include <pcg/pcg_random.hpp>   // pcg prng
#include "aeonia_types.hpp"  // integer types



using namespace Rcpp;




/*
 ==============================================================================
 ==============================================================================
 Mummy population
 ==============================================================================
 ==============================================================================
 */

class MummyPop {

    arma::vec M_0;          // initial mummy densities
    /*
     Proportion of mummies that will NOT take exactly 3 days to develop.
     As this value approaches 2/3, it will provide greater smoothing of
     wasp numbers through time.
     */
    double smooth;

    // Survival from predators
    double pred_surv;

    // Extinction threshold
    double extinct_N;

public:



    // Changing through time
    arma::vec M;            // Mummy density

    // Constructors
    inline MummyPop()
        : M_0(4, arma::fill::zeros),
          smooth(0),
          pred_surv(1),
          M(4, arma::fill::zeros) {};
    inline MummyPop(const double& smooth_,
                    const double& pred_surv_,
                    const double& extinct_N_)
        : M_0(4, arma::fill::zeros),
          smooth(smooth_),
          pred_surv(pred_surv_),
          extinct_N(extinct_N_),
          M(4, arma::fill::zeros) {};
    inline MummyPop(const arma::vec& M_0_,
                    const double& smooth_,
                    const double& pred_surv_,
                    const double& extinct_N_)
        : M_0(arma::join_vert(arma::vec(1, arma::fill::zeros), M_0_)),
          smooth(smooth_),
          pred_surv(pred_surv_),
          extinct_N(extinct_N_),
          M(arma::join_vert(arma::vec(1, arma::fill::zeros), M_0_)) {};
    inline MummyPop(const MummyPop& other)
        : M_0(other.M_0), smooth(other.smooth), pred_surv(other.pred_surv),
          extinct_N(other.extinct_N),
          M(other.M) {};
    inline MummyPop& operator=(const MummyPop& other) {
        M_0 = other.M_0;
        smooth = other.smooth;
        pred_surv = other.pred_surv;
        extinct_N = other.extinct_N;
        M = other.M;
        return *this;
    }


    // Update # mummies
    // `extinct` is # aphids that are newly "mummified"
    // returns # of new adult parasitoids (both male and female)
    double iterate(const double& new_mummies);

    // Clearing a field kills all mummies
    inline void clear() {
        M.fill(0);
    }
    // Clearing part of field
    inline void clear(const double& surv) {
        M *= surv;
    }

    inline void refresh_abunds(double mummy_M0) {
        double m0_sum = arma::accu(M_0);
        if (m0_sum != 1 && m0_sum > 0) mummy_M0 /= m0_sum;
        M_0 *= mummy_M0;
        M = M_0;
        return;
    }

    inline double total() const {
        return arma::accu(M);
    }


};







/*
 ==============================================================================
 ==============================================================================
 Wasp attack info
 ==============================================================================
 ==============================================================================
 */

class WaspAttack {

public:

    arma::vec rel_attack;    // relative wasp attack rates by aphid stage
    double a;                // overall parasitoid attack rate
    double k;                // aggregation parameter of the nbinom distribution
    double h;                // parasitoid attack rate handling time


    // Constructors
    WaspAttack()
        : rel_attack(), a(), k(), h() {};
    WaspAttack(const arma::vec& rel_attack_,
               const double& a_,
               const double& k_,
               const double& h_)
        : rel_attack(rel_attack_),
          a(a_),
          k(k_),
          h(h_) {};
    WaspAttack(const WaspAttack& other)
        : rel_attack(other.rel_attack),
          a(other.a),
          k(other.k),
          h(other.h) {};
    WaspAttack& operator=(const WaspAttack& other) {
        rel_attack = other.rel_attack;
        a = other.a;
        k = other.k;
        h = other.h;
        return *this;
    }

    /*
     Compute probabilities of aphids surviving wasp attacks and of
     aphids being successfully mummified.
     These don't necessarily add to 1 bc of the possibility of mutual mortality,
     which is especially common for superparasitism.

     Survivals are equation 6 from Meisner et al. (2014), where
     `rel_attack` is equivalent to p_i
     */
    void A_mats(const double& Y_m, const double& x, arma::vec& A_surv,
                const arma::vec& attack_surv) const;



};




/*
 ==============================================================================
 ==============================================================================
 Wasp population
 ==============================================================================
 ==============================================================================
 */



// Adult wasp and mummy population
class WaspPop {

    typedef std::normal_distribution<double> Norm;
    typedef std::binomial_distribution<uint32> Binom;
    typedef std::binomial_distribution<uint32>::param_type BinomParams;

    // For environmental stochasticity:
    Norm norm = Norm(0, 1);
    // For demographic stochasticity:
    Binom binom = Binom(1, 0.5);

public:

    WaspAttack attack;  // info for attack rates
    double Y_0;         // initial adult wasp density
    uint32 delay;       // when to add initial wasps
    double sex_ratio;   // proportion of female wasps
    double s_y;         // parasitoid adult daily survival
    double zeta;        // wasp response to aphid density
    double sigma_y;     // for process error
    bool demog_error;   // whether to include demographic stochasticity

    // Changing through time
    double Y;               // Wasp density


    // Constructors
    inline WaspPop()
        : attack(), Y_0(), delay(), sex_ratio(), s_y(), zeta(), sigma_y(),
          demog_error(), Y() {};
    inline WaspPop(const arma::vec& rel_attack_,
            const double& a_,
            const double& k_,
            const double& h_,
            const double& Y_0_,
            const uint32& delay_,
            const double& sex_ratio_,
            const double& s_y_,
            const double& zeta_,
            const double& sigma_y_,
            const bool& demog_error_)
        : attack(rel_attack_, a_, k_, h_),
          Y_0(Y_0_),
          delay(delay_),
          sex_ratio(sex_ratio_),
          s_y(s_y_),
          zeta(zeta_),
          sigma_y(sigma_y_),
          demog_error(demog_error_),
          Y((delay_ == 0) ? Y_0_ : 0.0) {};
    inline WaspPop(const WaspPop& other)
        : attack(other.attack),
          Y_0(other.Y_0),
          delay(other.delay),
          sex_ratio(other.sex_ratio),
          s_y(other.s_y),
          zeta(other.zeta),
          sigma_y(other.sigma_y),
          demog_error(other.demog_error),
          Y(other.Y) {};
    inline WaspPop& operator=(const WaspPop& other) {
        attack = other.attack;
        Y_0 = other.Y_0;
        delay = other.delay;
        sex_ratio = other.sex_ratio;
        s_y = other.s_y;
        zeta = other.zeta;
        sigma_y = other.sigma_y;
        demog_error = other.demog_error;
        Y = other.Y;
        return *this;
    }


    // Adjust starting abundances
    inline void refresh_abunds(const double& adult_Y0) {
        Y_0 = adult_Y0;
        // refresh starting conditions:
        if (delay == 0) {
            Y = Y_0;
        } else Y = 0;
        return;
    }


    // Fills matrix for Pr(aphids survive)
    inline void A_mats(arma::vec& A_surv,
                       const arma::vec& attack_surv,
                       const double& Yi,
                       const double& x) const {
        attack.A_mats(Yi, x, A_surv, attack_surv);
        return;
    }

    // Check for whether to add initial wasps:
    inline void add_Y_0_check(const uint32& t) {
        if (t == delay) Y += Y_0;
        return;
    }

    /*
     Update # adult wasps
     `old_mums` is # mummies at time t that are in the last mummy stage
     before merging
     */
    void iterate(const double& old_mums, pcg32& eng);


};









#endif
