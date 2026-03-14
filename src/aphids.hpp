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




inline double alate_prop_cpp(const double& z,
                             const double& alate_b0,
                             const double& alate_b1) {
    const double lap = alate_b0 + alate_b1 * z;
    double ap;
    inv_logit__(lap, ap);
    return ap;
}




// Necessary here to declare friendship
class AphidPop;


// Generic aphid "type" population: either apterous or alate for a particular clonal line
class AphidTypePop {

protected:

    arma::mat leslie_;       // Leslie matrix with survival and reproduction
    arma::vec X_0_;          // initial aphid abundances by stage


public:
    // Changing through time
    arma::vec X;                // Aphid density

    /*
     Constructors
     */
    inline AphidTypePop() : leslie_(), X_0_(), X() {};
    inline AphidTypePop(const arma::mat& leslie_mat,
                        const arma::vec& aphid_density_0)
        : leslie_(leslie_mat),
          X_0_(aphid_density_0),
          X(aphid_density_0) {};

    inline AphidTypePop(const AphidTypePop& other)
        : leslie_(other.leslie_),
          X_0_(other.X_0_),
          X(other.X) {};

    inline AphidTypePop& operator=(const AphidTypePop& other) {
        leslie_ = other.leslie_;
        X_0_ = other.X_0_;
        X = other.X;
        return *this;
    }


    /*
     Total aphids
     */
    inline double total() const {
        return arma::accu(X);
    }

    // Kill all aphids
    inline void clear() {
        X.fill(0);
        return;
    }
    // Kill some of aphids
    inline void clear(const double& surv) {
        X *= surv;
        return;
    }

    // Add process error:
    void process_error(const arma::vec& Xt,
                       const double& sigma_x,
                       const double& rho,
                       const bool& demog_error,
                       const double& aphids_sum,
                       std::normal_distribution<double>& norm_distr,
                       pcg32& eng);

    // Returning references to private members:
    inline const arma::mat& leslie() const {return leslie_;}
    inline const arma::vec& X_0() const {return X_0_;}


};



// Aphid "type" population for apterous of a particular clonal line
class ApterousPop : public AphidTypePop {

    friend class AphidPop;

public:

    // Parameters for logit(Pr(alates)) ~ b0 + b1 * N
    double alate_b0_;
    double alate_b1_;


    inline ApterousPop() : AphidTypePop(), alate_b0_(0), alate_b1_(0) {};
    inline ApterousPop(const arma::mat& leslie_mat,
                       const arma::vec& aphid_density_0,
                       const double& alate_b0,
                       const double& alate_b1)
        : AphidTypePop(leslie_mat, aphid_density_0),
          alate_b0_(alate_b0),
          alate_b1_(alate_b1){};
    inline ApterousPop(const ApterousPop& other)
        : AphidTypePop(other),
          alate_b0_(other.alate_b0_),
          alate_b1_(other.alate_b1_){};
    inline ApterousPop& operator=(const ApterousPop& other) {
        AphidTypePop::operator=(other);
        alate_b0_ = other.alate_b0_;
        alate_b1_ = other.alate_b1_;
        return *this;
    }


    // logit(Pr(alates)) ~ b0 + b1 * z, where `z` is # aphids (all lines)
    inline double alate_prop(const double& z) const {
        double ap = alate_prop_cpp(z, alate_b0_, alate_b1_);
        return ap;
    }

};





// Aphid "type" population for alates of a particular clonal line
class AlatePop : public AphidTypePop {

    friend class AphidPop;


public:

    inline AlatePop() : AphidTypePop() {};
    inline AlatePop(const arma::mat& leslie_mat,
             const arma::vec& aphid_density_0)
        : AphidTypePop(leslie_mat, aphid_density_0) {};

    inline AlatePop(const AlatePop& other)
        : AphidTypePop(other) {};

    inline AlatePop& operator=(const AlatePop& other) {
        AphidTypePop::operator=(other);
        return *this;
    }


};

// Aphid "type" population for parasitized (but alive) aphids
class ParasitizedPop : public AphidTypePop {

    friend class AphidPop;

protected:

    arma::vec s;        // vector of survival rates of parasitized aphids by day

public:

    inline ParasitizedPop() : AphidTypePop(), s() {};
    inline ParasitizedPop(const arma::mat& leslie_mat,
                          const uint32& living_days)
        : AphidTypePop(arma::mat(), arma::vec(living_days, arma::fill::zeros)),
          s(arma::diagvec(leslie_mat, -1)) {
        s.resize(living_days);
    };
    inline ParasitizedPop(const ParasitizedPop& other)
        : AphidTypePop(other),
          s(other.s) {};

    inline ParasitizedPop& operator=(const ParasitizedPop& other) {
        AphidTypePop::operator=(other);
        s = other.s;
        return *this;
    }



};






// Aphid population: both alates and apterous for one clonal line in a field
class AphidPop {

    typedef std::normal_distribution<double> Norm;
    typedef std::binomial_distribution<uint32> Binom;
    typedef std::binomial_distribution<uint32>::param_type BinomParams;

    // For demographic stochasticity:
    Norm norm = Norm(0, 1);
    // For sampling numbers of alate leaving plant:
    Binom binom = Binom(1, 0.5);


    /*
     Non-dispersal of `iterate`, updates densities without moving alates.
     Returns number of newly mummified aphids.
     */
    double no_disp_iterate(const arma::vec& A_surv, pcg32& eng);


public:

    double K;           // aphid density dependence
    double K_p;         // parasitized aphid density dependence
    double pseudo_surv; // survival from Pseudomonas infection
    double pred_surv;   // aphid survival from generalist predators and host plant death
    double extinct_N;   // extinction threshold
    bool demog_error;   // whether to include demographic stochasticity
    double sigma_x;     // environmental stochasticity
    double rho;         // environmental correlation among instars
    double fly_p;       // probability of an adult alate leaving patch each day

    uint32 adult_age;           // age at which aphids are adults

    /*
     Vector of length >=2 with survival probabilities of singly & multiply
     attacked aphids:
     */
    arma::vec attack_surv;


    std::string aphid_name;    // unique identifying name for this aphid line
    ApterousPop apterous;
    AlatePop alates;
    ParasitizedPop paras;
    bool extinct;



    /*
     Constructors.
     */
    inline AphidPop()
        : K(0),
          K_p(0),
          pseudo_surv(0),
          pred_surv(0),
          extinct_N(0),
          demog_error(0),
          sigma_x(0),
          rho(0),
          fly_p(0),
          adult_age(0),
          attack_surv(arma::vec()),
          aphid_name("NULL"),
          apterous(),
          alates(),
          paras(),
          extinct(true) {};

    // Make sure `leslie_mat` has 3 slices and `aphid_density_0` has two columns!
    inline AphidPop(const double& K_,
                    const double& K_p_,
                    const double& pseudo_surv_,
                    const double& pred_surv_,
                    const double& extinct_N_,
                    const bool& demog_error_,
                    const double& sigma_x_,
                    const double& rho_,
                    const double& fly_p_,
                    const arma::vec& attack_surv_,
                    const std::string& aphid_name_,
                    const arma::cube& leslie_mat,
                    const arma::mat& aphid_density_0,
                    const double& alate_b0,
                    const double& alate_b1,
                    const uint32& adult_age_,
                    const uint32& living_days)
        : K(K_),
          K_p(K_p_),
          pseudo_surv(pseudo_surv_),
          pred_surv(pred_surv_),
          extinct_N(extinct_N_),
          demog_error(demog_error_),
          sigma_x(sigma_x_),
          rho(rho_),
          fly_p(fly_p_),
          adult_age(adult_age_),
          attack_surv(attack_surv_),
          aphid_name(aphid_name_),
          apterous(leslie_mat.slice(0), aphid_density_0.col(0), alate_b0, alate_b1),
          alates(leslie_mat.slice(1), aphid_density_0.col(1)),
          paras(leslie_mat.slice(2), living_days),
          extinct(false) {};

    inline AphidPop(const AphidPop& other)
        : K(other.K),
          K_p(other.K_p),
          pseudo_surv(other.pseudo_surv),
          pred_surv(other.pred_surv),
          extinct_N(other.extinct_N),
          demog_error(other.demog_error),
          sigma_x(other.sigma_x),
          rho(other.rho),
          fly_p(other.fly_p),
          adult_age(other.adult_age),
          attack_surv(other.attack_surv),
          aphid_name(other.aphid_name),
          apterous(other.apterous),
          alates(other.alates),
          paras(other.paras),
          extinct(other.extinct) {};

    inline AphidPop& operator=(const AphidPop& other) {

        K = other.K;
        K_p = other.K_p;
        pseudo_surv = other.pseudo_surv;
        pred_surv = other.pred_surv;
        extinct_N = other.extinct_N;
        demog_error = other.demog_error;
        sigma_x = other.sigma_x;
        rho = other.rho;
        fly_p = other.fly_p;
        adult_age = other.adult_age;
        attack_surv = other.attack_surv;
        aphid_name = other.aphid_name;
        apterous = other.apterous;
        alates = other.alates;
        paras = other.paras;
        extinct = other.extinct;

        return *this;

    };


    /*
     Total aphids
     */
    inline double total() const {
        double ta = apterous.total() + alates.total() +
            paras.total();
        return ta;
    }
    // Total unparasitized
    inline double total_unpar() const {
        double ta = apterous.total() + alates.total();
        return ta;
    }



    // Kill all aphids
    inline void clear() {
        apterous.clear();
        alates.clear();
        paras.clear();
        extinct = true;
        return;
    }
    // Kill some aphids
    inline void clear(const double& surv) {
        apterous.clear(surv);
        alates.clear(surv);
        paras.clear(surv);
        return;
    }

    // Adjust starting abundances
    void refresh_abunds(double N0, double W0);

    /*
     Update new aphid abundances, update # alates leaving,
     return the # newly mummified aphids
     */
    double iterate(arma::uvec& n_alates, const arma::vec& A_surv, pcg32& eng);
    // Overloaded for not doing any dispersing (used in `test_insect_pops`)
    double iterate(const arma::vec& A_surv, pcg32& eng);

    // /*
    //  Adult and juvenile numbers for alates and apterous:
    //  */
    // inline double total_adult_apterous() const {
    //     return arma::accu(apterous.X.tail(apterous.X.n_elem - adult_age));
    // }
    // inline double total_juven_apterous() const {
    //     return arma::accu(apterous.X.head(adult_age));
    // }
    // inline double total_adult_alates() const {
    //     return arma::accu(alates.X.tail(alates.X.n_elem - adult_age));
    // }
    // inline double total_juven_alates() const {
    //     return arma::accu(alates.X.head(adult_age));
    // }

    inline arma::vec unparas_X() const {
        return arma::join_cols(apterous.X, alates.X);
    }
    inline arma::vec unparas_X_0() const {
        return arma::join_cols(apterous.X_0(), alates.X_0());
    }

    inline uint32 n_stages() const {
        return apterous.X.n_elem + alates.X.n_elem;
    }
    inline uint32 n_age_stages() const {
        return apterous.X.n_elem;
    }


};






#endif
