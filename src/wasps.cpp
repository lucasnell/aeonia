
#include "wasps.hpp"            // wasp types



using namespace Rcpp;



// Update # mummies
// `extinct` is # aphids that are newly "mummified"
// returns # of new adult parasitoids (both male and female)
double MummyPop::iterate(const double& new_mummies) {

    double old_mums = M.back();

    // Go backwards through stages to avoid conflicts...
    for (uint32 i = M.n_elem-1; i > 0; i--) {
        M(i) = pred_surv * M(i-1);
    }

    // Add newly mummified over three days:
    if (smooth > 0) {
        M(0) = new_mummies * smooth / 2;
        M(1) += (new_mummies * (1 - smooth));
        M(2) += (new_mummies * smooth / 2);
    } else M(1) += new_mummies;

    if (arma::accu(M) < extinct_N) M.zeros();

    return old_mums;

}





/*
 Compute probabilities of aphids surviving wasp attacks and of
 aphids being successfully mummified.
 These don't necessarily add to 1 bc of the possibility of mutual mortality,
 which is especially common for superparasitism.

 Survivals are equation 6 from Meisner et al. (2014), where
 `rel_attack` is equivalent to p_i
 */
void WaspAttack::A_mats(const double& Y_m,
                        const double& x,
                        arma::vec& A_surv,
                        const arma::vec& attack_surv) const {

    uint32 n_max_attacks = attack_surv.n_elem;
    uint32 n_stages = rel_attack.n_elem;

    if (A_surv.n_elem != n_stages) A_surv.set_size(n_stages);

    if (Y_m == 0) {
        A_surv.ones();
        return;
    }

    // Mean of negative binomial distribution:
    arma::vec A_bar = rel_attack * (a * Y_m) / (h * x + 1);

    // Probabilities of being attacked by stage and number of attacks,
    // where the last probability is the prob of being attacked
    // **at least** `n_max_attacks` times:
    arma::mat attack_probs(n_stages, n_max_attacks+1U, arma::fill::none);

    // Probability of being attacked zero times:
    attack_probs.col(0) = arma::pow((1 + A_bar / k), -k);

    // Now we (optionally) expand to being attacked >0 times.
    // The last term will be 1 - (sum of all other probs) because it refers
    // to the prob of being attacked **at least** `n_max_attacks` times.
    if (n_max_attacks == 1) {
        attack_probs.col(1) = 1 - attack_probs.col(0);
    } else if (n_max_attacks > 1) {
        // Two terms that get used for all:
        arma::vec Aa = 1 + k / A_bar;
        arma::vec Ab = 1 + A_bar / k;
        double prod = 1;
        // Will sum all other cols into this vector in the loop below:
        arma::vec A_prob_sums = attack_probs.col(0);
        // Do for all in `attack_surv` and add to `A_`:
        for (uint32 i = 1; i < n_max_attacks; i++) {
            double nt = i; // number of times attacked
            prod *= ((k - 1) / nt + 1);
            attack_probs.col(i) = arma::pow(Aa, -nt) % arma::pow(Ab, -k);
            attack_probs.col(i) *= prod;
            // To avoid an extra loop (see below):
            A_prob_sums += attack_probs.col(i);
        }
        // Now fill in the last column which is 1 - sum(other cols):
        attack_probs.col(n_max_attacks) = (1 - A_prob_sums);
    }

    /*
     Converting attack probabilities to survivals and successful
     mummifications
     */
    for (uint32 j = 0; j < n_stages; j++) {
        A_surv(j) = attack_probs(j, 0);
        for (uint32 i = 0; i < n_max_attacks; i++) {
            A_surv(j) += attack_probs(j, i+1) * attack_surv(i);
        }
    }


    return;
}





/*
 Update # adult wasps
 `old_mums` is # mummies at time t that are in the last mummy stage
 before merging
 */
void WaspPop::iterate(const double& old_mums, pcg32& eng) {
    double max_Y = old_mums + Y;
    if (max_Y == 0) return;
    if (demog_error) {
        // Number of adult females that survives is binomial with `s_y`
        // as probability and number of adult females as number of trials.
        uint32 n = static_cast<uint32>(std::round(Y));
        std::binomial_distribution<uint32>::param_type prms(n, s_y);
        binom.param(prms);
        Y = static_cast<double>(binom(eng));
        // Number of mummies that are female is binomial with `sex_ratio`
        // as probability and number of old mummies as number of trials.
        n = static_cast<uint32>(std::round(old_mums));
        prms = std::binomial_distribution<uint32>::param_type(n, sex_ratio);
        binom.param(prms);
        Y += static_cast<double>(binom(eng));
    } else {
        Y *= s_y;
        Y += (sex_ratio * old_mums);
    }
    if (sigma_y > 0) {
        Y *= std::exp(norm(eng) * sigma_y);
        // make sure it doesn't exceed what's possible:
        if (Y > max_Y) Y = max_Y;
    }
    return;
}



/*
 ==============================================================================
 ==============================================================================
 Pointer makers
 ==============================================================================
 ==============================================================================
 */



//[[Rcpp::export]]
SEXP make_wasps_ptr(const arma::vec& rel_attack,
                    const double& a,
                    const double& k,
                    const double& h,
                    const double& sex_ratio,
                    const double& s_y,
                    const double& zeta,
                    const double& sigma_y,
                    const bool& demog_error) {

    // these will also be adjusted in the simulation fxn so are not included here
    uint32 delay = 0;
    double Y_0 = 0;

    XPtr<WaspPop> wasps_xptr(
            new WaspPop(rel_attack, a, k, h, Y_0, delay, sex_ratio, s_y, zeta,
                        sigma_y, demog_error), true);

    return wasps_xptr;

}


//[[Rcpp::export]]
SEXP make_mummies_ptr(const double& pred_surv,
                      const double& mummy_smooth,
                      const double& extinct_N,
                      const uint32& mummy_dev_time) {

    arma::vec M_0(mummy_dev_time, arma::fill::ones);

    XPtr<MummyPop> mummies_xptr(new MummyPop(M_0, mummy_smooth, pred_surv,
                                             extinct_N), true);

    return mummies_xptr;

}
