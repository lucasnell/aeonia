#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <random>               // normal distribution
#include <pcg/pcg_random.hpp>   // pcg prng
#include "aeonia_types.hpp"     // integer types
#include "aphids.hpp"           // aphid classes
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;





// Add process error
void AphidTypePop::process_error(const arma::vec& Xt,
                                 const double& sigma_x,
                                 const double& rho,
                                 const bool& demog_error,
                                 const double& aphids_sum,
                                 std::normal_distribution<double>& norm_distr,
                                 pcg32& eng) {

    if (!demog_error && sigma_x == 0) return;

    uint32 n_stages = X.n_elem;

    // Variance for all process error:
    double sigma2 = sigma_x * sigma_x;
    if (demog_error) sigma2 += std::min(0.5, 1 / (1 + aphids_sum));

    arma::mat Se = sigma2 *
        (rho * arma::mat(n_stages, n_stages, arma::fill::ones) +
        (1 - rho) * arma::mat(n_stages, n_stages, arma::fill::eye));

    // chol doesn't work with zeros on diagonal
    arma::uvec non_zero = arma::find(Se.diag() > 0);

    /*
     Cholesky decomposition of Se so output has correct variance-covariance
     matrix:
     "a vector of independent normal random variables,
     when multiplied by the transpose of the Cholesky deposition of [Se] will
     have covariance matrix equal to [Se]."
     */
    arma::mat chol_decomp = arma::chol(Se(non_zero,non_zero)).t();

    // Random numbers from distribution N(0,1)
    arma::vec E(non_zero.n_elem);
    for (uint32 i = 0; i < E.n_elem; i++) E(i) = norm_distr(eng);

    // Making each element of E have correct variance-covariance matrix
    E = chol_decomp * E;

    // Plugging in errors into the X[t+1] vector
    for (uint32 i = 0; i < non_zero.n_elem; i++) {
        X(non_zero(i)) *= std::exp(E(i));
    }

    /*
     Because we used normal distributions to approximate demographic and environmental
     stochasticity, it is possible for aphids and parasitoids to
     "spontaneously appear" when the estimate of e(t) is large. To disallow this
     possibility, the number of aphids and parasitized aphids in a given age class
     on day t was not allowed to exceed the number in the preceding age class on
     day t – 1.
     */
    for (uint32 i = 1; i < n_stages; i++) {
        if (X(i) > Xt(i-1)) X(i) = Xt(i-1);
    }

    return;

}






/*
 Non-dispersal of `iterate`, updates densities without moving alates.
 Returns number of newly mummified aphids.
 */
double AphidPop::no_disp_iterate(const arma::vec& A_surv, pcg32& eng) {

    double new_mummies = 0; // newly mummified

    double z = total();

    if (z <= 0) return new_mummies;

    double S = 1 / (1 + z / K);
    double S_y = 1 / (1 + z / K_p);

    arma::vec A_mumm = 1 - A_surv;

    // making adult alates not able to be parasitized:
    arma::vec A_surv_ala = A_surv;
    arma::vec A_mumm_ala = A_mumm;
    for (uint32 i = adult_age; i < alates.leslie_.n_cols; i++) {
        // if (i >= A_surv_ala.n_elem) RcppThread::Rcout << "i too long" << std::endl;
        // else if (i >= A_mumm_ala.n_elem) RcppThread::Rcout << "i too long 2" << std::endl;
        // else {
        A_surv_ala[i] = 1;
        A_mumm_ala[i] = 0;
        // }
    }

    // Starting abundances (used in `process_error`):
    arma::vec apterous_Xt, alates_Xt, paras_Xt;
    bool do_error = demog_error || sigma_x > 0;
    if (do_error) {
        apterous_Xt = apterous.X;
        alates_Xt = alates.X;
        paras_Xt = paras.X;
    }

    // Basic updates for non-parasitized aphids:
    arma::mat LX_apt = apterous.leslie_ * apterous.X;
    arma::mat LX_ala = alates.leslie_ * alates.X;
    apterous.X = (pseudo_surv * pred_surv * S * A_surv) % LX_apt;
    alates.X = (pseudo_surv * pred_surv * S * A_surv_ala) % LX_ala;

    double np = 0; // newly parasitized
    np += pseudo_surv * pred_surv * S_y * arma::as_scalar(A_mumm.t() * LX_apt);
    np += pseudo_surv * pred_surv * S_y * arma::as_scalar(A_mumm_ala.t() * LX_ala);

    new_mummies += pred_surv * paras.X.back();  // newly mummified

    // alive but parasitized
    if (paras.X.n_elem > 1) {
        for (uint32 i = paras.X.n_elem - 1; i > 0; i--) {
            paras.X(i) = pseudo_surv * pred_surv * paras.s(i) * S_y * paras.X(i-1);
        }
    }
    paras.X.front() = np;


    // Process error
    if (do_error) {

        // Total aphids for this line:
        double aphids_sum = arma::accu(apterous_Xt) + arma::accu(alates_Xt) +
            arma::accu(paras_Xt);

        apterous.process_error(apterous_Xt, sigma_x, rho, demog_error, aphids_sum,
                               norm, eng);
        alates.process_error(alates_Xt, sigma_x, rho, demog_error, aphids_sum,
                             norm, eng);
        paras.process_error(paras_Xt, sigma_x, rho, demog_error, aphids_sum,
                            norm, eng);
    }

    // # offspring from apterous aphids that are alates:
    double alate_prop = apterous.alate_prop(z);
    double new_alates = apterous.X.front() * alate_prop;

    /*
     All offspring from alates are assumed to be apterous,
     so the only way to get new alates is from apterous aphids.
     */
    apterous.X.front() -= new_alates;
    apterous.X.front() += alates.X.front(); // <-- we assume alates make apterous
    alates.X.front() = new_alates;

    return new_mummies;
}



// Adjust starting abundances
void AphidPop::refresh_abunds(double N0, double W0) {

    double X_0_sum;
    if (N0 > 0) {
        X_0_sum = arma::accu(apterous.X_0_);
        if (X_0_sum <= 0) stop("apterous X_0 sums to <= 0");
        if (X_0_sum != 1) N0 /= X_0_sum;
        apterous.X_0_ *= N0;
    } else apterous.X_0_.zeros();

    if (W0 > 0) {
        X_0_sum = arma::accu(alates.X_0_);
        if (X_0_sum <= 0) stop("alates X_0 sums to <= 0");
        if (X_0_sum != 1) W0 /= X_0_sum;
        alates.X_0_ *= W0;
    } else alates.X_0_.zeros();

    // refresh starting conditions:
    apterous.X = apterous.X_0_;
    alates.X = alates.X_0_;
    paras.X.zeros(); // << because parasitized always start at zero abundance

    return;
}


/*
 Update new aphid abundances, update # alates leaving,
 return the # newly mummified aphids
 */
double AphidPop::iterate(arma::uvec& n_alates,
                         const arma::vec& A_surv,
                         pcg32& eng) {

    double new_mummies = no_disp_iterate(A_surv, eng);

    if (n_alates.n_elem != (alates.X.n_elem - adult_age))
        n_alates.set_size(alates.X.n_elem - adult_age);

    if (fly_p > 0) {
        uint32 binom_n;
        for (uint32 i = adult_age; i < alates.X.n_elem; i++) {
            uint32& n_alates_i(n_alates(i - adult_age));
            n_alates_i = 0U;
            binom_n = std::floor(alates.X(i));
            if (binom_n > 0U) {
                binom.param(BinomParams(binom_n, fly_p));
                n_alates_i = binom(eng);
                alates.X(i) -= static_cast<double>(n_alates_i);
            }
        }
    }

    return new_mummies;

}
// Overloaded for not doing any dispersing (used in `test_insect_pops`)
double AphidPop::iterate(const arma::vec& A_surv,
                         pcg32& eng) {
    double new_mummies = no_disp_iterate(A_surv, eng);
    return new_mummies;
}






//[[Rcpp::export]]
SEXP make_aphids_ptr(const double& K,
                    const double& K_p,
                    const double& pseudo_surv,
                    const double& pred_surv,
                    const double& extinct_N,
                    const bool& demog_error,
                    const double& sigma_x,
                    const double& rho,
                    const double& fly_p,
                    const arma::vec& attack_surv,
                    const std::string& aphid_name,
                    const arma::cube& leslie_mat,
                    const arma::mat& aphid_density_0,
                    const double& alate_b0,
                    const double& alate_b1,
                    const uint32& adult_age,
                    const uint32& living_days) {

    XPtr<AphidPop> aphid_xptr(new AphidPop(K, K_p, pseudo_surv, pred_surv,
                                           extinct_N, demog_error, sigma_x,
                                           rho, fly_p, attack_surv,
                                           aphid_name, leslie_mat,
                                           aphid_density_0, alate_b0, alate_b1,
                                           adult_age, living_days), true);

    return aphid_xptr;

}


//[[Rcpp::export]]
uint32 get_aphid_n_age_stages(SEXP aphid_ptr) {
    XPtr<AphidPop> aphid_xptr(aphid_ptr);
    return aphid_xptr->n_age_stages();
}



//' Calculate the proportion of offspring that are alates.
//'
//' @param z Numeric vector of total aphid abundances (including non-winged,
//'     winged, and parasitized).
//' @inheritParams make_insects_ptr
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector alate_prop(const NumericVector& z,
                         const double& alate_b0 = -5,
                         const double& alate_b1 = 0.0022) {

    if (is_true(any(z < 0.0))) stop("any(z < 0)");

    uint32 n = z.size();
    NumericVector out(n);

    for (uint32 i = 0; i < n; i++) {
        out[i] = alate_prop_cpp(z[i], alate_b0, alate_b1);
    }

    return out;

}
