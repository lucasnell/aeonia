
/*
 This contains code for aphid and predator population dynamics in a single patch.
 */

#include <RcppArmadillo.h>
#include <vector>


#include "aeonia_types.hpp"     // integer types
#include "insect-pops.hpp"
#include "pcg.hpp"  // pcg types
#include "util.hpp"  // retrieve_dataset function


using namespace Rcpp;



void check_insect_args(const uint32& max_t,
                       const double& N0,
                       const double& W0,
                       const double& Y0,
                       const double& B,
                       const double& disaster_p,
                       const double& disaster_s,
                       const double& extinct_N,
                       const double& sigma_x,
                       const double& surv_j,
                       const double& surv_a,
                       const double& recruit,
                       const double& fecund,
                       const double& K,
                       const double& K_p,
                       const double& s_p,
                       const arma::vec& R_,
                       const double& theta_m,
                       const double& theta_p,
                       const double& m,
                       const double& alate_1,
                       const double& a,
                       const double& h,
                       const double& k,
                       const double& s_y,
                       const double& fly_p,
                       const double& wasp_disp_m0) {

    if (max_t > (uint32)1e9) stop("max_t > 1e9");
    if (N0 < 0) stop("N0 < 0");
    if (W0 < 0) stop("W0 < 0");
    if (Y0 < 0) stop("Y0 < 0");
    if (B < 0 || B > 1) stop("B < 0 || B > 1");
    if (disaster_p < 0 || disaster_p > 1) stop("disaster_p < 0 || disaster_p > 1");
    if (disaster_s < 0 || disaster_s > 1) stop("disaster_s < 0 || disaster_s > 1");
    if (extinct_N < 0) stop("extinct_N < 0");
    if (sigma_x < 0) stop("sigma_x < 0");
    if (surv_j <= 0 || surv_j > 1) stop("surv_j <= 0 || surv_j > 1");
    if (surv_a <= 0 || surv_a > 1) stop("surv_a <= 0 || surv_a > 1");
    if (recruit <= 0 || recruit > 1) stop("recruit <= 0 || recruit > 1");
    if (fecund <= 0) stop("fecund <= 0");
    if (K <= 0) stop("K <= 0");
    if (K_p <= 0) stop("K_p <= 0");
    if (s_p < 0 || s_p > 1) stop("s_p < 0 || s_p > 1");
    if (arma::any(R_ < 0) || arma::any(R_ > 1)) stop("any(R < 0 | R > 1)");
    if (R_.n_elem != 4) stop("length(R) != 4");
    if (theta_m < 0 || theta_m > 1) stop("theta_m < 0 || theta_m > 1");
    if (theta_p < 0 || theta_p > 1) stop("theta_p < 0 || theta_p > 1");
    if (m < 0 || m >= 1) stop("m < 0 || m >= 1");
    if (alate_1 < 0) stop("alate_1 < 0");
    if (a < 0) stop("a < 0");
    if (h < 0) stop("h < 0");
    if (k < 0) stop("k < 0");
    if (s_y < 0 || s_y > 1) stop("s_y < 0 || s_y > 1");
    if (fly_p < 0 || fly_p > 1) stop("fly_p < 0 || fly_p > 1");
    if (wasp_disp_m0 < 0 || wasp_disp_m0 > 1) stop("wasp_disp_m0 < 0 || wasp_disp_m0 > 1");

    return;

}


/*
 Fill parameters from `aeonia::pop_info` if they're not provided.
 */
void fill_pop_info(double& surv_j,
                   double& surv_a,
                   double& recruit,
                   double& fecund,
                   double& K,
                   double& K_p,
                   double& s_p,
                   Rcpp::NumericVector& R,
                   double& theta_m,
                   double& theta_p,
                   double& m,
                   double& alate_0,
                   double& alate_1,
                   double& a,
                   double& h,
                   double& k,
                   double& s_y) {

    List pop_info = retrieve_dataset<List>("pop_info");

    if (NumericVector::is_na(surv_j)) surv_j = pop_info["surv_j"];
    if (NumericVector::is_na(surv_a)) surv_a = pop_info["surv_a"];
    if (NumericVector::is_na(recruit)) recruit = pop_info["recruit"];
    if (NumericVector::is_na(fecund)) fecund = pop_info["fecund"];
    if (NumericVector::is_na(K)) K = pop_info["K"];
    if (NumericVector::is_na(K_p)) K_p = pop_info["K_p"];
    if (NumericVector::is_na(s_p)) s_p = pop_info["s_p"];
    if (R.size() == 0) R = pop_info["R"];
    if (NumericVector::is_na(theta_m)) theta_m = pop_info["theta_m"];
    if (NumericVector::is_na(theta_p)) theta_p = pop_info["theta_p"];
    if (NumericVector::is_na(m)) m = pop_info["m"];
    if (NumericVector::is_na(alate_0)) alate_0 = pop_info["alate_0"];
    if (NumericVector::is_na(alate_1)) alate_1 = pop_info["alate_1"];
    if (NumericVector::is_na(a)) a = pop_info["a"];
    if (NumericVector::is_na(h)) h = pop_info["h"];
    if (NumericVector::is_na(k)) k = pop_info["k"];
    if (NumericVector::is_na(s_y)) s_y = pop_info["s_y"];

    return;
}




//' Create a pointer object in which to store insect population info.
//'
//' This pointer is used as an argument to [sim_plantscape()].
//'
//'
//' @param B Single numeric indicating the effect of *Pseudomonas* on aphid
//'     population growth.
//' @param fly_p Single numeric indicating the proportion of alates that fly
//'     off plants each day.
//' @param wasp_disp_m0 Proportion of adult wasps from each field that
//'     are added to the dispersal pool when there are no aphids present.
//'     Defaults to `0`.
//' @param wasp_disp_m1 Effect of aphid density on wasp emigration from a patch.
//'     Emigration is `wasp_disp_m0 * exp(-wasp_disp_m1 * log(z))`, where `z` is
//'     the total number of living aphids in the patch.
//'     Defaults to `0`.
//' @param disaster_p Single numeric indicating the probability of disaster
//'     each day. Defaults to `0`.
//' @param disaster_s Single numeric indicating disaster survival.
//'     Defaults to `0`.
//' @param disaster_p Single numeric indicating the probability of disaster
//'     each day. Defaults to `0`.
//' @param extinct_N Single numeric indicating the extinction threshold.
//'     Defaults to `0`.
//' @param demog_error Single logical for whether to include demographic
//'     stochasticity. Defaults to `FALSE`.
//' @param sigma_x Single numeric indicating the standard deviation
//'     for the lognormal distribution used to generate environmental
//'     stochasticity.
//'     Defaults to `0`.
//' @param surv_j Single numeric indicating aphid juvenile survival.
//'     Defaults to `NA`, which results in `pop_info$surv_j` being used.
//'     This is from a previous study.
//' @param surv_a Single numeric indicating aphid adult survival.
//'     Defaults to `NA`, which results in `pop_info$surv_a` being used.
//'     This is from a previous study.
//' @param recruit Single numeric indicating aphid recruitment.
//'     Defaults to `NA`, which results in `pop_info$recruit` being used.
//'     This is from a previous study.
//' @param fecund Single numeric indicating aphid fecundity.
//'     Defaults to `NA`, which results in `pop_info$fecund` being used.
//'     This is from a previous study.
//' @param K Single numeric indicating the carrying capacity (not including
//'     parasitoids or *Pseudomonas*) of the aphid population.
//'     Defaults to `NA`, which results in `pop_info$K` being used.
//'     This is from a previous study.
//' @param m Single numeric indicating aphid mortality.
//'     Defaults to `NA`, which results in `pop_info$m` being used.
//'     This is from a previous study.
//' @param alate_0 Single numeric.
//'     The proportion of winged offspring from apterous aphids is
//'     `inv_logit(alate_0 + alate_1 * z)` where `z` is the total number of
//'     aphids on that plant.
//'     Defaults to `NA`, which results in `pop_info$alate_0` being used.
//'     This is from a previous study.
//' @param alate_1 Single numeric affecting how strongly aphid density
//'     influences alate production. See `alate_0` above for the equation.
//'     Defaults to `NA`, which results in `pop_info$alate_1` being used.
//'     This is from a previous study.
//' @param a Single numeric indicating the parasitoid attack rate.
//' @param h Single numeric indicating the parasitoid handling time.
//' @param k Single numeric indicating the parasitoid aggregation parameter.
//' @param s_y Single numeric indicating the adult parasitoid daily survival.
//'
//' @export
//'
//[[Rcpp::export]]
SEXP make_insect_ptr(const double& B,
                     const double& fly_p,
                     const double& wasp_disp_m0 = 0,
                     const double& wasp_disp_m1 = 0,
                     const double& disaster_p = 0,
                     const double& disaster_s = 0,
                     const double& extinct_N = 0,
                     const bool& demog_error = false,
                     const double& sigma_x = 0,
                     double surv_j = NA_REAL,
                     double surv_a = NA_REAL,
                     double recruit = NA_REAL,
                     double fecund = NA_REAL,
                     double K = NA_REAL,
                     double K_p = NA_REAL,
                     double s_p = NA_REAL,
                     NumericVector R = NumericVector::create(),
                     double theta_m = NA_REAL,
                     double theta_p = NA_REAL,
                     double m = NA_REAL,
                     double alate_0 = NA_REAL,
                     double alate_1 = NA_REAL,
                     double a = NA_REAL,
                     double h = NA_REAL,
                     double k = NA_REAL,
                     double s_y = NA_REAL) {

    fill_pop_info(surv_j, surv_a, recruit, fecund, K, K_p, s_p, R,
                  theta_m, theta_p, m, alate_0, alate_1, a, h, k, s_y);

    arma::vec R_ = Rcpp::as<arma::vec>(R);

    // Densities will be set later:
    double N0 = 0;
    double W0 = 0;
    double Y0 = 0;
    uint32 max_t = 100;

    check_insect_args(max_t, N0, W0, Y0, B, disaster_p, disaster_s,
                      extinct_N, sigma_x, surv_j, surv_a, recruit, fecund,
                      K, K_p, s_p, R_, theta_m, theta_p, m, alate_1,
                      a, h, k, s_y,
                      fly_p, wasp_disp_m0);

    XPtr<InsectPops> insect_xptr(new InsectPops(surv_j, surv_a, recruit, fecund,
                                                K, K_p, s_p, R_,
                                                theta_m, theta_p, m, B,
                                                disaster_p, disaster_s, extinct_N,
                                                demog_error, sigma_x,
                                                a, h, k, s_y,
                                                alate_0, alate_1, fly_p,
                                                wasp_disp_m0, wasp_disp_m1,
                                                N0, W0, Y0), true);

    return insect_xptr;
}





//' Test population dynamics for insects for a set of parameters.
//'
//' This function is just for one plant, so the probability that an
//' alate leaves the plant (argument `fly_p` in [make_insect_ptr()]) is always 1.
//'
//' @inheritParams make_insect_ptr
//' @inheritParams sim_plantscape
//'
//' @export
//'
//[[Rcpp::export]]
DataFrame test_insect_pops(const uint32& max_t,
                           const double& N0,
                           const double& W0,
                           const double& Y0,
                           const double& B = 0,
                           const double& disaster_p = 0,
                           const double& disaster_s = 0,
                           const double& extinct_N = 0,
                           const bool& demog_error = false,
                           const double& sigma_x = 0,
                           double surv_j = NA_REAL,
                           double surv_a = NA_REAL,
                           double recruit = NA_REAL,
                           double fecund = NA_REAL,
                           double K = NA_REAL,
                           double K_p = NA_REAL,
                           double s_p = NA_REAL,
                           NumericVector R = NumericVector::create(),
                           double theta_m = NA_REAL,
                           double theta_p = NA_REAL,
                           double m = NA_REAL,
                           double alate_0 = NA_REAL,
                           double alate_1 = NA_REAL,
                           double a = NA_REAL,
                           double h = NA_REAL,
                           double k = NA_REAL,
                           double s_y = NA_REAL) {

    fill_pop_info(surv_j, surv_a, recruit, fecund, K, K_p, s_p, R,
                  theta_m, theta_p, m, alate_0, alate_1, a, h, k, s_y);

    arma::vec R_ = Rcpp::as<arma::vec>(R);

    double fly_p = 0;
    double wasp_disp_m0 = 0;
    double wasp_disp_m1 = 0;

    check_insect_args(max_t, N0, W0, Y0, B, disaster_p, disaster_s,
                      extinct_N, sigma_x, surv_j, surv_a, recruit, fecund,
                      K, K_p, s_p, R_, theta_m, theta_p, m, alate_1,
                      a, h, k, s_y, fly_p, wasp_disp_m0);

    InsectPops insects(surv_j, surv_a, recruit, fecund, K, K_p, s_p, R_,
                       theta_m, theta_p, m, B,
                       disaster_p, disaster_s, extinct_N, demog_error, sigma_x,
                       a, h, k, s_y,
                       alate_0, alate_1, fly_p, wasp_disp_m0, wasp_disp_m1,
                       N0, W0, Y0);

    std::vector<uint32> time;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> parasitized;
    std::vector<double> mummies;
    std::vector<double> wasps;
    time.reserve(max_t+1);
    aphids.reserve(max_t+1);
    alates.reserve(max_t+1);
    parasitized.reserve(max_t+1);
    mummies.reserve(max_t+1);
    wasps.reserve(max_t+1);
    time.push_back(0);
    aphids.push_back(N0);
    alates.push_back(W0);
    wasps.push_back(Y0);

    pcg32 eng;
    seed_pcg(eng);

    for (uint32 t = 0; t < max_t; t++) {
        insects.iterate(eng);
        time.push_back(t+1);
        aphids.push_back(insects.aphids());
        alates.push_back(insects.alates());
        parasitized.push_back(insects.P);
        mummies.push_back(insects.M);
        wasps.push_back(insects.Y);
    }

    DataFrame out_df = DataFrame::create(_["time"] = time,
                                         _["aphids"] = aphids,
                                         _["alates"] = alates,
                                         _["parasitized"] = parasitized,
                                         _["mummies"] = mummies,
                                         _["wasps"] = wasps);

    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out_df;

}


