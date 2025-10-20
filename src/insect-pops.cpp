
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
                       const double& pseudo_surv,
                       const double& extinct_N,
                       const double& sigma_x,
                       const double& surv_j,
                       const double& surv_a,
                       const double& recruit,
                       const double& fecund,
                       const double& K,
                       const double& K_p_mult,
                       const double& s_p,
                       const arma::vec& R_,
                       const double& trans_ma,
                       const double& trans_pm,
                       const double& pred_surv,
                       const double& alate_infl,
                       const double& alate_slope,
                       const double& a,
                       const double& h,
                       const double& k,
                       const double& s_y,
                       const double& fly_p,
                       const double& zeta) {

    if (max_t > (uint32)1e9) stop("max_t > 1e9");
    if (N0 < 0) stop("N0 < 0");
    if (W0 < 0) stop("W0 < 0");
    if (Y0 < 0) stop("Y0 < 0");
    if (pseudo_surv < 0 || pseudo_surv > 1) stop("pseudo_surv < 0 || pseudo_surv > 1");
    if (extinct_N < 0) stop("extinct_N < 0");
    if (sigma_x < 0) stop("sigma_x < 0");
    if (surv_j <= 0 || surv_j > 1) stop("surv_j <= 0 || surv_j > 1");
    if (surv_a <= 0 || surv_a > 1) stop("surv_a <= 0 || surv_a > 1");
    if (recruit <= 0 || recruit > 1) stop("recruit <= 0 || recruit > 1");
    if (fecund <= 0) stop("fecund <= 0");
    if (K <= 0) stop("K <= 0");
    if (K_p_mult <= 0) stop("K_p_mult <= 0");
    if (s_p < 0 || s_p > 1) stop("s_p < 0 || s_p > 1");
    if (arma::any(R_ < 0) || arma::any(R_ > 1)) stop("any(R < 0 | R > 1)");
    if (R_.n_elem != 4) stop("length(R) != 4");
    if (trans_ma < 0 || trans_ma > 1) stop("trans_ma < 0 || trans_ma > 1");
    if (trans_pm < 0 || trans_pm > 1) stop("trans_pm < 0 || trans_pm > 1");
    if (pred_surv < 0 || pred_surv >= 1) stop("pred_surv < 0 || pred_surv >= 1");
    if (alate_infl < 0) stop("alate_infl < 0");
    if (alate_slope < 0) stop("alate_slope < 0");
    // One of these can be true, but not both
    // (because it results in Inf * 0 which is NaN):
    if (alate_infl == arma::datum::inf && alate_slope == 0) {
        stop("alate_infl == Inf && alate_slope == 0");
    }
    if (a < 0) stop("a < 0");
    if (h < 0) stop("h < 0");
    if (k < 0) stop("k < 0");
    if (s_y < 0 || s_y > 1) stop("s_y < 0 || s_y > 1");
    if (fly_p < 0 || fly_p > 1) stop("fly_p < 0 || fly_p > 1");
    if (zeta < 0 || zeta > 1) stop("zeta < 0 || zeta > 1");

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
                   double& K_p_mult,
                   double& s_p,
                   Rcpp::NumericVector& R,
                   double& trans_ma,
                   double& trans_pm,
                   double& pred_surv,
                   double& alate_infl,
                   double& alate_slope,
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
    if (NumericVector::is_na(K_p_mult)) K_p_mult = pop_info["K_p_mult"];
    if (NumericVector::is_na(s_p)) s_p = pop_info["s_p"];
    if (R.size() == 0) R = pop_info["R"];
    if (NumericVector::is_na(trans_ma)) trans_ma = pop_info["trans_ma"];
    if (NumericVector::is_na(trans_pm)) trans_pm = pop_info["trans_pm"];
    if (NumericVector::is_na(pred_surv)) pred_surv = pop_info["pred_surv"];
    if (NumericVector::is_na(alate_infl)) alate_infl = pop_info["alate_infl"];
    if (NumericVector::is_na(alate_slope)) alate_slope = pop_info["alate_slope"];
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
//' For parameters that defaults to `NA`, keeping the default results in
//' `pop_info$Z` (where `Z` is the name of the variable) being used.
//' These parameters are mostly from previous studies, but some were
//' estimated for this simpler model from the full one.
//' See `data-raw/pop_info.R` for how these were generated.
//'
//'
//' @param pseudo_surv Single numeric indicating aphid survival from *Pseudomonas*.
//' @param fly_p Single numeric indicating the proportion of alates that fly
//'     off plants each day.
//' @param zeta Constant between 0 and 1 that affects the extent to which
//'     parasitoids respond to aphid density, where 0 results in an even
//'     distribution of parasitoids, and 1 results in a linear relationship
//'     between aphid density and parasitoids.
//'     Defaults to `0`.
//' @param extinct_N Single numeric indicating the extinction threshold.
//'     Defaults to `0`.
//' @param demog_error Single logical for whether to include demographic
//'     stochasticity. Defaults to `FALSE`.
//' @param sigma_x Single numeric indicating the standard deviation
//'     for the lognormal distribution used to generate environmental
//'     stochasticity.
//'     Defaults to `0`.
//' @param surv_j Single numeric indicating aphid juvenile survival.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param surv_a Single numeric indicating aphid adult survival.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param recruit Single numeric indicating aphid recruitment.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param fecund Single numeric indicating aphid fecundity.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param K Single numeric indicating the density dependence of the
//'     aphid population.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param K_p_mult Single numeric indicating the multiplier for parasitized
//'     aphid density dependence (`K_p = K * K_p_mult`).
//'     Defaults to `NA`. See 'Details' for more info.
//' @param s_p Single numeric indicating parasitized aphid daily survival.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param R Single numeric indicating
//'     Defaults to `NA`. See 'Details' for more info.
//' @param trans_ma Single numeric indicating the proportion of mummies that
//'     transition to adult parasitoids each day.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param trans_pm Single numeric indicating the proportion of parasitized
//'     aphids that transition to mummies each day.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param pred_surv Single numeric indicating aphid and mummy mortality due
//'     to generalist predators.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param alate_infl Single numeric for the inflection point for the sigmoid
//'     relationship between aphid density and alate offspring proportion.
//'     The proportion of winged offspring from apterous aphids is
//'     `1 / {1 + 10^((alate_infl - z) * alate_slope)}` where `z` is the total number of
//'     aphids on that plant.
//'     To turn off alate production, set this parameter to `Inf`.
//'     Note: Do NOT set both `alate_infl = Inf` and `alate_slope = 0`
//'     because that'll trigger an error that's put in place because
//'     setting these parameters in this way would result in `NaN` values.
//'     Must be > 0.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param alate_slope Single numeric for the slope for the sigmoid
//'     relationship between aphid density and alate offspring proportion.
//'     See `alate_infl` above for the equation.
//'     Must be > 0.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param a Single numeric indicating the parasitoid attack rate.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param h Single numeric indicating the parasitoid handling time.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param k Single numeric indicating the parasitoid aggregation parameter.
//'     Defaults to `NA`. See 'Details' for more info.
//' @param s_y Single numeric indicating the adult parasitoid daily survival.
//'     Defaults to `NA`. See 'Details' for more info.
//'
//' @export
//'
//' @return An `externalptr` object that points to a C++ object that can
//' be pass to [sim_plantscape()].
//'
//[[Rcpp::export]]
SEXP make_insect_ptr(const double& pseudo_surv,
                     const double& fly_p,
                     const double& zeta = 0,
                     const double& extinct_N = 0,
                     const bool& demog_error = false,
                     const double& sigma_x = 0,
                     double surv_j = NA_REAL,
                     double surv_a = NA_REAL,
                     double recruit = NA_REAL,
                     double fecund = NA_REAL,
                     double K = NA_REAL,
                     double K_p_mult = NA_REAL,
                     double s_p = NA_REAL,
                     NumericVector R = NumericVector::create(),
                     double trans_ma = NA_REAL,
                     double trans_pm = NA_REAL,
                     double pred_surv = NA_REAL,
                     double alate_infl = NA_REAL,
                     double alate_slope = NA_REAL,
                     double a = NA_REAL,
                     double h = NA_REAL,
                     double k = NA_REAL,
                     double s_y = NA_REAL) {

    fill_pop_info(surv_j, surv_a, recruit, fecund, K, K_p_mult, s_p, R,
                  trans_ma, trans_pm, pred_surv, alate_infl, alate_slope, a, h, k, s_y);

    arma::vec R_ = Rcpp::as<arma::vec>(R);

    // Densities will be set later:
    double N0 = 0;
    double W0 = 0;
    double Y0 = 0;
    uint32 max_t = 100;

    check_insect_args(max_t, N0, W0, Y0, pseudo_surv,
                      extinct_N, sigma_x, surv_j, surv_a, recruit, fecund,
                      K, K_p_mult, s_p, R_, trans_ma, trans_pm, pred_surv,
                      alate_infl, alate_slope,
                      a, h, k, s_y,
                      fly_p, zeta);

    XPtr<InsectPops> insect_xptr(new InsectPops(surv_j, surv_a, recruit, fecund,
                                                K, K_p_mult, s_p, R_,
                                                trans_ma, trans_pm,
                                                pred_surv, pseudo_surv,
                                                extinct_N, demog_error,
                                                sigma_x, a, h, k, alate_infl,
                                                alate_slope, fly_p, N0, W0,
                                                s_y, zeta, Y0), true);

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
                           const double& pseudo_surv = 1,
                           const double& extinct_N = 0,
                           const bool& demog_error = false,
                           const double& sigma_x = 0,
                           double surv_j = NA_REAL,
                           double surv_a = NA_REAL,
                           double recruit = NA_REAL,
                           double fecund = NA_REAL,
                           double K = NA_REAL,
                           double K_p_mult = NA_REAL,
                           double s_p = NA_REAL,
                           NumericVector R = NumericVector::create(),
                           double trans_ma = NA_REAL,
                           double trans_pm = NA_REAL,
                           double pred_surv = NA_REAL,
                           double alate_infl = NA_REAL,
                           double alate_slope = NA_REAL,
                           double a = NA_REAL,
                           double h = NA_REAL,
                           double k = NA_REAL,
                           double s_y = NA_REAL) {

    fill_pop_info(surv_j, surv_a, recruit, fecund, K, K_p_mult, s_p, R,
                  trans_ma, trans_pm, pred_surv, alate_infl, alate_slope, a, h, k, s_y);

    arma::vec R_ = Rcpp::as<arma::vec>(R);

    double fly_p = 0;
    double zeta = 0;

    check_insect_args(max_t, N0, W0, Y0, pseudo_surv,
                      extinct_N, sigma_x, surv_j, surv_a, recruit, fecund,
                      K, K_p_mult, s_p, R_, trans_ma, trans_pm, pred_surv,
                      alate_infl, alate_slope,
                      a, h, k, s_y, fly_p, zeta);

    InsectPops insects(surv_j, surv_a, recruit, fecund, K, K_p_mult, s_p, R_,
                       trans_ma, trans_pm, pred_surv, pseudo_surv,
                       extinct_N, demog_error, sigma_x, a, h, k, alate_infl,
                       alate_slope, fly_p, N0, W0, s_y, zeta, Y0);

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
    parasitized.push_back(0);
    mummies.push_back(0);
    wasps.push_back(Y0);

    pcg32 eng;
    seed_pcg(eng);

    for (uint32 t = 0; t < max_t; t++) {
        insects.iterate(eng);
        time.push_back(t+1);
        aphids.push_back(insects.aphids.aphids());
        alates.push_back(insects.aphids.alates());
        parasitized.push_back(insects.aphids.P);
        mummies.push_back(insects.aphids.M);
        wasps.push_back(insects.wasps.Y);
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


