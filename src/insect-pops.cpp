
/*
 This contains code for aphid and predator population dynamics in a single patch.
 */

#include <RcppArmadillo.h>
#include <vector>


#include "aeonia_types.hpp"     // integer types
#include "insect-pops.hpp"
#include "util.hpp"  // retrieve_dataset function


using namespace Rcpp;


/*
 Fill parameters from `aeonia::pop_info` if they're not provided.
 */
void fill_pop_info(double& surv_j,
                   double& surv_a,
                   double& recruit,
                   double& fecund,
                   double& K,
                   double& alate_0,
                   double& alate_1,
                   double& a,
                   double& h,
                   double& k,
                   double& s) {

    List pop_info = retrieve_dataset<List>("pop_info");

    if (NumericVector::is_na(surv_j)) surv_j = pop_info["surv_j"];
    if (NumericVector::is_na(surv_a)) surv_a = pop_info["surv_a"];
    if (NumericVector::is_na(recruit)) recruit = pop_info["recruit"];
    if (NumericVector::is_na(fecund)) fecund = pop_info["fecund"];
    if (NumericVector::is_na(K)) K = pop_info["K"];
    if (NumericVector::is_na(alate_0)) alate_0 = pop_info["alate_0"];
    if (NumericVector::is_na(alate_1)) alate_1 = pop_info["alate_1"];
    if (NumericVector::is_na(a)) a = pop_info["a"];
    if (NumericVector::is_na(h)) h = pop_info["h"];
    if (NumericVector::is_na(k)) k = pop_info["k"];
    if (NumericVector::is_na(s)) s = pop_info["s"];

    return;
}




//' Create a pointer object in which to store insect population info.
//'
//' This pointer is used as an argument to [sim_plantscape()].
//'
//'
//' @param B Single numeric indicating the effect of *Pseudomonas* on aphid
//'     population growth.
//' @param a Single numeric indicating the natural enemy attack rate.
//' @param h Single numeric indicating the natural enemy handling time.
//' @param k Single numeric indicating the natural enemy aggregation parameter.
//' @param s Single numeric indicating the natural enemy daily survival.
//' @param fly_p Single numeric indicating the proportion of alates that fly
//'     off plants each day.
//' @param wasp_d_p Single numeric indicating the proportion of adult
//'     parasitoids that are added to the dispersal pool each day.
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
//'     natural enemies or *Pseudomonas*) of the aphid population.
//'     Defaults to `NA`, which results in `pop_info$K` being used.
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
//'
//' @export
//'
//[[Rcpp::export]]
SEXP make_insect_ptr(const double& B,
                     const double& fly_p,
                     const double& wasp_d_p = 0,
                     double surv_j = NA_REAL,
                     double surv_a = NA_REAL,
                     double recruit = NA_REAL,
                     double fecund = NA_REAL,
                     double K = NA_REAL,
                     double alate_0 = NA_REAL,
                     double alate_1 = NA_REAL,
                     double a = NA_REAL,
                     double h = NA_REAL,
                     double k = NA_REAL,
                     double s = NA_REAL) {

    fill_pop_info(surv_j, surv_a, recruit, fecund, K, alate_0, alate_1,
                  a, h, k, s);

    if (B < 0 || B > 1) stop("B < 0 || B > 1");
    if (a < 0) stop("a < 0");
    if (h < 0) stop("h < 0");
    if (k < 0) stop("k < 0");
    if (s < 0 || s > 1) stop("s < 0 || s > 1");
    if (fly_p < 0 || fly_p > 1) stop("fly_p < 0 || fly_p > 1");
    if (wasp_d_p < 0 || wasp_d_p > 1) stop("wasp_d_p < 0 || wasp_d_p > 1");
    if (surv_j <= 0 || surv_j > 1) stop("surv_j <= 0 || surv_j > 1");
    if (surv_a <= 0 || surv_a > 1) stop("surv_a <= 0 || surv_a > 1");
    if (recruit <= 0 || recruit > 1) stop("recruit <= 0 || recruit > 1");
    if (fecund <= 0) stop("fecund <= 0");
    if (K <= 0) stop("K <= 0");
    if (alate_1 < 0) stop("alate_1 < 0");
    // Densities will be set later:
    double A0 = 0;
    double W0 = 0;
    double P0 = 0;

    XPtr<InsectPops> insect_xptr(new InsectPops(surv_j, surv_a, recruit, fecund,
                                                K, B, a, h, k, s,
                                                alate_0, alate_1, fly_p, wasp_d_p,
                                                A0, W0, P0), true);

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
                           const double& A0,
                           const double& W0,
                           const double& P0,
                           const double& B,
                           double surv_j = NA_REAL,
                           double surv_a = NA_REAL,
                           double recruit = NA_REAL,
                           double fecund = NA_REAL,
                           double K = NA_REAL,
                           double alate_0 = NA_REAL,
                           double alate_1 = NA_REAL,
                           double a = NA_REAL,
                           double h = NA_REAL,
                           double k = NA_REAL,
                           double s = NA_REAL) {

    fill_pop_info(surv_j, surv_a, recruit, fecund, K, alate_0, alate_1,
                  a, h, k, s);

    if (max_t > (uint32)1e9) stop("max_t > 1e9");
    if (A0 < 0) stop("A0 < 0");
    if (W0 < 0) stop("W0 < 0");
    if (P0 < 0) stop("P0 < 0");
    if (B < 0 || B > 1) stop("B < 0 || B > 1");
    if (a < 0) stop("a < 0");
    if (h < 0) stop("h < 0");
    if (k < 0) stop("k < 0");
    if (s < 0 || s > 1) stop("s < 0 || s > 1");
    if (surv_j <= 0 || surv_j > 1) stop("surv_j <= 0 || surv_j > 1");
    if (surv_a <= 0 || surv_a > 1) stop("surv_a <= 0 || surv_a > 1");
    if (recruit <= 0 || recruit > 1) stop("recruit <= 0 || recruit > 1");
    if (fecund <= 0) stop("fecund <= 0");
    if (K <= 0) stop("K <= 0");
    if (alate_1 < 0) stop("alate_1 < 0");

    double fly_p = 0;
    double wasp_d_p = 0;

    InsectPops insects(surv_j, surv_a, recruit, fecund, K, B, a, h, k, s,
                       alate_0, alate_1, fly_p, wasp_d_p, A0, W0, P0);

    std::vector<uint32> time;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> enemies;
    time.reserve(max_t+1);
    aphids.reserve(max_t+1);
    alates.reserve(max_t+1);
    enemies.reserve(max_t+1);
    time.push_back(0);
    aphids.push_back(A0);
    alates.push_back(W0);
    enemies.push_back(P0);

    for (uint32 t = 0; t < max_t; t++) {
        insects.iterate();
        time.push_back(t+1);
        aphids.push_back(insects.A());
        alates.push_back(insects.W());
        enemies.push_back(insects.P);
    }

    DataFrame out_df = DataFrame::create(_["time"] = time, _["aphids"] = aphids,
                                         _["alates"] = alates, _["enemies"] = enemies);

    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out_df;

}


