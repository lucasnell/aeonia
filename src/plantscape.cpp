#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng

#ifndef RCPPTHREAD_OVERRIDE_COUT
#define RCPPTHREAD_OVERRIDE_COUT 1    // std::cout override
#endif
#ifndef RCPPTHREAD_OVERRIDE_CERR
#define RCPPTHREAD_OVERRIDE_CERR 1    // std::cerr override
#endif
// #ifndef RCPPTHREAD_OVERRIDE_THREAD
// #define RCPPTHREAD_OVERRIDE_THREAD 1  // std::thread override
// #endif
#include <RcppThread.h>         // multithreading


#include "aeonia-types.hpp"     // integer types
#include "plantscape.hpp"       // PlantScape class
#include "pcg.hpp"              // mt_seeds fxn
#include "util.hpp"             // thread_check fxn



using namespace Rcpp;






//' Create a pointer object in which to store insect population info.
//'
//'
//'
//' @param r Single numeric indicating the growth rate of the aphid population.
//' @param K Single numeric indicating the carrying capacity (no including
//'     predators or *Pseudomonas*) of the aphid population.
//' @param B Single numeric indicating the effect of *Pseudomonas* on aphid
//'     population growth.
//' @param pred_a Single numeric indicating the predator attack rate.
//' @param pred_h Single numeric indicating the predator handling time.
//' @param pred_c Single numeric indicating the predator conversion efficiency.
//' @param pred_m Single numeric indicating the predator mortality rate.
//' @param alate_0 Single numeric.
//'     The proportion of winged offspring from apterous aphids is
//'     `inv_logit(alate_0 + alate_1 * z)` where `z` is the total number of
//'     aphids on that plant.
//' @param alate_1 Single numeric affecting how strongly aphid density
//'     influences alate production. See `alate_0` above for the equation.
//' @param fly_p Single numeric indicating the proportion of alates that fly
//'     off plants each day.
//'
//' @export
//'
//[[Rcpp::export]]
SEXP make_insect_ptr(const double& r,
                     const double& K,
                     const double& B,
                     const double& pred_a,
                     const double& pred_h,
                     const double& pred_c,
                     const double& pred_m,
                     const double& alate_0,
                     const double& alate_1,
                     const double& fly_p) {

    if (r <= 0) stop("r <= 0");
    if (K <= 0) stop("K <= 0");
    if (B < 0) stop("B < 0");
    if (pred_a < 0) stop("pred_a < 0");
    if (pred_h < 0) stop("pred_h < 0");
    if (pred_c < 0) stop("pred_c < 0");
    if (pred_m < 0) stop("pred_m < 0");
    if (alate_1 < 0) stop("alate_1 < 0");
    if (fly_p < 0 || fly_p > 1) stop("fly_p < 0 || fly_p > 1");
    // Densities will be set later:
    double A0 = 0;
    double W0 = 0;
    double P0 = 0;

    XPtr<InsectPops> insect_xptr(new InsectPops(r, K, B, pred_a, pred_h, pred_c,
                                                pred_m, alate_0, alate_1, fly_p,
                                                A0, W0, P0), true);

    return insect_xptr;
}



//' Simulation plant landscapes with virus spread.
//'
//' @param landscapes Integer cube with the types of each plant in each
//'     landscape. It's assumed that rows are x the dimension,
//'     columns are the y dimension, and "slices" (i.e., `dim(landscapes)[3]`)
//'     indicate separate landscapes.
//' @param max_sim_t Single integer giving the maximum time the simulations run.
//' @param insect_ptr External pointer to a C++ object with insect population
//'     information, output from function [make_insect_ptr()].
//' @param A0 Numeric matrix indicating the starting aphid (non-winged) population
//'     density for each plant.
//'     To indicate separate densities for each plant, the matrix should
//'     have the same number of rows and columns as `landscapes`.
//'     The matrix can also be 1x1, in which case it's assumed that all plants
//'     start with the same density of aphids.
//' @param W0 Numeric matrix indicating the starting winged aphid population
//'     density for each plant.
//'     To indicate separate densities for each plant, the matrix should
//'     have the same number of rows and columns as `landscapes`.
//'     The matrix can also be 1x1, in which case it's assumed that all plants
//'     start with the same density of winged aphids.
//' @param P0 Numeric matrix indicating the starting predator population
//'     density for each plant.
//'     To indicate separate densities for each plant, the matrix should
//'     have the same number of rows and columns as `landscapes`.
//'     The matrix can also be 1x1, in which case it's assumed that all plants
//'     start with the same density of predators.
//' @param max_fly_t Single integer indicating the maximum number of visits
//'     an alate can have before it's forced to settle. This is mostly
//'     to avoid computationally troublesome situations when alates simply do
//'     not ever settle, like if `w` is very low.
//' @param radius Max distance that alates will travel between plants.
//'     Defaults to `7.336451`, which is based on previous work.
//'     See "Radius" section below for details.
//' @param alpha Effect of virus infection on alate alighting.
//'     Values `> 0` cause alates to be attracted to virus-infected plants,
//'     while values `< 0` cause them to be repelled by virus-infected plants.
//' @param beta Effect of *Pseudomonas* infection on alate alighting.
//'     Values `> 0` cause alates to be attracted to *Pseudomonas*-infected plants,
//'     while values `< 0` cause them to be repelled by *Pseudomonas*-infected plants.
//' @param epsilon Effect of virus infection on alate acceptance.
//'     Values `> 1` cause alates to be more likely to stay and feed
//'     (indefinitely) on virus-infected plants,
//'     while values `< 1` cause them to be less likely to stay and feed on
//'     virus-infected plants.
//'     Values must be `> 0`, and `epsilon * w` must be `< 1`.
//' @param w Probability that an alate accepts a plant, meaning that
//'     it stays to feed on it indefinitely.
//'     Must be `> 0` and `< 1`. Defaults to `0.2`.
//' @param delta_a Single numeric indicating the probability that an
//'     uninoculated alate is loaded with a virus if it interacts with an
//'     inoculated plant.
//' @param delta_p Single numeric indicating the probability that an
//'     uninoculated plant is loaded with a virus if it interacts with an
//'     inoculated alate.
//' @param total_exp_days Single integer indicating the number of days required
//'     for a plant to transition from exposed (inoculated with virus but
//'     not able to pass it on) to infectious (able to infect other plants).
//'     A value of `0` means that an alate inoculating a plant causes that
//'     plant to be infectious the same day.
//' @param out_by_plant Single logical for whether to split output by plant
//'     instead of summing across the entire landscape.
//' @param show_progress Single logical for whether to show progress bar.
//'     Defaults to `FALSE`.
//' @param n_threads Single integer for the number of threads to use.
//'     Ignored if `dim(landscapes)[3] == 1`.
//'     Defaults to `1L`.
//'
//[[Rcpp::export]]
DataFrame sim_plantscape(const arma::ucube& landscapes,
                         const uint32& max_sim_t,
                         SEXP insect_ptr,
                         const arma::mat& A0,
                         const arma::mat& W0,
                         const arma::mat& P0,
                         const uint32& max_fly_t,
                         const double& radius,
                         const double& alpha,
                         const double& beta,
                         const double& epsilon,
                         const double& w,
                         const double& delta_a,
                         const double& delta_p,
                         const uint32& total_exp_days,
                         const bool& out_by_plant,
                         const bool& show_progress = false,
                         uint32 n_threads = 1) {

    uint32 n_x = landscapes.n_rows;
    uint32 n_y = landscapes.n_cols;
    uint32 n_reps = landscapes.n_slices;

    if (n_x < (uint32)2) stop("nrow(landscapes) < 2");
    if (n_y < (uint32)2) stop("ncol(landscapes) < 2");
    if (n_reps < (uint32)1) stop("dim(landscapes)[3] < 1");
    if (arma::any(arma::vectorise(landscapes) > 3))
        stop("landscapes cannot contain values > 3");

    if (A0.n_rows != n_x && A0.n_rows != 1) stop("nrow(A0) must be 1 or nrow(landscapes)");
    if (A0.n_cols != n_y && A0.n_cols != 1) stop("ncol(A0) must be 1 or ncol(landscapes)");
    if (A0.n_cols == 1 && A0.n_rows != 1) stop("if ncol(A0) == 1, then nrow(A0) must be 1");
    if (A0.n_cols != 1 && A0.n_rows == 1) stop("if nrow(A0) == 1, then ncol(A0) must be 1");
    if (arma::size(A0) != arma::size(W0)) stop("A0, W0, and P0 must be same size");
    if (arma::size(A0) != arma::size(P0)) stop("A0, W0, and P0 must be same size");

    if (max_sim_t == 0 || max_sim_t > 1e6) stop("max_sim_t == 0 || max_sim_t > 1e6");
    if (max_fly_t == 0 || max_fly_t > 1e6) stop("max_fly_t == 0 || max_fly_t > 1e6");
    if (radius < 1) stop("radius < 1");
    if (epsilon < 0) stop("epsilon < 0");
    if (w < 0 || w > 1) stop("w < 0 || w > 1");
    if ((w*epsilon) > 1) stop("w*epsilon > 1");
    if (delta_a < 0 || delta_a > 1) stop("delta_a < 0 || delta_a > 1");
    if (delta_p < 0 || delta_p > 1) stop("delta_p < 0 || delta_p > 1");
    if (total_exp_days > 1e6) stop("total_exp_days > 1e6");

    thread_check(n_threads); // Check that # threads isn't too high

    // Base for all insect populations to start with:
    XPtr<InsectPops> insects_xptr(insects_ptr);
    const InsectPops& insects0(*insects_xptr);

    std::vector<std::vector<uint64>> seeds = mt_seeds(n_reps);
    std::vector<PlantScape> plantscapes;
    plantscapes.reserve(n_reps);

    for (uint32 i = 0; i < n_reps; i++) {
        plantscapes.push_back(PlantScape(max_sim_t, out_by_plant, max_fly_t,
                                         landscapes.slice(i), radius, alpha,
                                         beta, epsilon, w, delta_a, delta_p,
                                         total_exp_days, insects0, A0, W0, P0,
                                         seeds[i]));
    }

    RcppThread::ProgressBar prog_bar(n_reps * max_sim_t, 1);

    if (n_threads > 1U && n_reps > 1U) {
        RcppThread::parallelFor(0, n_reps, [&] (uint32 i) {
            plantscapes[i].run(prog_bar, show_progress);
        }, n_threads);
    } else {
        for (uint32 i = 0; i < n_reps; i++) {
            plantscapes[i].run(prog_bar, show_progress);
        }
    }

    // Produce output dataframe:
    DataFrame out_df;
    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<double> virus;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> preds;
    if (out_by_plant) {
        std::vector<uint32> plant_x;
        std::vector<uint32> plant_y;
        uint32 n_rows = n_reps * (max_sim_t + (uint32)1U) * n_x * n_y;

        rep.reserve(n_rows);
        time.reserve(n_rows);
        plant_x.reserve(n_rows);
        plant_y.reserve(n_rows);
        virus.reserve(n_rows);
        aphids.reserve(n_rows);
        alates.reserve(n_rows);
        preds.reserve(n_rows);

        for (uint32 r = 0; r < n_reps; r++) {
            const arma::mat& rep_out(plantscapes[r].output);
            for (uint32 t = 0; t < rep_out.n_slices; t++) {
                for (uint32 i = 0; i < rep_out.n_rows; i++) {
                    rep.push_back(r+1);
                    time.push_back(t+1);
                    plant_x.push_back(rep_out(i,0,t));
                    plant_y.push_back(rep_out(i,1,t));
                    virus.push_back(rep_out(i,2,t));
                    aphids.push_back(rep_out(i,3,t));
                    alates.push_back(rep_out(i,4,t));
                    preds.push_back(rep_out(i,5,t));
                }
            }
        }

        out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                                   _["x"] = plant_x, _["y"] = plant_y,
                                   _["virus"] = virus, _["aphids"] = aphids,
                                   _["alates"] = alates, _["preds"] = preds);

    } else {
        uint32 n_rows = n_reps * (max_sim_t + (uint32)1U);
        rep.reserve(n_rows);
        time.reserve(n_rows);
        virus.reserve(n_rows);
        aphids.reserve(n_rows);
        alates.reserve(n_rows);
        preds.reserve(n_rows);
        for (uint32 r = 0; r < n_reps; r++) {
            const arma::mat& rep_out(plantscapes[r].output);
            for (uint32 t = 0; t < rep_out.n_slices; t++) {
                rep.push_back(r+1);
                time.push_back(t+1);
                virus.push_back(rep_out(0,0,t));
                aphids.push_back(rep_out(0,1,t));
                alates.push_back(rep_out(0,2,t));
                preds.push_back(rep_out(0,3,t));
            }
        }
        out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                                   _["virus"] = virus, _["aphids"] = aphids,
                                   _["alates"] = alates, _["preds"] = preds);
    }


    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out_df;

}




