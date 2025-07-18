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


#include "aeonia_types.hpp"     // integer types
#include "plantscape.hpp"       // PlantScape class
#include "pcg.hpp"              // mt_seeds fxn
#include "util.hpp"             // thread_check fxn



using namespace Rcpp;








//' Simulation plant landscapes with virus spread.
//'
//' @details # Radius
//' From "The Role of Aphid Behaviour in the Epidemiology of Potato Virus Y:
//' a Simulation Study" by Thomas Nemecek (1993; p. 72), dispersal distances
//' follow a Weibull distribution with shape = 0.6569 and scale = 9.613.
//'
//' The default for the `radius` argument uses the median of this
//' distribution.
//' I'm dividing by 0.75 to convert from meters to plant locations that are
//' 0.75 meters apart (typical spacing for pea):
//' `radius = qweibull(0.5, 0.6569, 9.613) / 0.75`.
//'
//' @param landscapes Integer cube with the types of each plant in each
//'     landscape. It's assumed that rows are x the dimension,
//'     columns are the y dimension, and "slices" (i.e., `dim(landscapes)[3]`)
//'     indicate separate landscapes.
//' @param max_t Single integer giving the maximum time the simulations run.
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
//'     Defaults to `7`, which is based on the paper "Cucumber mosaic virus
//'     isolates seedborne in *Phaseolus vulgaris*: serology, host-pathogen
//'     relationships, and seed transmission" (Davis & Hampton, 1986).
//' @param w Probability that an alate accepts a plant, meaning that
//'     it stays to feed on it indefinitely.
//'     Must be `> 1e-4` and `< 1`.
//'     The lower limit is because a very small value of `w` causes alates to
//'     fly to so many plants that it becomes computationally problematic.
//'     This is the same reason that I check for (and return and error) if
//'     `w*epsilon < 1e-4`.
//'     Defaults to `0.2`.
//' @param radius Max distance that alates will travel between plants.
//'     Defaults to `7.336451`, which is based on previous work.
//'     See "Radius" section below for details.
//' @param out_by_plant Single logical for whether to split output by plant
//'     instead of summing across the entire landscape.
//'     Defaults to `TRUE`.
//' @param infect_stop Single logical for whether to stop simulations
//'     when all plants are infected with virus.
//'     Defaults to `TRUE`.
//' @param show_progress Single logical for whether to show progress bar.
//'     Defaults to `FALSE`.
//' @param n_threads Single integer for the number of threads to use.
//'     Ignored if `dim(landscapes)[3] == 1`.
//'     Defaults to `1L`.
//'
//'
//' @export
//'
//[[Rcpp::export]]
DataFrame sim_plantscape(const arma::ucube& landscapes,
                         const uint32& max_t,
                         SEXP insect_ptr,
                         const arma::mat& A0,
                         const arma::mat& W0,
                         const arma::mat& P0,
                         const double& alpha,
                         const double& beta,
                         const double& epsilon,
                         const double& delta_a,
                         const double& delta_p,
                         const uint32& total_exp_days = 7,
                         const double& w = 0.2,
                         const double& radius = 7.336451,
                         const bool& out_by_plant = true,
                         const bool& infect_stop = true,
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

    if (max_t == 0 || max_t > 1e6) stop("max_t == 0 || max_t > 1e6");
    if (radius < 1) stop("radius < 1");
    if (epsilon < 0) stop("epsilon < 0");
    if (w < 0.0001 || w > 1) stop("w < 0.0001 || w > 1");
    if ((w*epsilon) > 1) stop("w*epsilon > 1");
    if ((w*epsilon) < 0.0001) stop("w*epsilon < 0.0001");
    if (delta_a < 0 || delta_a > 1) stop("delta_a < 0 || delta_a > 1");
    if (delta_p < 0 || delta_p > 1) stop("delta_p < 0 || delta_p > 1");
    if (total_exp_days > 1e6) stop("total_exp_days > 1e6");

    thread_check(n_threads); // Check that # threads isn't too high

    // Set a maximum on # plants an alate can fly to such that the probability
    // that it reaches this threshold < 1e-9:
    uint32 max_fly_t = 100;
    double max_leave_p = 1 - std::min(w, epsilon * w); // max Pr(leave plant)
    while (std::pow(max_leave_p, max_fly_t) > (double)1e-9) {
        max_fly_t *= 10;
        if (max_fly_t > (uint32)1000000) stop("INTERNAL ERROR: max_fly_t too high");
    }

    // Make sure output object won't be too big for R:
    uint32 n_rows = n_reps * (max_t + (uint32)1U);
    if (out_by_plant) n_rows *= (n_x * n_y);
    if (n_rows > (uint32)2147483647)
        stop("This combo of parameters will produce too large of an output for R");

    // Base for all insect populations to start with:
    XPtr<InsectPops> insect_xptr(insect_ptr);
    const InsectPops& insects0(*insect_xptr);

    std::vector<std::vector<uint64>> seeds = mt_seeds(n_reps);
    std::vector<PlantScape> plantscapes;
    plantscapes.reserve(n_reps);

    for (uint32 i = 0; i < n_reps; i++) {
        plantscapes.push_back(PlantScape(max_t, out_by_plant, max_fly_t,
                                         landscapes.slice(i), radius, alpha,
                                         beta, epsilon, w, delta_a, delta_p,
                                         total_exp_days, insects0, A0, W0, P0,
                                         seeds[i]));
    }

    RcppThread::ProgressBar prog_bar(n_reps * max_t, 1);

    if (n_threads > 1U && n_reps > 1U) {
        RcppThread::parallelFor(0, n_reps, [&] (uint32 i) {
            plantscapes[i].run(infect_stop, prog_bar, show_progress);
        }, n_threads);
    } else {
        for (uint32 i = 0; i < n_reps; i++) {
            plantscapes[i].run(infect_stop, prog_bar, show_progress);
        }
    }

    // Produce output dataframe:
    DataFrame out_df;
    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<double> virus;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> enemies;

    // // If sims stop when all plants are infected, update `n_rows`:
    // if (infect_stop) {
    //     n_rows = 0;
    //     for (uint32 r = 0; r < n_reps; r++) {
    //         for (uint32 t = 0; t < plantscapes[r].output.size(); t++) {
    //             n_rows += plantscapes[r].output[t].n_rows;
    //         }
    //     }
    // }

    if (out_by_plant) {

        std::vector<uint32> plant_x;
        std::vector<uint32> plant_y;

        rep.reserve(n_rows);
        time.reserve(n_rows);
        plant_x.reserve(n_rows);
        plant_y.reserve(n_rows);
        virus.reserve(n_rows);
        aphids.reserve(n_rows);
        alates.reserve(n_rows);
        enemies.reserve(n_rows);

        for (uint32 r = 0; r < n_reps; r++) {
            for (uint32 t = 0; t < plantscapes[r].output.size(); t++) {
                const arma::mat& rep_out(plantscapes[r].output[t]);
                for (uint32 i = 0; i < rep_out.n_rows; i++) {
                    rep.push_back(r+1);
                    time.push_back(t+1);
                    plant_x.push_back(rep_out(i,0));
                    plant_y.push_back(rep_out(i,1));
                    virus.push_back(rep_out(i,2));
                    aphids.push_back(rep_out(i,3));
                    alates.push_back(rep_out(i,4));
                    enemies.push_back(rep_out(i,5));
                }
            }
        }

        out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                                   _["x"] = plant_x, _["y"] = plant_y,
                                   _["virus"] = virus, _["aphids"] = aphids,
                                   _["alates"] = alates, _["enemies"] = enemies);

    } else {

        rep.reserve(n_rows);
        time.reserve(n_rows);
        virus.reserve(n_rows);
        aphids.reserve(n_rows);
        alates.reserve(n_rows);
        enemies.reserve(n_rows);

        for (uint32 r = 0; r < n_reps; r++) {
            for (uint32 t = 0; t < plantscapes[r].output.size(); t++) {
                const arma::mat& rep_out(plantscapes[r].output[t]);
                rep.push_back(r+1);
                time.push_back(t+1);
                virus.push_back(rep_out(0,0));
                aphids.push_back(rep_out(0,1));
                alates.push_back(rep_out(0,2));
                enemies.push_back(rep_out(0,3));
            }
        }

        out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                                   _["virus"] = virus, _["aphids"] = aphids,
                                   _["alates"] = alates, _["enemies"] = enemies);
    }


    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out_df;

}




