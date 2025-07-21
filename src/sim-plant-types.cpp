
#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <pcg/pcg_random.hpp>   // pcg prng


#include "aeonia_types.hpp"     // integer types
#include "pcg.hpp"              // runif_01, seed_rng functions
#include "sim-plant-types.hpp"  // DimensionConverter and LocationSampler classes
#include "util.hpp"             // thread_check



using namespace Rcpp;



//' Simulate field(s) of plant types.
//'
//' Simulate locations of plants of different types
//' (virus-infected or -uninfected, *Pseudomonas*-infected or uninfected)
//' along an evenly spaced grid of integers, where the placement of one
//' target can affect subsequent placement of other targets.
//' Target locations are drawn from all combinations of `1` to
//' `n_x` and `1` to `n_y`.
//'
//'
//' @param n_x Single integer indicating x dimension of search area.
//'     Locations will be drawn from `1` to `n_x`.
//'     Must be at least `2` and less than `1e6`.
//' @param n_y Single integer indicating y dimension of search area.
//'     Locations will be drawn from `1` to `n_y`.
//'     Must be at least `2` and less than `1e6`.
//' @param wt_mat Square numeric matrix indicating how sample weighting on
//'     neighboring locations is affected by a target of each type being
//'     placed in a particular spot.
//'     Item `wt_mat[i,j]` indicates the effect of target type `i` on
//'     subsequent samplings of target type `j`.
//'     Values above 1 cause neighboring locations to be more likely to be
//'     sampled later, while values below 1 cause them to be less likely
//'     sampled.
//'     For locations that have been adjusted using `wt_mat` multiple times,
//'     the weights are multiplied by each other
//'     (e.g., `w *= wt_mat[1,3]` at time `t`, then
//'     `w *= wt_mat[1,2]` at time `t+1`).
//'     All weights start with values of `1`.
//'     The matrix should have the same number of rows and columns as the
//'     number of target types.
//' @param n_virus Single integer indicating the number of plants that should
//'     be virus infected. This value must be `>=1` and `<= n_x * n_y`.
//' @param n_pseudo Single integer indicating the number of plants that should
//'     be *Pseudomonas* infected.
//'     This value must be `>=1` and `<= n_x * n_y`.
//' @param n_lands Single integer indicating the number of independent
//'     landscapes to simulate. Each landscape will have a number of samples
//'     per type according to `n_virus` and `n_pseudo`.
//'     Must be `> 0` and `< 1e6`.
//'     Defaults to `1`.
//' @param show_progress Single logical for whether to show progress bar.
//'     Defaults to `FALSE`.
//' @param n_threads Single integer for the number of threads to use.
//'     Multithreading happens across different landscapes, so this is
//'     ignored if `n_lands == 1`.
//'     Defaults to `1L`.
//'
//' @return A cube where rows indicate the x coordinate, columns indicate the
//'     y coordinate, and slices indicate the landscape number.
//'     For each item,
//'     `0` indicates nothing,
//'     `1` indicates a virus infected plant (no *Pseudomonas*),
//'     `2` indicates a *Pseudomonas*-infected plant (no virus), and
//'     `3` indicates a plant with both virus and *Pseudomonas*.
//'     See [land_cube2list()] for converting this to a list of length `n_lands`,
//'     where each item is a [`tibble`][tibble::tbl_df] with columns
//'     `x`, `y`, and `type`.
//'
//'
//' @export
//'
//[[Rcpp::export]]
arma::ucube sim_plant_types(const uint32& n_x,
                            const uint32& n_y,
                            const arma::mat& wt_mat,
                            const uint32& n_virus,
                            const uint32& n_pseudo,
                            const uint32& n_lands = 1,
                            const bool& show_progress = false,
                            uint32 n_threads = 1) {

    if (n_x < 1) stop("n_x must be >= 1");
    if (n_x > 1e6) stop("n_x must be <= 1e6");
    if (n_y < 1) stop("n_y must be >= 1");
    if (n_y > 1e6) stop("n_y must be <= 1e6");
    if (n_x * n_y < 2) stop("n_x * n_y must be >= 2");
    if (wt_mat.n_rows != 2 || wt_mat.n_cols != 2) stop("wt_mat must be 2x2");
    if (arma::any(arma::vectorise(wt_mat) < 0)) stop("wt_mat cannot contain values < 0");
    if ((n_virus + n_pseudo) < 1) {
        std::string err_msg = "n_virus + n_pseudo cannot be < 1; ";
        err_msg += "just use an array of zeros instead of using this function.";
        stop(err_msg.c_str());
    }
    if (n_virus > n_x * n_y) stop("n_virus cannot be > n_x * n_y");
    if (n_pseudo > n_x * n_y) stop("n_pseudo cannot be > n_x * n_y");
    if (n_lands < 1) stop("n_lands must be >= 1");
    if (n_lands >= 1e6) stop("n_lands must be < 1e6");
    thread_check(n_threads); // Check that # threads isn't too high

    RcppThread::ProgressBar prog_bar(n_lands * (n_virus + n_pseudo), 1);

    std::vector<OnePlantTypeSimmer> simmers;
    simmers.reserve(n_lands);
    for (uint32 i = 0; i < n_lands; i++) {
        simmers.push_back(OnePlantTypeSimmer(wt_mat, n_virus, n_pseudo, n_x, n_y));
    }


    if (n_threads > 1U && n_lands > 1U) {
        auto job = [&] (OnePlantTypeSimmer& simmer) {
            simmer.run(prog_bar, show_progress);
        };
        RcppThread::parallelForEach(simmers, job, n_threads);
    } else {
        for (uint32 i = 0; i < n_lands; i++) {
            simmers[i].run(prog_bar, show_progress);
        }
    }

    arma::ucube out(n_x, n_y, n_lands);
    for (uint32 i = 0; i < n_lands; i++) {
        simmers[i].fill_output(out, i);
    }

    return out;

}



//' Convert a landscape cube to a list of dataframes.
//'
//' @param land_cube Cube containing landscape info, output from [sim_plant_types()].
//'
//' @return A list of length `dim(land_cube)[3]`, where each item is a
//'     [`tibble`][tibble::tbl_df] with columns
//'     `x`, `y`, and `type`.
//'
//[[Rcpp::export]]
List land_cube2list(const arma::ucube& land_cube) {

    uint32 n_x = land_cube.n_rows;
    uint32 n_y = land_cube.n_cols;
    uint32 n_lands = land_cube.n_slices;

    if (n_x < 2) stop("nrow(land_cube) must be >= 2");
    if (n_x > 1e6) stop("nrow(land_cube) must be <= 1e6");
    if (n_y < 2) stop("ncol(land_cube) must be >= 2");
    if (n_y > 1e6) stop("ncol(land_cube) must be <= 1e6");
    if (n_lands < 1) stop("dim(land_cube)[3] must be >= 1");
    if (n_lands >= 1e6) stop("dim(land_cube)[3] must be < 1e6");

    std::vector<uint32> x;
    std::vector<uint32> y;
    std::vector<uint32> type;
    x.reserve(n_x * n_y);
    y.reserve(n_x * n_y);
    type.reserve(n_x * n_y);

    List out(n_lands);

    DataFrame out_k;

    // Do first landscape while filling x and y (and setting size for type):
    for (uint32 j = 0; j < n_y; j++) {
        for (uint32 i = 0; i < n_x; i++) {
            x.push_back(i+1);
            y.push_back(j+1);
            type.push_back(land_cube(i,j,0));
        }
    }
    out_k = DataFrame::create(_["x"] = x, _["y"] = y, _["type"] = type);
    out_k.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});
    out[0] = out_k;

    if (n_lands > 1) {
        for (uint32 k = 1; k < n_lands; k++) {
            uint32 t = 0;
            for (uint32 j = 0; j < n_y; j++) {
                for (uint32 i = 0; i < n_x; i++) {
                    type[t] = land_cube(i,j,k);
                    t++;
                }
            }
            out_k = DataFrame::create(_["x"] = x, _["y"] = y,
                                       _["type"] = type);
            out_k.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});
            out[k] = out_k;
        }
    }

    return out;

}
