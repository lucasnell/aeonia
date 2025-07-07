#ifndef __AEONIA_ALATE_DISPERSAL_H
#define __AEONIA_ALATE_DISPERSAL_H

/*
 This contains code for alate dispersal.
 */

#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "aeonia_types.hpp"     // integer types
#include "convert-dims.hpp"     // XY
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;






/*
 Info to create flight paths for alates.
 */
class AlateFlightInfo {

    uint32 max_t;

    arma::umat landscape;
    arma::imat neigh_dxdy;

    double alpha;
    double beta;
    double epsilon;
    double w;

    pcg32 eng;  // Random number generator

    arma::vec weights;      // sampling weights
    arma::vec cs_probs;     // cumulative sum of sampling probabilities

    bool virus;
    bool pseudo;


    // Iterate one time step. Adjusts `x` and `y` for new coordinates, and
    // returns true if sims should stop bc alate settles.
    bool iterate(uint32& x, uint32& y) {

        set_virus(x, y);

        // Sample for whether alate will stay to feed at this plant:
        double feed_p = w;
        if (virus) feed_p *= epsilon;
        if (runif_01(eng) < feed_p) {
            return true;
        }

        set_pseudo(x, y);

        // Sample for new location (this function also updates x and y):
        sample_location(x, y);

        return false;

    }

    // Is current plant infected with virus?
    void set_virus(const uint32& x, const uint32& y) {
        set_bit_bool(virus, 0U, landscape(x, y));
        return;
    }
    // Does current plant have Pseudomonas?
    void set_pseudo(const uint32& x, const uint32& y) {
        set_bit_bool(pseudo, 1U, landscape(x, y));
        return;
    }

    // Get the kth bit of n and cast to bool:
    inline bool set_bit_bool(bool& b, const uint32& k, const uint32& n) {
        b = (n & ( 1 << k )) >> k;
        return;
    }

    // Sample new location if alate doesn't stay to feed.
    // Assumes that `virus` and `pseudo` fields are already set.
    void sample_location(uint32& x, uint32& y) {

        // First calculate new sampling weights based on current location:
        uint32 xi, yi;
        double wt_sum = 0;
        for (uint32 i = 0; i < weights.n_elem; i++) {
            const int32& dx(neigh_dxdy(i,0));
            const int32& dy(neigh_dxdy(i,1));
            // This check prevents integer overflow:
            if ((dx < 0 && std::abs(dx) > x) || (dy < 0 && std::abs(dy) > y)) {
                weights(i) = 0.0;
            } else {
                xi = x + dx;
                yi = y + dy;
                weights(i) = 0.0;
                if (xi < landscape.n_rows && yi < landscape.n_cols) {
                    if (virus) weights(i) += alpha;
                    if (pseudo) weights(i) += beta;
                    weights(i) = std::exp(weights(i));
                    wt_sum += weights(i);
                }
            }
        }

        // Now make `cs_probs` into a vector that's the cumulative sum of
        // weights / sum(weights). The last value in `cs_probs` should always be 1.
        cs_probs(0) = weights(0) / wt_sum;
        for (uint32 i = 1; i < cs_probs.n_elem; i++) {
            cs_probs(i) = cs_probs(i-1) + weights(i) / wt_sum;
        }

        // Now sample for a new location:
        double u = runif_01(eng);
        uint32 k = 0;
        uint32 n = cs_probs.n_elem;
        while (k < n && cs_probs(k) < u) k++;

        // Use `k` to update `x` and `y`:
        x += neigh_dxdy(k,0);
        y += neigh_dxdy(k,1);

        return;

    }


public:

    AlateFlightInfo(const uint32& max_t_
                    const arma::umat& landscape_,
                    const double& radius,
                    const double& alpha_,
                    const double& beta_,
                    const double& epsilon_,
                    const double& w_)
        : max_t(max_t_),
          landscape(landscape_),
          neigh_dxdy(),
          dim_conv(landscape_.n_rows, landscape_.n_cols),
          alpha(alpha_),
          beta(beta_),
          epsilon(epsilon_),
          w(w_),
          weights(),
          cs_probs() {

        // Fill `neigh_dxdy` based on radius:
        uint32 total_rows = 0;
        int32 fl_radius = std::floor(radius);
        std::vector<arma::imat> dxdy_vec;
        dxdy_vec.reserve(fl_radius * 2U + 1U);
        int32 max_dy;
        arma::imat dxdy_i;
        double radius2 = radius * radius;
        for (int32 dx = -fl_radius; dx <= fl_radius; dx++) {
            max_dy = std::floor(std::sqrt(radius2 - static_cast<double>(dx * dx)));
            dxdy_i.set_size(max_dy * 2U + 1U, 2U);
            uint32 i = 0;
            for (int32 dy = -max_dy; dy <= max_dy; dy++) {
                dxdy_i(i,0) = dx;
                dxdy_i(i,1) = dy;
                i++;
            }
            dxdy_vec.push_back(dxdy_i);
            total_rows += dxdy_i.n_rows;
        }

        neigh_dxdy.set_size(total_rows, 2U);
        uint32 k = 0;
        for (const arma::imat& dxdy : dxdy_vec) {
            for (uint32 j = 0; j < dxdy.n_rows; j++) {
                neigh_dxdy(k,0) = dxdy(j,0);
                neigh_dxdy(k,1) = dxdy(j,1);
                k++;
            }
        }

        weights.set_size(total_rows);
        cs_probs.set_size(total_rows);
    }


    // Have this alate fly and eventually settle.
    // This function also handles reserving sufficient memory for `out` such
    // that it should never be reallocated.
    void fly(std::vector<XY>& out, const uint32& x0, const uint32& y0) {
        uint32 x = x0;
        uint32 y = y0;
        if (out.capacity() < (max_t+1U)) out.reserve(max_t+1U);
        out.clear();
        out.push_back(XY(x,y));
        bool stop = false;
        uint32 t = 0;
        while (t < max_t && !stop) {
            stop = this->iterate(x, y);
            // Add new coordinates to output:
            out.push_back(XY(x, y));
            t++;
        }
        return;
    }



};


//
// void check_a_s_args(const uint32& max_t,
//                     const arma::imat& plant_xy,
//                     const arma::ivec& plant_types,
//                     const double& alpha,
//                     const double& beta,
//                     const double& epsilon,
//                     const double& w,
//                     const double& radius,
//                     Nullable<IntegerVector> xy0,
//                     const bool& randomize_xy0,
//                     arma::imat& landscape,
//                     uint32& x0,
//                     uint32& y0,
//                     const uint32& n_alates,
//                     const std::string& summarize,
//                     uint32& n_threads) {
//
//     // Check that # threads isn't too high:
//     thread_check(n_threads);
//
//     if (max_t == 0) stop("max_t == 0");
//     if (max_t > 1000000000U) stop("max_t > 1e9");
//
//     if (plant_xy.n_cols != 2) stop("ncol(plant_xy) != 2");
//     if (plant_xy.n_rows < 2) stop("nrow(plant_xy) < 2");
//     if (plant_xy.min() < 1) stop("min(plant_xy) < 1");
//     if (plant_xy.max() > 1000000000U) stop("max(plant_xy) > 1e9");
//
//     if (plant_types.min() < 1) stop("min(plant_types) < 1");
//     if (plant_types.max() > 4) stop("max(plant_types) > 4");
//     if (plant_types.n_elem != plant_xy.n_rows)
//         stop("length(plant_types) != nrow(plant_xy)");
//
//     if (w < 0 || w > 1) stop("w < 0 || w > 1");
//     if (epsilon < 0) stop("epsilon < 0");
//     if ((w*epsilon) > 1) stop("w*epsilon > 1");
//
//     if (radius < 1) stop("radius < 1");
//
//     if (n_alates == 0) stop("n_alates == 0");
//     if (n_alates > 1000000000U) stop("n_alates > 1e9");
//
//     if (summarize != "none" && summarize != "types" && summarize != "xy") {
//         stop("summarize must be \"none\", \"types\", or \"xy\".");
//     }
//
//     // Fill landscape and check that all locations are filled:
//     int32 n_x = plant_xy.col(0).max();
//     int32 n_y = plant_xy.col(1).max();
//     // Quick check for whether number of rows is consistent with all locations
//     // being referenced:
//     if (plant_xy.n_rows != (n_x * n_y))
//         stop("nrow(plant_xy) != (max(plant_xy[,1]) * max(plant_xy[,2]))");
//     // Start with all 100's so I can later check if any haven't been filled:
//     landscape.set_size(n_x, n_y);
//     landscape.fill(100);
//     for (uint32 i = 0; i < plant_xy.n_rows; i++) {
//         const int32& x(plant_xy(i,0));
//         const int32& y(plant_xy(i,1));
//         landscape(x-1U, y-1U) = plant_types(i);
//     }
//     if (landscape.max() > 4) {
//         stop("plant_xy and plant_types must provide a type for all locations");
//     }
//
//     // Fill starting x and y if provided:
//     if (xy0.isNotNull()) {
//         if (randomize_xy0) {
//             std::string msg("If `randomize_xy0 = TRUE`, this will negate ");
//             msg += "your provided starting values for `xy0`";
//             stop(msg.c_str());
//         }
//         IntegerVector xy0_vec(xy0);
//         if (xy0_vec.size() != 2) {
//             stop("xy0_vec must be NULL or a length-2 numeric vector");
//         }
//         if (xy0_vec(0) < 1 || xy0_vec(0) > n_x) {
//             stop("xy0[1] must be 1-based index for a value in plant_xy[,1]");
//         }
//         if (xy0_vec(1) < 1 || xy0_vec(1) > n_y) {
//             stop("xy0[2] must be 1-based index for a value in plant_xy[,2]");
//         }
//         x0 = xy0_vec(0) - 1;
//         y0 = xy0_vec(1) - 1;
//     }
//
//     return;
//
// }
//
//
//
//
//
// /*
//  Fill the following objects:
//  - binary variable for whether a plant has virus (`v`)
//  - binary variable for whether a plant has *Pseudomonas* (`b`)
//  - change in x and y for neighboring plants (those within radius) (`neigh_dxdy`)
//  */
// void fill_v_b_neigh(const arma::imat& landscape,
//                     const double& radius,
//                     arma::mat& v,
//                     arma::mat& b,
//                     arma::imat& neigh_dxdy) {
//
//     // Fill v and b based on landscape:
//     v.set_size(arma::size(landscape));
//     b.set_size(arma::size(landscape));
//     for (uint32 i = 0; i < landscape.n_rows; i++) {
//         for (uint32 j = 0; j < landscape.n_cols; j++) {
//             const int32& type(landscape(i,j));
//             if (type == 1 || type == 4) {
//                 v(i,j) = 1;
//             } else v(i,j) = 0;
//             if (type == 2 || type == 4) {
//                 b(i,j) = 1;
//             } else b(i,j) = 0;
//         }
//     }
//
//     // Fill `neigh_dxdy` based on radius:
//     uint32 total_rows = 0;
//     int32 fl_radius = std::floor(radius);
//     std::vector<arma::imat> dxdy_vec;
//     dxdy_vec.reserve(fl_radius * 2U + 1U);
//     int32 max_dy;
//     arma::imat dxdy_i;
//     double radius2 = radius * radius;
//     for (int32 dx = -fl_radius; dx <= fl_radius; dx++) {
//         max_dy = std::floor(std::sqrt(radius2 - static_cast<double>(dx * dx)));
//         dxdy_i.set_size(max_dy * 2U + 1U, 2U);
//         uint32 i = 0;
//         for (int32 dy = -max_dy; dy <= max_dy; dy++) {
//             dxdy_i(i,0) = dx;
//             dxdy_i(i,1) = dy;
//             i++;
//         }
//         dxdy_vec.push_back(dxdy_i);
//         total_rows += dxdy_i.n_rows;
//     }
//
//     neigh_dxdy.set_size(total_rows, 2U);
//     uint32 k = 0;
//     for (const arma::imat& dxdy : dxdy_vec) {
//         for (uint32 j = 0; j < dxdy.n_rows; j++) {
//             neigh_dxdy(k,0) = dxdy(j,0);
//             neigh_dxdy(k,1) = dxdy(j,1);
//             k++;
//         }
//     }
//
//     return;
//
// }
//
//
//
//
// DataFrame create_output(const std::vector<OneAlate>& alates,
//                         const std::string& summarize) {
//
//     DataFrame out;
//
//
//     if (summarize == "none") {
//
//         uint32 n_rows = 0;
//         for (const OneAlate& ala : alates) n_rows += ala.size();
//
//         std::vector<int> alate;
//         std::vector<int> time;
//         std::vector<int> x;
//         std::vector<int> y;
//         std::vector<int> to;
//         std::vector<int> from;
//
//         alate.reserve(n_rows);
//         time.reserve(n_rows);
//         x.reserve(n_rows);
//         y.reserve(n_rows);
//         to.reserve(n_rows);
//         from.reserve(n_rows);
//
//         for (uint32 i = 0; i < alates.size(); i++) {
//             for (uint32 t = 0; t < alates[i].size(); t++) {
//                 alate.push_back(i+1);
//                 time.push_back(t);
//                 x.push_back(alates[i].x[t]+1);
//                 y.push_back(alates[i].y[t]+1);
//                 to.push_back(alates[i].to[t]);
//                 from.push_back(alates[i].from[t]);
//             }
//         }
//
//         out = DataFrame::create(
//             _["alate"] = alate,
//             _["time"] = time,
//             _["x"] = x,
//             _["y"] = y,
//             _["to"] = to,
//             _["from"] = from);
//
//     } else if (summarize == "types") {
//
//         std::vector<int> to;
//         std::vector<int> from;
//         std::vector<int> n;
//         uint32 n_rows = 16; // 16 is (number of types)^2
//
//         to.reserve(n_rows);
//         from.reserve(n_rows);
//         n.reserve(n_rows);
//
//         for (uint32 t = 1; t <= 4; t++) {
//             for (uint32 f = 1; f <= 4; f++) {
//                 to.push_back(t);
//                 from.push_back(f);
//                 n.push_back(0);
//             }
//         }
//
//         for (uint32 i = 0; i < alates.size(); i++) {
//             if (alates[i].size() < 2) continue;
//             // Now count numbers of `from` --> `to` hits and store in `n`:
//             uint32 to_t, from_t, k_t;
//             // (below, starting at 2nd step bc 1st has no reasonable `from`)
//             for (uint32 t = 1; t < alates[i].size(); t++) {
//                 if (alates[i].to[t] < 1) stop("alates[i].to[t] < 1");
//                 if (alates[i].from[t] < 1) stop("alates[i].from[t] < 1");
//                 to_t = alates[i].to[t] - 1;
//                 from_t = alates[i].from[t] - 1;
//                 k_t = to_t * 4U + from_t;
//                 if (k_t >= n_rows) stop("k_t >= 16");
//                 n[k_t] += 1;
//             }
//         }
//
//         out = DataFrame::create(
//             _["to"] = to,
//             _["from"] = from,
//             _["n"] = n);
//
//     } else {
//
//         uint32 n_rows = alates.size();
//
//         std::vector<int> alate;
//         std::vector<int> x;
//         std::vector<int> y;
//
//         alate.reserve(n_rows);
//         x.reserve(n_rows);
//         y.reserve(n_rows);
//
//         for (uint32 i = 0; i < alates.size(); i++) {
//             alate.push_back(i+1);
//             x.push_back(alates[i].x.back()+1);
//             y.push_back(alates[i].y.back()+1);
//         }
//
//         out = DataFrame::create(
//             _["alate"] = alate,
//             _["x"] = x,
//             _["y"] = y);
//
//     }
//
//     out.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});
//
//     return out;
//
// }
//
//
//
//
//
//
//
//
// //' Simulate aphid alates seeking host plants.
// //'
// //'
// //' @details # Radius
// //' From "The Role of Aphid Behaviour in the Epidemiology of Potato Virus Y:
// //' a Simulation Study" by Thomas Nemecek (1993; p. 72), dispersal distances
// //' follow a Weibull distribution with shape = 0.6569 and scale = 9.613.
// //'
// //' The default for the `radius` argument uses the median of this
// //' distribution.
// //' I'm dividing by 0.75 to convert from meters to plant locations that are
// //' 0.75 meters apart (typical spacing for pea):
// //' `radius = qweibull(0.5, 0.6569, 9.613) / 0.75`.
// //'
// //'
// //' @param max_t Single integer indicating the number of time steps to simulate.
// //' @param plant_xy Two-column integer matrix containing x and y coordinates
// //'     for plants. Note that all locations from 1 to the max in each dimension
// //'     must be represented by a plant.
// //' @param plant_types Integer vector with the same number of items as the
// //'     number of rows in `plant_xy` indicating the type of plant each
// //'     plant is. Values in this vector must be 1 (virus-infected),
// //'     2 (*Pseudomonas*-infected), 3 (none), or 4 (both).
// //'     Not all values must be present in the landscape.
// //'     Note that all locations from 1 to the max in each dimension
// //'     must be represented by a plant.
// //' @param alpha Effect of virus infection on alate alighting.
// //'     Values `> 0` cause alates to be attracted to virus-infected plants,
// //'     while values `< 0` cause them to be repelled by virus-infected plants.
// //' @param beta Effect of *Pseudomonas* infection on alate alighting.
// //'     Values `> 0` cause alates to be attracted to *Pseudomonas*-infected plants,
// //'     while values `< 0` cause them to be repelled by *Pseudomonas*-infected plants.
// //' @param epsilon Effect of virus infection on alate acceptance.
// //'     Values `> 1` cause alates to be more likely to stay and feed
// //'     (indefinitely) on virus-infected plants,
// //'     while values `< 1` cause them to be less likely to stay and feed on
// //'     virus-infected plants.
// //'     Values must be `> 0`, and `epsilon * w` must be `< 1`.
// //' @param w Probability that an alate accepts a plant, meaning that
// //'     it stays to feed on it indefinitely.
// //'     Must be `> 0` and `< 1`. Defaults to `0.2`.
// //' @param radius Max distance that alates will travel between plants.
// //'     Defaults to `7.336451`, which is based on previous work.
// //'     See "Radius" section below for details.
// //' @param xy0 Numeric vector of length 2 indicating starting x and y
// //'     coordinates for alate(s). When `xy0 = NULL` and
// //'     `randomize_xy0 = FALSE`, alate(s) start at the middle of the
// //'     search space (x = 0, y = 0).
// //'     This function throws an error when `xy0` is provided and
// //'     `randomize_xy0 = TRUE` because these are conflicting.
// //'     Defaults to `NULL`.
// //' @param randomize_xy0 Single logical for whether to randomize starting
// //'     coordinates for alate(s). If `TRUE`, x and y coordinates are
// //'     randomly chosen from the landscape.
// //'     If `FALSE`, alates start in the middle of the landscape.
// //'     Defaults to `TRUE`.
// //' @param n_alates Single integer indicating the number of independent
// //'     alates to simulate.
// //'     Defaults to `1L`.
// //' @param summarize Single string for whether and how to summarize output.
// //'     Options are `"none"`, `"types"`, or `"xy"`.
// //'     See below for details on how this changes the output.
// //'     Defaults to `"none"`.
// //' @param show_progress Single logical for whether to show progress bar.
// //'     Defaults to `FALSE`.
// //' @param n_threads Single integer for the number of threads to use.
// //'     Ignored if `n_alates == 1`.
// //'     Defaults to `1L`.
// //'
// //' @returns If `summarize = "none"`, a tibble with the columns
// //'     `alate` (alate number),
// //'     `time` (time), `x` (x coordinate), `y` (y coordinate),
// //'    `to` (the plant type the alate traveled to), and
// //'     `from` (the plant type the alate traveled from).
// //'     If `summarize = "types"`, then it outputs a tibble with the columns
// //'     `to`, `from`, and
// //'     `n` (the number of times all searchers traveled from this plant type to
// //'     this other type).
// //'     If `summarize = "xy"`, then it outputs a tibble with the columns
// //'     `alate`, `x` (the final x coordinate), and
// //'     `y` (the final y coordinate).
// //'
// //'
// //' @export
// //'
// //[[Rcpp::export]]
// DataFrame alate_search_sims(const uint32& max_t,
//                             const arma::imat& plant_xy,
//                             const arma::ivec& plant_types,
//                             const double& alpha,
//                             const double& beta,
//                             const double& epsilon = 1,
//                             const double& w = 0.2,
//                             const double& radius = 7.336451,
//                             Nullable<IntegerVector> xy0 = R_NilValue,
//                             const bool& randomize_xy0 = true,
//                             const uint32& n_alates = 1,
//                             const std::string& summarize = "none",
//                             const bool& show_progress = false,
//                             uint32 n_threads = 1) {
//
//     uint32 x0 = plant_xy.col(0).max() / 2U;
//     uint32 y0 = plant_xy.col(1).max() / 2U;
//     // Landscape where x is rows, y is cols:
//     arma::imat landscape;
//
//     // Check arguments, fill `landscape`, and optionally set x0 and y0:
//     check_a_s_args(max_t, plant_xy, plant_types, alpha, beta, epsilon, w, radius,
//                    xy0, randomize_xy0, landscape, x0, y0, n_alates, summarize, n_threads);
//
//     arma::mat v;  // Binary variable for whether a plant has virus
//     arma::mat b;  // Binary variable for whether a plant has *Pseudomonas*
//     arma::imat neigh_dxdy; // change in x and y for neighbors (those within radius)
//     // Fill these objects:
//     fill_v_b_neigh(landscape, radius, v, b, neigh_dxdy);
//
//     // Create a OneAlate object per alate.
//     // This also generates a seeded RNG for each, using R's runif(...) for
//     // reproducibility.
//     std::vector<OneAlate> alates;
//     alates.reserve(n_alates);
//     for (uint32 i = 0; i < n_alates; i++) {
//         alates.push_back(OneAlate(max_t, x0, y0, landscape, v, b,
//                                   neigh_dxdy, alpha, beta, epsilon, w,
//                                   randomize_xy0));
//     }
//
//
//     RcppThread::ProgressBar prog_bar(n_alates, 1);
//
//     if (n_threads > 1U && n_alates > 1U) {
//         RcppThread::parallelFor(0, n_alates, [&] (uint32 i) {
//             alates[i].run(prog_bar, show_progress);
//         }, n_threads);
//     } else {
//         for (uint32 i = 0; i < n_alates; i++) {
//             alates[i].run(prog_bar, show_progress);
//         }
//     }
//
//
//     DataFrame out = create_output(alates, summarize);
//
//     return out;
// }



#endif
