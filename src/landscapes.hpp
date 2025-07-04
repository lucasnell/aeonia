#ifndef __AEONIA_LANDSCAPES_H
#define __AEONIA_LANDSCAPES_H

#include <RcppArmadillo.h>
#include <vector>
#include <deque>
#include <string>
#include <limits>
#include <numeric>
#include <pcg/pcg_random.hpp>   // pcg prng


#include "aeonia_types.hpp"
#include "pcg.hpp"              // runif_ab fxn


using namespace Rcpp;





/*
 Class to convert back and forth between 2D to 1D indices.
 It assumes that x and y coordinates are sorted first by y coordinates, then
 by x coordinates.
 It also assumes 0-based indices.
 An example matrix of x and y coordinates might start like this:

 #>      x y
 #> [0,] 0 0
 #> [1,] 1 0
 #> [2,] 2 0
 #> [3,] 3 0
 #> [4,] 4 0
 #> [5,] 5 0

 */
class DimensionConverter {

    std::vector<uint32> neigh_x;
    std::vector<uint32> neigh_y;

    uint32 x_size;
    uint32 y_size;

public:

    DimensionConverter(const uint32& x_size_, const uint32& y_size_)
        : neigh_x(), neigh_y(), x_size(x_size_), y_size(y_size_) {
        neigh_x.reserve(3);
        neigh_y.reserve(3);
    }

    // Convert from 1D to 2D:
    void to_2d(uint32& x, uint32& y, const uint32& k) const {
        x = k - y_size * (k / y_size);
        y = k / y_size;
        return;
    }
    // Overloaded for signed ints (for use with Rcpp::IntegerVector)
    void to_2d(int& x, int& y, const uint32& k) const {
        x = k - y_size * (k / y_size);
        y = k / y_size;
        return;
    }
    // Convert from 2D to 1D:
    void to_1d(uint32& k, const uint32& x, const uint32& y) const {
        k = (y * x_size + x);
        return;
    }
    /*
     Return indices (in 1D) for all neighbors based on a 1D input coordinate.
     NOTE:
        - It clears `indices` before adding to it
        - It also returns an index for the focal point
     */
    void get_neighbors(std::vector<uint32>& indices,
                       const uint32& k) {
        uint32 x0 = k - y_size * (k / y_size);
        uint32 y0 = k / y_size;
        indices.clear();
        neigh_x.clear();
        neigh_y.clear();
        if (x0 > 0) neigh_x.push_back(x0-1);
        neigh_x.push_back(x0);
        if (x0 < x_size-1) neigh_x.push_back(x0+1);
        if (y0 > 0) neigh_y.push_back(y0-1);
        neigh_y.push_back(y0);
        if (y0 < x_size-1) neigh_y.push_back(y0+1);
        uint32 k_out;
        for (const uint32& x : neigh_x) {
            for (const uint32& y : neigh_y) {
                to_1d(k_out, x, y);
                indices.push_back(k_out);
            }
        }
        return;
    }


};






/*
 ==============================================================================
 ==============================================================================
 ==============================================================================
 ==============================================================================
 */

/*
 Sample locations from a vector of probabilities, where probabilities
 get updated as samples are chosen.
 For one plant type (virus- or Pseudomonas-infected)
 */
class LocationSampler {

    // This is to avoid infinite weights.
    // It's still a massive number (1.797693e+302 on my machine):
    const double max_wt = std::numeric_limits<double>::max() / 1e6;

    arma::vec weights;
    arma::vec cs_probs;
    uint32 n;

    bool needs_recalc;

    // Make `cs_probs` into a vector that's the cumulative sum of
    // weights / sum(weights). The last value in `cs_probs` should always be 1.
    void calc_cumsum() {
        double p_sum = arma::accu(weights);
        if (p_sum <= 0) {
            weights.ones();
            p_sum = weights.n_elem;
        }
        cs_probs(0) = weights(0) / p_sum;
        for (uint32 i = 1; i < n; i++) {
            cs_probs(i) = cs_probs(i-1) + weights(i) / p_sum;
        }
        needs_recalc = false;
        return;
    }

public:

    LocationSampler(const uint32& n_)
        : weights(n_, arma::fill::ones),
          cs_probs(n_, arma::fill::none),
          n(n_),
          needs_recalc(true) {
        calc_cumsum();
    }

    /*
     Update weights by multiplying by a given value, but only if the weight
     wasn't already changed to zero (which happens if its location was
     already sampled).
     These functions do NOT re-calculate cumulative sums because
     that is done later inside `sample`.
     They both also change `needs_recalc` to true if they actually change
     one or more probabilities.
     */
    void update_weights(const std::vector<uint32>& indices,
                        const double& wt_val) {
        if (wt_val == 1) return;
        uint32 n_changed = 0;
        bool not_exceeded;
        for (const uint32& k : indices) {
            /*
             Only change weight if it's greater than zero and either
             the weight will decrease or the weight will increase but
             hasn't already exceeded the maximum allowed value.
             */
            not_exceeded = weights(k) < max_wt && wt_val > 1;
            if (weights(k) > 0 && (not_exceeded || wt_val < 1)) {
                weights(k) *= wt_val;
                if (weights(k) > max_wt) weights(k) = max_wt;
                n_changed++;
            }
        }
        if (n_changed > 0) needs_recalc = true;
        return;
    }
    void update_weights(const uint32& k,
                        const double& wt_val) {
        bool not_exceeded = weights(k) < max_wt && wt_val > 1;
        if (wt_val != 1 && weights(k) > 0 && (not_exceeded || wt_val < 1)) {
            weights(k) *= wt_val;
            if (weights(k) > max_wt) weights(k) = max_wt;
            needs_recalc  = true;
        }
        return;
    }

    // Check to see if probabilities needs re-calculated, then
    // do weighted sampling for an index from 0 to (n-1):
    uint32 sample(pcg32& eng) {
        if (needs_recalc) calc_cumsum();
        double u = runif_01(eng);
        uint32 k = 0;
        while (k < n && cs_probs(k) < u) k++;
        return k;
    }



};







/*
 ==============================================================================
 ==============================================================================
 ==============================================================================
 ==============================================================================
 */

/*
 Do all sampling for a single landscape.
 */
class LandSimmer {

    arma::mat wt_mat;
    std::vector<uint32> n_samples;
    uint32 n_plants;
    uint32 x_size;
    uint32 y_size;

    // One location sampler for each type (virus and Pseudomonas):
    std::vector<LocationSampler> samplers;

    // Convert coordinates between 1D and 2D:
    DimensionConverter dim_conv;

    // Object collecting sampled type for each plant:
    std::vector<uint32> out_types;

    // Random number generator
    pcg32 eng;


    // Use bit assignment to quickly assign types:
    void put_type(const uint32& type, const uint32& k) {
        out_types[k] = out_types[k] | ((uint32)1 << type);
        return;
    }


public:

    LandSimmer(const arma::mat& wt_mat_,
               const uint32& n_virus_,
               const uint32& n_pseudo_,
               const uint32& x_size_,
               const uint32& y_size_)
    : wt_mat(wt_mat_),
      n_samples({n_virus_, n_pseudo_}),
      n_plants(x_size_ * y_size_),
      x_size(x_size_),
      y_size(y_size_),
      samplers(2, LocationSampler(n_plants)),
      dim_conv(x_size, y_size),
      out_types(n_plants, 0U),
      eng() {
        seed_pcg(eng);
      }


    void run(RcppThread::ProgressBar& prog_bar, const bool& show_progress) {

        uint32 total_samps = std::accumulate(n_samples.begin(), n_samples.end(), 0U);
        // Because `put_type` is assigning individual bits, and `uint32` only
        // have 32 bits:
        if (n_samples.size() > 32)
            stop("INTERNAL ERROR: Cannot bit-assign with >32 types");

        // Vector of vector types that will be shuffled to avoid having one
        // type always sampled first:
        std::vector<uint32> type_to_samp;
        type_to_samp.reserve(total_samps);
        for (uint32 i = 0; i < n_samples.size(); i++) {
            for (uint32 j = 0; j < n_samples[i]; j++) {
                type_to_samp.push_back(i);
            }
        }
        // Fast shuffle:
        for (uint32 i = type_to_samp.size(); i > 1; i--) {
            uint32 j = runif_01(eng) * i;
            std::swap(type_to_samp[i-1], type_to_samp[j]);
        }
        if (type_to_samp.size() != total_samps)
            stop("INTERNAL ERROR: type_to_samp.size() != total_samps");

        uint32 x, y, k;
        std::vector<uint32> neighbors;
        neighbors.reserve(9); // highest number of neighbors possible
        uint32 n = 0; // for error checking


        for (const uint32& i : type_to_samp) {

            k = samplers[i].sample(eng);
            dim_conv.to_2d(x, y, k); // assign new x and y based on k

            // add to output:
            put_type(i, k);

            // Adjust sampling probabilities:
            dim_conv.get_neighbors(neighbors, k); // fill neighbors vector
            for (uint32 j = 0; j < n_samples.size(); j++) {
                samplers[j].update_weights(neighbors, wt_mat(i,j));
            }
            /*
             Note: You don't have to update `k`th prob to zero after the call
             to `update_weights` on `neighbors` even though the latter also
             updates `k` because `update_weights` never updates weights that
             are already set to zero.
             */
            samplers[i].update_weights(k, 0.0);

            n++;

            if (show_progress) prog_bar++;
            if (n % 10 == 0) RcppThread::checkUserInterrupt();

        }

        return;
    }


    // Fill an output cube with the matrix for this landscape's output.
    // `s` refers to the slice index for this landscape
    void fill_output(arma::ucube& out, const uint32& s) {

        uint32 x, y;

        for (uint32 k = 0; k < out_types.size(); k++) {
            dim_conv.to_2d(x, y, k);
            out(x, y, s) = out_types[k];
        }

        return;

    }



};



#endif
