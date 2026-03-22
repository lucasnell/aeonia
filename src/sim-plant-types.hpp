#ifndef __AEONIA_LANDSCAPES_H
#define __AEONIA_LANDSCAPES_H

#include "aeonia_types.hpp"


#include <RcppArmadillo.h>
#include <vector>
#include <deque>
#include <string>
#include <limits>
#include <numeric>
#include <pcg/pcg_random.hpp>   // pcg prng
#include <RcppThread.h>         // multithreading


#include "convert-dims.hpp"     // DimensionConverter class
#include "util.hpp"             // thread_check fxn
#include "pcg.hpp"              // runif_ab fxn


using namespace Rcpp;


/*
 Convert two Nullable NumericMatrix objects of starting virus and pseudomonas
 xy positions to a single object for use in OnePlantTypeSimmer constructor.
 If both `starts` are NULL, then an empty `arma::mat` is returned.
 */
arma::umat get_xy_starts(const Rcpp::Nullable<IntegerMatrix>& virus_starts,
                         const Rcpp::Nullable<IntegerMatrix>& pseudo_starts,
                         const uint32& n_virus,
                         const uint32& n_pseudo,
                         const uint32& n_x,
                         const uint32& n_y) {
    arma::imat tmp1, tmp2;
    if (virus_starts.isNotNull()) {
        tmp1 = Rcpp::as<arma::imat>(virus_starts);
        if (tmp1.n_cols != 2) stop("ncol(virus_starts) must be 2");
        if (tmp1.n_rows > n_virus) stop("nrow(virus_starts) must be <= n_virus");
        if (arma::any(arma::vectorise(tmp1) < 1)) stop("items in virus_starts must be >= 1");
        if (arma::any(tmp1.col(0) > n_x)) stop("items in virus_starts[,1] must be <= n_x");
        if (arma::any(tmp1.col(1) > n_y)) stop("items in virus_starts[,2] must be <= n_y");
    }
    if (pseudo_starts.isNotNull()) {
        tmp2 = Rcpp::as<arma::imat>(pseudo_starts);
        if (tmp2.n_cols != 2) stop("ncol(pseudo_starts) must be 2");
        if (tmp2.n_rows > n_pseudo) stop("nrow(pseudo_starts) must be <= n_pseudo");
        if (arma::any(arma::vectorise(tmp2) < 1)) stop("items in pseudo_starts must be >= 1");
        if (arma::any(tmp2.col(0) > n_x)) stop("items in pseudo_starts[,1] must be <= n_x");
        if (arma::any(tmp2.col(1) > n_y)) stop("items in pseudo_starts[,2] must be <= n_y");
    }

    arma::umat xy0(tmp1.n_rows + tmp2.n_rows, 3U, arma::fill::none);
    for (uint32 i = 0; i < tmp1.n_rows; i++) {
        xy0(i, 0) = (uint32)0;
        xy0(i, 1) = static_cast<uint32>(tmp1(i, 0)) - (uint32)1;
        xy0(i, 2) = static_cast<uint32>(tmp1(i, 1)) - (uint32)1;
    }
    uint32 j;
    for (uint32 i = 0; i < tmp2.n_rows; i++) {
        j = tmp1.n_rows + i;
        xy0(j, 0) = (uint32)1;
        xy0(j, 1) = static_cast<uint32>(tmp2(i, 0)) - (uint32)1;
        xy0(j, 2) = static_cast<uint32>(tmp2(i, 1)) - (uint32)1;
    }

    return xy0;
}




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
        cs_probs.at(0) = weights.at(0) / p_sum;
        for (uint32 i = 1; i < n; i++) {
            cs_probs.at(i) = cs_probs.at(i-1) + weights.at(i) / p_sum;
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
            not_exceeded = weights.at(k) < max_wt && wt_val > 1;
            if (weights.at(k) > 0 && (not_exceeded || wt_val < 1)) {
                weights.at(k) *= wt_val;
                if (weights.at(k) > max_wt) weights.at(k) = max_wt;
                n_changed++;
            }
        }
        if (n_changed > 0) needs_recalc = true;
        return;
    }
    void update_weights(const uint32& k,
                        const double& wt_val) {
        bool not_exceeded = weights.at(k) < max_wt && wt_val > 1;
        if (wt_val != 1 && weights.at(k) > 0 && (not_exceeded || wt_val < 1)) {
            weights.at(k) *= wt_val;
            if (weights.at(k) > max_wt) weights.at(k) = max_wt;
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
        while (k < n && cs_probs.at(k) < u) k++;
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
class OnePlantTypeSimmer {

    arma::mat wt_mat;
    arma::uvec n_samples;
    uint32 n_plants;
    uint32 x_size;
    uint32 y_size;

    // One location sampler for each type (virus and Pseudomonas):
    std::array<LocationSampler, 2U> samplers;

    // Convert coordinates between 1D and 2D:
    DimensionConverter dim_conv;

    // Object collecting sampled type for each plant:
    arma::uvec out_types;

    // vector for neighbors:
    std::vector<uint32> neighbors;

    // Random number generator
    pcg32 eng;


    // Use bit assignment to quickly assign types:
    inline void put_type(const uint32& type, const uint32& k) {
        out_types.at(k) = out_types.at(k) | ((uint32)1 << type);
        return;
    }

    /*
     Add to output and adjust both sampling probabilities and samplers.
     `i` is the type (0 for virus, 1 for Pseudomonas),
     `k` is the location where that type is being placed
     */
    void out_and_adjust(const uint32& i, const uint32& k) {

        // add to output:
        put_type(i, k);

        // Adjust sampling probabilities:
        dim_conv.nextdoor_neighbors(neighbors, k); // fill neighbors vector
        for (uint32 j = 0; j < n_samples.n_elem; j++) {
            samplers[j].update_weights(neighbors, wt_mat.at(i,j));
        }
        /*
         Note: You don't have to update `k`th prob to zero after the call
         to `update_weights` on `neighbors` even though the latter also
         updates `k` because `update_weights` never updates weights that
         are already set to zero.
         */
        samplers[i].update_weights(k, 0.0);

        return;

    }



public:

    OnePlantTypeSimmer(const arma::mat& wt_mat_,
                       const uint32& n_virus_,
                       const uint32& n_pseudo_,
                       const uint32& x_size_,
                       const uint32& y_size_,
                       const arma::umat& virus_pseudo_xy0)
    : wt_mat(wt_mat_),
      n_samples({n_virus_, n_pseudo_}),
      n_plants(x_size_ * y_size_),
      x_size(x_size_),
      y_size(y_size_),
      samplers({LocationSampler(n_plants), LocationSampler(n_plants)}),
      dim_conv(x_size, y_size),
      out_types(n_plants, arma::fill::zeros),
      neighbors(),
      eng() {

        seed_pcg(eng);
        neighbors.reserve(9); // highest number of neighbors possible

        if (!virus_pseudo_xy0.is_empty()) {

            uint32 k;

            for (uint32 r = 0; r < virus_pseudo_xy0.n_rows; r++) {

                const uint32& i(virus_pseudo_xy0(r, 0));
                const uint32& x(virus_pseudo_xy0(r, 1));
                const uint32& y(virus_pseudo_xy0(r, 2));

                // assign new k based on x and y:
                dim_conv.to_1d(k, x, y);
                // Add to output and adjust both sampling probs and samplers
                out_and_adjust(i, k);
                // Adjust # samples:
                n_samples.at(i)--;

            }
        }

    }


    void run(RcppThread::ProgressBar& prog_bar, const bool& show_progress) {

        uint32 total_samps = arma::accu(n_samples);
        // Because `put_type` is assigning individual bits, and `uint32` only
        // have 32 bits:
        if (n_samples.n_elem > 32)
            stop("INTERNAL ERROR: Cannot bit-assign with >32 types");

        // Vector of vector types that will be shuffled to avoid having one
        // type always sampled first:
        std::vector<uint32> type_to_samp;
        type_to_samp.reserve(total_samps);
        for (uint32 i = 0; i < n_samples.n_elem; i++) {
            for (uint32 j = 0; j < n_samples.at(i); j++) {
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

        // uint32 x, y;
        uint32 k;
        uint32 n = 0; // for error checking

        for (const uint32& i : type_to_samp) {

            k = samplers[i].sample(eng);
            // dim_conv.to_2d(x, y, k); // assign new x and y based on k

            // Add to output and adjust both sampling probs and samplers
            out_and_adjust(i, k);

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

        for (uint32 k = 0; k < out_types.n_elem; k++) {
            dim_conv.to_2d(x, y, k);
            out.at(x, y, s) = out_types.at(k);
        }

        return;

    }



};



#endif
