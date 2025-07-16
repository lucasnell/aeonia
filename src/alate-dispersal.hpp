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
#include "convert-dims.hpp"     // XY and get_bit_bool
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;






/*
 Info to create flight paths for alates.
 */
class AlateFlightInfo {

    uint32 max_fly_t;

    arma::imat neigh_dxdy;

    double alpha;
    double beta;
    double epsilon;
    double w;

    arma::vec weights;      // sampling weights
    arma::vec cs_probs;     // cumulative sum of sampling probabilities

    // Bounds of landscape:
    uint32 n_x;
    uint32 n_y;


    /*
     Iterate one time step (after the first one). Adjusts `x` and `y` for new
     coordinates, and returns true if sims should keep going bc alate has
     not settled.
     It also adds new coordinates to `path`.
     */
    bool iterate(uint32& x, uint32& y, pcg32& eng) {

        // Sample for whether alate will stay to feed at this plant:
        double feed_p = w;
        if (virus[x][y]) feed_p *= epsilon;
        if (runif_01(eng) < feed_p) {
            return false;
        }

        // Sample for new location (this function also updates x and y):
        sample_location(x, y, eng);

        // add new coordinates to `path`:
        path.push_back(XY(x, y));

        return true;

    }

    /*
     Do first iteration where the alate ALWAYS leaves the plant.
     This function also clears `path` for a fresh start.
     */
    void first_iterate(uint32& x, uint32& y, pcg32& eng) {
        sample_location(x, y, eng);
        path.clear();
        path.push_back(XY(x, y));
        return;
    }

    // Sample new location if alate doesn't stay to feed.
    // Assumes that `virus` and `pseudo` fields are already set.
    void sample_location(uint32& x, uint32& y, pcg32& eng) {

        // First calculate new sampling weights based on current location:
        double wt_tmp;
        double wt_sum = 0;
        uint32 xi, yi;
        for (uint32 i = 0; i < weights.n_elem; i++) {
            const int32& dx(neigh_dxdy(i,0));
            const int32& dy(neigh_dxdy(i,1));
            xi = x+dx;
            yi = y+dy;
            /*
             These first two checks prevent (1) integer overflow and
             (2) going past bounds of landscape matrix.
             The checks need to happen in this order because if the first check
             is true, `xi` and `yi` are massive numbers near 4e9.
             I'm marking these weights as negative in case the weights sum
             to zero. See below (immediately after this for loop) for more info.
             */
            if ((dx < 0 && std::abs(dx) > x) || (dy < 0 && std::abs(dy) > y)) {
                weights(i) = -1.0;
            } else if ((x + dx) >= n_x || (y + dy) >= n_y) {
                weights(i) = -1.0;
            } else {
                wt_tmp = 0.0;
                if (virus[xi][yi]) wt_tmp += alpha;
                if (pseudo[xi][yi]) wt_tmp += beta;
                weights(i) = std::exp(wt_tmp);
                wt_sum += weights(i);
            }
        }
        /*
         For very low values of alpha or beta (e.g., -1e6), weights can be
         zero. If all weights sum to zero, then we have to choose a
         plant randomly but don't want to choose one outside the landscape!
         That's why we marked locations outside the landscape as negative.
         Here, if sum(weights) == 0, then we set all weights equal to zero to
         one and ignore the ones set to negative values.
         In all cases, we set the negative values to zero after the check.
         */
        if (wt_sum <= 0) {
            for (double& w : weights) {
                if (w == 0) {
                    w = 1;
                    wt_sum += 1;
                } else w = 0;
            }
        } else {
            for (double& w : weights) {
                if (w < 0) w = 0;
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

    // virus and pseudomonas presence for each plant in the landscape:
    std::vector<std::vector<bool>> virus;
    std::vector<std::vector<bool>> pseudo;
    std::vector<XY> path; // flight path to be accessed after running `fly`

    AlateFlightInfo(const uint32& max_fly_t_,
                    const arma::umat& landscape_,
                    const double& radius,
                    const double& alpha_,
                    const double& beta_,
                    const double& epsilon_,
                    const double& w_)
        : max_fly_t(max_fly_t_),
          neigh_dxdy(),
          alpha(alpha_),
          beta(beta_),
          epsilon(epsilon_),
          w(w_),
          weights(),
          cs_probs(),
          n_x(landscape_.n_rows),
          n_y(landscape_.n_cols),
          virus(),
          pseudo(),
          path() {

        path.reserve(max_fly_t+1U);

        // Fill `virus` and `pseudo` landscapes:
        virus.reserve(n_x);
        pseudo.reserve(n_x);
        bool virus_xy, pseudo_xy;
        for (uint32 x = 0; x < n_x; x++) {
            virus.push_back(std::vector<bool>());
            pseudo.push_back(std::vector<bool>());
            std::vector<bool>& virus_x(virus.back());
            std::vector<bool>& pseudo_x(pseudo.back());
            virus_x.reserve(n_y);
            pseudo_x.reserve(n_y);
            for (uint32 y = 0; y < n_y; y++) {
                virus_xy = get_bit_bool(0U, landscape_(x, y));
                pseudo_xy = get_bit_bool(1U, landscape_(x, y));
                virus_x.push_back(virus_xy);
                pseudo_x.push_back(pseudo_xy);
            }
        }

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


    /*
     Have this alate fly and eventually settle.
     The higher-level function can then use the public `path` field from
     this object to access the plants the alate visited.
     Note: the `path` vector does NOT include the plant the alate started at.
     */
    void fly(const uint32& x0,
             const uint32& y0,
             pcg32& eng) {
        uint32 x = x0;
        uint32 y = y0;
        // Do first iteration where alate has to leave plant:
        first_iterate(x, y, eng);
        // Do remaining iterations:
        bool keep_going = false;
        uint32 t = 1;
        while (t < max_fly_t && keep_going) {
            keep_going = this->iterate(x, y, eng);
            t++;
        }
        return;
    }



};





#endif
