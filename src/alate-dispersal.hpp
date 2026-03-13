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
#include "convert-dims.hpp"     // DimensionConverter, XY, and get_bit_bool
#include "aphids.hpp"        // AphidPop class
// #include "one-plant.hpp"        // OnePlant class
#include "pcg.hpp"              // pcg type, runif_01 fxn



using namespace Rcpp;






class AliasSampler {
public:
    AliasSampler() : W0(), Prob(), Alias(), n(0) {};
    AliasSampler(const std::vector<double>& weights)
        : W0(weights), Prob(weights.size()), Alias(weights.size()), n(weights.size()) {
        construct();
    }
    AliasSampler(arma::vec weights)
        : W0(arma::conv_to<std::vector<double>>::from(weights)),
          Prob(weights.n_elem), Alias(weights.n_elem), n(weights.n_elem) {
        construct();
    }
    // Copy constructor
    AliasSampler(const AliasSampler& other)
        : W0(other.W0), Prob(other.Prob), Alias(other.Alias), n(other.n) {}

    // Actual alias sampling
    inline uint32 sample(pcg32& eng) const {
        // Fair dice roll from n-sided die
        uint32 i = runif_01(eng) * n;
        // uniform in range (0,1)
        double u = runif_01(eng);
        if (u < Prob[i]) return(i);
        return Alias[i];
    };

    // Reconstruct using an entirely new vector:
    void reconstruct(const std::vector<double>& new_wts) {
        n = new_wts.size();
        W0 = new_wts;
        construct();
        return;
    }


private:

    std::vector<double> W0; // weights to start
    std::vector<double> Prob;
    std::vector<uint32> Alias;
    uint32 n;


    void construct() {

        if (Prob.size() != n) Prob.resize(n);
        if (Alias.size() != n) Alias.resize(n);

        // make sure they sum to `n`:
        double mult = static_cast<double>(n) / std::accumulate(W0.begin(), W0.end(), 0.0);
        std::vector<double> p;
        p.reserve(W0.size());
        for (double& w : W0) p.push_back(w * mult);

        std::deque<uint32> Small;
        std::deque<uint32> Large;
        for (uint32 i = 0; i < n; i++) {
            if (p[i] < 1.0) {
                Small.push_back(i);
            } else Large.push_back(i);
        }

        uint32 l, g;
        while (!Small.empty() && !Large.empty()) {
            l = Small.front();
            Small.pop_front();
            g = Large.front();
            Large.pop_front();
            Prob[l] = p[l];
            Alias[l] = g;
            p[g] = (p[g] + p[l]) - 1.0;
            if (p[g] < 1.0) {
                Small.push_back(g);
            } else Large.push_back(g);
        }
        while (!Large.empty()) {
            g = Large.front();
            Large.pop_front();
            Prob[g] = 1.0;
        }
        while (!Small.empty()) {
            l = Small.front();
            Small.pop_front();
            Prob[l] = 1.0;
        }

        return;
    }
};




/*
 ==============================================================================*
 ==============================================================================*
 Info to create flight paths for alates.
 ==============================================================================*
 ==============================================================================*
 */
class AlateFlightInfo {

    DimensionConverter dim_conv;

    /*
     Logical for whether any plants have changed to infectious, and
     vector to keep track of which plants need to have their samplers updated.
     */
    bool any_changed;
    std::vector<bool> update_sampler;

    // Vector of which other plants alates can travel to from each plant:
    std::vector<std::vector<uint32>> neighbors;

    // Vector of sampling weights for each plant:
    std::vector<double> land_wts;

    // Alias sampler for each plant:
    std::vector<AliasSampler> samplers;

    double virus_attract;
    double pseudo_repel;
    double epsilon;
    double w;

    uint32 max_fly_t;

    // Bounds of landscape:
    uint32 n_x;
    uint32 n_y;
    uint32 n_plants; // Total plants (n_x * n_y)
    uint32 n_neigh; // Max neighbors



    /*
     Fill statuses (`virus` and `pseudo`), plus `samplers` and related objects
     (notably `land_wts` and `neighbors`).
     */
    void fill_status_samplers(const arma::umat& landscape_, const double& radius);

    // Sample new location if alate doesn't stay to feed.
    // Assumes that `virus` and `pseudo` fields are already set.
    uint32 sample(const uint32& k, pcg32& eng) {
        return neighbors[k][samplers[k].sample(eng)];
    }

    inline void sample_inoculation(const double& p_load_alate,
                                   const double& p_load_plant,
                                   double& u,
                                   bool& has_virus,
                                   bool& infectious,
                                   bool& exposed,
                                   uint32& exp_days,
                                   pcg32& eng);


public:


    // virus (infectious only) and pseudomonas presence for each plant in
    // the landscape:
    vMatrix<bool> virus;
    vMatrix<bool> pseudo;



    AlateFlightInfo(const uint32& max_fly_t_,
                    const arma::umat& landscape_,
                    const double& radius,
                    const double& virus_attract_,
                    const double& pseudo_repel_,
                    const double& epsilon_,
                    const double& w_)
        : dim_conv(landscape_.n_rows, landscape_.n_cols),
          any_changed(false),
          update_sampler(landscape_.n_elem, false),
          neighbors(),
          land_wts(),
          samplers(),
          virus_attract(virus_attract_),
          pseudo_repel(pseudo_repel_),
          epsilon(epsilon_),
          w(w_),
          max_fly_t(max_fly_t_),
          n_x(landscape_.n_rows),
          n_y(landscape_.n_cols),
          n_plants(landscape_.n_elem),
          n_neigh(),
          virus(),
          pseudo() {

        fill_status_samplers(landscape_, radius);

    }

    // Let this object know that a plant was newly infected so that it can
    // update the landscape (`virus` and `land_wts`)
    // and let samplers know to update:
    void newly_infected(const uint32& x, const uint32& y) {
        virus(x,y) = true;
        uint32 k = dim_conv.to_1d(x, y);
        land_wts[k] *= virus_attract;
        any_changed = true;
        for (const uint32& l : neighbors[k]) update_sampler[l] = true;
        return;
    }




    /*
     Have all alates across the landscape fly and eventually settle.
     It also samples for whether virus spread happens, and adjusts the
     `OnePlant` objects accordingly.
     */
    void infest(const double& p_load_alate,
                const double& p_load_plant,
                std::vector<XY>& alate_plants,
                vMatrix<bool>& exposed,
                vMatrix<bool>& infectious,
                vMatrix<uint32>& exp_days,
                vMatrix<AphidPop>& aphids,
                vMatrix<arma::uvec>& n_alates,
                arma::umat& dispersals,
                pcg32& eng);


    void to_2d(uint32& x, uint32& y, const uint32& k) const {
        dim_conv.to_2d(x, y, k);
        return;
    }
    void to_1d(uint32& k, const uint32& x, const uint32& y) const {
        dim_conv.to_1d(k, x, y);
        return;
    }


};





#endif
