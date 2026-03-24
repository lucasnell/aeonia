#ifndef __AEONIA_ALATE_DISPERSAL_H
#define __AEONIA_ALATE_DISPERSAL_H

/*
 This contains code for alate dispersal.
 */
#include "aeonia_types.hpp"     // integer types


#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "convert-dims.hpp"     // DimensionConverter, XY, and get_bit_bool
#include "aphids.hpp"        // AphidPop class
#include "pcg.hpp"              // pcg type, runif_01 fxn



using namespace Rcpp;






class AliasSampler {
public:
    AliasSampler() : W0(), Prob(), Alias(), n(0) {};
    AliasSampler(const std::vector<double>& weights)
        : W0(weights),
          Prob(weights.size(), arma::fill::none),
          Alias(weights.size(), arma::fill::none),
          n(weights.size()) {
        construct();
    }
    AliasSampler(const arma::vec& weights)
        : W0(weights),
          Prob(weights.n_elem, arma::fill::none),
          Alias(weights.n_elem, arma::fill::none),
          n(weights.n_elem) {
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
        if (u < Prob.at(i)) return(i);
        return Alias.at(i);
    };

    // Reconstruct using an entirely new vector:
    void reconstruct(const std::vector<double>& new_wts) {
        n = new_wts.size();
        if (W0.n_elem != n) W0.set_size(n);
        for (uint32 i = 0; i < n; i++) W0.at(i) = new_wts[i];
        construct();
        return;
    }
    void reconstruct(const arma::vec& new_wts) {
        n = new_wts.n_elem;
        W0 = new_wts;
        construct();
        return;
    }


private:

    arma::vec W0; // weights to start
    arma::vec Prob;
    arma::uvec Alias;
    uint32 n;


    void construct() {

        if (Prob.n_elem != n) Prob.set_size(n);
        if (Alias.n_elem != n) Alias.set_size(n);

        // make sure they sum to `n`:
        double mult = static_cast<double>(n) / arma::accu(W0);
        arma::vec p = W0 * mult;

        std::deque<uint32> Small;
        std::deque<uint32> Large;
        for (uint32 i = 0; i < n; i++) {
            if (p.at(i) < 1.0) {
                Small.push_back(i);
            } else Large.push_back(i);
        }

        uint32 l, g;
        while (!Small.empty() && !Large.empty()) {
            l = Small.front();
            Small.pop_front();
            g = Large.front();
            Large.pop_front();
            Prob.at(l) = p.at(l);
            Alias.at(l) = g;
            p.at(g) = (p.at(g) + p.at(l)) - 1.0;
            if (p.at(g) < 1.0) {
                Small.push_back(g);
            } else Large.push_back(g);
        }
        while (!Large.empty()) {
            g = Large.front();
            Large.pop_front();
            Prob.at(g) = 1.0;
        }
        while (!Small.empty()) {
            l = Small.front();
            Small.pop_front();
            Prob.at(l) = 1.0;
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
    std::vector<arma::uvec> neighbors_;

    // Vector of sampling weights for each plant:
    arma::vec land_wts;

    // Alias sampler for each plant:
    std::vector<AliasSampler> samplers_;

    // 2D versions of neighbors and samplers:
    Span2D<AliasSampler> samplers;
    Span2D<arma::uvec> neighbors;

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


    // Sample new location if alate doesn't stay to feed.
    // Assumes that `virus` and `pseudo` fields are already set.
    inline uint32 sample(const uint32& k, pcg32& eng) {
        return neighbors_[k].at(samplers_[k].sample(eng));
    }
    // Same for x and y input coordinates:
    inline uint32 sample(const uint32& x, const uint32& y, pcg32& eng) {
        return neighbors[x,y].at(samplers[x,y].sample(eng));
    }

    inline void sample_inoculation(const double& p_load_alate,
                                   const double& p_load_plant,
                                   double& u,
                                   bool& has_virus,
                                   uint16& infectious,
                                   uint16& exposed,
                                   uint32& exp_days,
                                   pcg32& eng);


public:


    // virus (infectious only) and pseudomonas presence for each plant in
    // the landscape:
    arma::Mat<uint16> virus;
    arma::Mat<uint16> pseudo;



    AlateFlightInfo(const uint32& max_fly_t_,
                    const arma::umat& landscape_,
                    const double& radius,
                    const double& virus_attract_,
                    const double& pseudo_repel_,
                    const double& epsilon_,
                    const double& w_);

    // Let this object know that a plant was newly infected so that it can
    // update the landscape (`virus` and `land_wts`)
    // and let samplers know to update:
    inline void newly_infected(const uint32& x, const uint32& y) {
        virus.at(x, y) = true;
        uint32 k = dim_conv.to_1d(x, y);
        land_wts.at(k) *= virus_attract;
        any_changed = true;
        for (const uint32& l : neighbors_[k]) update_sampler[l] = true;
        return;
    }




    /*
     Have all alates across the landscape fly and eventually settle.
     It also samples for whether virus spread happens, and adjusts the
     `OnePlant` objects accordingly.
     */
    void infest(const double& p_load_alate,
                const double& p_load_plant,
                std::vector<std::array<uint32, 2U>>& alate_plants,
                arma::Mat<uint16>& exposed,
                arma::Mat<uint16>& infectious,
                arma::umat& exp_days,
                Span2D<AphidPop>& aphids,
                arma::ucube& n_alates,
                arma::umat& dispersals,
                pcg32& eng);


    inline void to_2d(uint32& x, uint32& y, const uint32& k) const {
        dim_conv.to_2d(x, y, k);
        return;
    }
    inline void to_1d(uint32& k, const uint32& x, const uint32& y) const {
        dim_conv.to_1d(k, x, y);
        return;
    }


};





#endif
