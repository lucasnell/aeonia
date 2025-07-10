#ifndef __AEONIA_PLANTSCAPE_H
#define __AEONIA_PLANTSCAPE_H

/*
 This contains code for landscapes of plants, including
 plant disease dynamics and *Pseudomonas* presence.
 */

#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "aeonia_types.hpp"     // integer types
#include "convert-dims.hpp"     // XY and get_bit_bool
#include "insect-pops.hpp"      // InsectPops class
#include "alate-dispersal.hpp"  // AlateFlightInfo class
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;







/*
 ==============================================================================*
 ==============================================================================*
 OnePlant class
 ==============================================================================*
 ==============================================================================*
 */


struct OnePlant {

    // virus: exposed but not yet infectious?
    bool exposed;
    // virus: infectious?
    bool infectious;
    // contains Pseudomonas?
    bool pseudo;
    // days since exposure (ignored if not exposed):
    uint32 exp_days;

    // Insect populations:
    InsectPops insects;

    OnePlant(const bool& infectious_,
             const bool& pseudo_,
             const uint32& total_exp_days_,
             const InsectPops& insects_)
        : exposed(false),
          infectious(infectious_),
          pseudo(pseudo_),
          exp_days(0),
          total_exp_days(total_exp_days_),
          insects(insects_) {
        if (!pseudo) insects.set_B(0.0);
    }

    // iterate and output number of alates moving from this plant:
    uint32 iterate(pcg32& eng) {
        uint32 n_alates = insects.iterate(eng);
        if (exposed) {
            exp_days++;
            if (exp_days >= total_exp_days) {
                infectious = true;
                exposed = false;
                exp_days = 0;
            }
        }
        return n_alates;
    }



private:

    uint32 total_exp_days; // days since exposure required to transition to infected

};





/*
 ==============================================================================*
 ==============================================================================*
 PlantScape class
 ==============================================================================*
 ==============================================================================*
 */

class PlantScape {

    AlateFlightInfo flight;
    std::vector<std::vector<OnePlant>> plants;

    // Probability that an uninoculated alate is loaded with a virus if it
    // probes an infectious plant:
    const double& delta_a;
    // Probability that an uninfected plant is loaded with a virus if it
    // is probed by a virus-bearing aphid:
    const double& delta_p;

    // for storing numbers of alates per plant:
    arma::umat n_alates;
    // for storing indices of plants that have produced >=1 alate:
    std::vector<XY> alate_plants;

    uint32 n_x;
    uint32 n_y;
    pcg32 eng;

    // Index within `output` we're on and whether to output
    // separately by plant:
    uint32 out_idx;
    bool out_by_plant;


    /*
     ==========================================================================*
     After updating populations, infectiousness, and numbers of alates
     (including `alate_plants`, which stores coordinates for plants that
     have produced at least one alate), this function...
        1) Moves alates around in random order among plants to
           avoid strange patterns due to the order of plants
        2) Simulates virus-inoculation of alates and plants.
     ==========================================================================*
     */
    void spread_virus() {
        if (alate_plants.empty()) return;
        // Fast-shuffle `alate_plants` first:
        for (uint32 i = alate_plants.size(); i > 1; i--) {
            uint32 j = runif_01(eng) * i;
            std::swap(alate_plants[i-1], alate_plants[j]);
        }
        double u;
        bool has_virus;
        // Now go through in this random order to sample virus spread:
        for (const XY& alate_coords : alate_plants) {

            const uint32& x(alate_coords.x);
            const uint32& y(alate_coords.y);
            const uint32& n_alates_xy(n_alates(x,y));
            OnePlant& first_plant(plants[x][y]);

            for (uint32 a = 0; a < n_alates_xy; a++) {

                // Note: `flight.path` does NOT include the starting plant.
                flight.fly(x, y, eng);
                // If it starts on an infectious plant, then sample for whether
                // the alate is virus-bearing:
                if (first_plant.infectious) {
                    u = runif_01(rng);
                    has_virus = u < delta_a;
                } else has_virus = false;
                for (uint32 i = 0; i < flight.path.size(); i++) {
                    const XY& coord(flight.path[i]);
                    const OnePlant& new_plant(plants[coord.x][coord.y]);
                    /*
                     If alate is non-virus-bearing and plant is
                     infectious, sample for whether alate is inoculated:
                     */
                    if (!has_virus && new_plant.infectious) {
                        u = runif_01(rng);
                        if (u < delta_a) has_virus = true;
                    }
                    /*
                     If alate is virus-bearing and plant is not infectious
                     (or exposed), sample for whether plant is exposed.
                     If `total_exp_days` == 0, then it gets immediately
                     converted to infectious, in which case subsequent
                     alates can be inoculated by this plant the same day.
                     */
                    if (has_virus && !new_plant.infectious && !new_plant.exposed) {
                        u = runif_01(rng);
                        if (u < delta_p) {
                            if (new_plant.total_exp_days > 0U) {
                                new_plant.exposed = true;
                                new_plant.infectious = false;
                            } else {
                                new_plant.exposed = false;
                                new_plant.infectious = true;
                            }
                            new_plant.exp_days = 0U;
                        }
                    }
                }
                // Adding alate to winged population of plant settled on:
                const XY& final_coord(flight.path.back());
                plants[final_coord.x][final_coord.y].W += 1;

            }

        }

        return;

    }


    /*
     ==========================================================================*
     Write the current state of the PlantScape to the `output` field.
     `out_by_plant` is for whether to have separate output by plant.
     This isn't recommended for long time series of many plants!

     If `out_by_plant == true`, `n_x * n_y` rows are output for each time point.
     The columns are...
         1) plant x
         2) plant y
         3) plant infectious with virus (0 or 1)
         4) aphid (non-winged) density
         5) alate density
         6) predator density

     If `out_by_plant == false`, only one row is output per time point.
     The columns are...
         1) number of plants infectious with virus
         2) total aphid (non-winged) density summed across all plants
         3) total alate density summed across all plants
         4) total predator density summed across all plants
     ==========================================================================*
     */
    void fill_output() {

        if (out_idx >= output.n_slices)
            stop("INTERNAL ERROR: out_idx >= output.n_slices");

        if (out_by_plant) {
            uint32 k = 0;
            for (uint32 x = 0; x < n_x; x++) {
                for (uint32 y = 0; y < n_y; y++) {
                    const OnePlant& plant(plants[x][y]);
                    output(k,0,out_idx) = static_cast<double>(x+1U);
                    output(k,1,out_idx) = static_cast<double>(y+1U);
                    output(k,2,out_idx) = static_cast<double>(plant.infectious);
                    output(k,3,out_idx) = plant.insects.A;
                    output(k,4,out_idx) = plant.insects.W;
                    output(k,5,out_idx) = plant.insects.P;
                    k++;
                }
            }
        } else {
            for (uint32 x = 0; x < n_x; x++) {
                for (uint32 y = 0; y < n_y; y++) {
                    const OnePlant& plant(plants[x][y]);
                    output(0,0,out_idx) += static_cast<double>(plant.infectious);
                    output(0,1,out_idx) += plant.insects.A;
                    output(0,2,out_idx) += plant.insects.W;
                    output(0,3,out_idx) += plant.insects.P;
                }
            }
        }

        out_idx++;

        return;
    }




    /*
     ==========================================================================*
     Iterate for one time point:
     ==========================================================================*
     */
    void iterate() {

        alate_plants.clear();

        // Go through once, calculating and extracting alates, and updating
        // population dynamics and infectiousness
        bool infectious0;
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                OnePlant& plant_xy(plants[x][y]);
                infectious0 = plant_xy.infectious;
                n_alates(x,y) = plant_xy.iterate(pcg32& eng);
                if (n_alates(x,y) > 0) alate_plants.push_back(XY(x,y));
                // If newly infectious, update landscape:
                if (!infectious0 && plant_xy.infectious) {
                    // I do this by setting the first bit to 1:
                    bit_set1(0U, flight.landscape(x,y));
                }
            }
        }

        // Now go back through and simulate virus spread:
        spread_virus();

        // Fill `output` with current conditions:
        fill_output();

        return;

    }



public:

    // Output object:
    arma::cube output;


    PlantScape(const uint32& max_sim_t_,
               const bool& out_by_plant_,
               const uint32& max_fly_t_,
               const arma::umat& landscape_,
               const double& radius_,
               const double& alpha_,
               const double& beta_,
               const double& epsilon_,
               const double& w_,
               const double& delta_a_,
               const double& delta_p_,
               const uint32& total_exp_days_,
               const InsectPops& insects_,
               const arma::mat& A0,
               const arma::mat& W0,
               const arma::mat& P0,
               const std::vector<uint64>& seeds)
        : flight(max_fly_t_, landscape_, radius_, alpha_, beta_, epsilon_, w_),
          plants(),
          delta_a(delta_a_),
          delta_p(delta_p_),
          n_alates(arma::size(landscape_), arma::fill::none),
          alate_plants(),
          n_x(landscape_.n_rows),
          n_y(landscape_.n_cols),
          eng(),
          out_idx(0U),
          out_by_plant(out_by_plant_),
          output() {

        alate_plants.reserve(landscape_.n_elem);

        seed_pcg(eng, seeds);

        bool infectious, pseudo;

        plants.reserve(n_x);
        for (uint32 x = 0; x < n_x; x++) {
            plants.push_back(std::vector<OnePlant>());
            std::vector<OnePlant>& plants_x(plants.back());
            plants_x.reserve(n_y);
            for (uint32 y = 0; y < n_y; y++) {
                infectious = get_bit_bool(0U, landscape_(x, y));
                pseudo = get_bit_bool(1U, landscape_(x, y));
                plants_x.push_back(OnePlant(infectious, pseudo, total_exp_days_,
                                            insects_));
                InsectPops& insects(plants_x.back().insects);
                if (A0.n_elem == 1) {
                    insects.A = A0(0,0);
                    insects.W = W0(0,0);
                    insects.P = P0(0,0);
                } else {
                    insects.A = A0(x,y);
                    insects.W = W0(x,y);
                    insects.P = P0(x,y);
                }
            }
        }

        if (out_by_plant_) {
            output.set_size(n_x * n_y, 6, max_sim_t_+1U);
        } else {
            output.zeros(1, 4, max_sim_t_+1U);
        }
        // fill starting conditions:
        fill_output();

    }


    /*
     ==========================================================================*
     Run this PlantScape:
     ==========================================================================*
     */
    void run(RcppThread::ProgressBar& prog_bar, const bool& show_progress) {
        for (uint32 t = 1; t < output.n_slices; t++) {
            this->iterate();
            if (show_progress) prog_bar++;
            if (t % 10 == 0) RcppThread::checkUserInterrupt();
        }
        return;
    }



};




#endif
