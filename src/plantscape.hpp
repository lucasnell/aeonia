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
          insects(insects_),
          total_exp_days(total_exp_days_) {
        if (!pseudo) insects.set_B(0.0);
    }

    // iterate and set number of alates moving from this plant:
    void iterate(uint32& n_alates, pcg32& eng) {
        insects.iterate(n_alates, eng);
        if (exposed) {
            exp_days++;
            if (exp_days >= total_exp_days) {
                infectious = true;
                exposed = false;
                exp_days = 0;
            }
        }
        return;
    }

    friend class PlantScape;


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

    // Whether to output separately by plant:
    bool out_by_plant;
    // Max time points to simulate:
    uint32 max_t;


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
                    u = runif_01(eng);
                    has_virus = u < delta_a;
                } else has_virus = false;
                for (uint32 i = 0; i < flight.path.size(); i++) {
                    const XY& coord(flight.path[i]);
                    OnePlant& new_plant(plants[coord.x][coord.y]);
                    /*
                     If alate is non-virus-bearing and plant is
                     infectious, sample for whether alate is inoculated:
                     */
                    if (!has_virus && new_plant.infectious) {
                        u = runif_01(eng);
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
                        u = runif_01(eng);
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
                plants[final_coord.x][final_coord.y].insects.winged_adults() += 1;

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

        if (out_by_plant) {
            output.push_back(arma::mat(n_x * n_y, 6, arma::fill::none));
            arma::mat& output_t(output.back());
            uint32 k = 0;
            for (uint32 x = 0; x < n_x; x++) {
                for (uint32 y = 0; y < n_y; y++) {
                    const OnePlant& plant(plants[x][y]);
                    output_t(k,0) = static_cast<double>(x+1U);
                    output_t(k,1) = static_cast<double>(y+1U);
                    output_t(k,2) = static_cast<double>(plant.infectious);
                    output_t(k,3) = plant.insects.A();
                    output_t(k,4) = plant.insects.W();
                    output_t(k,5) = plant.insects.P;
                    k++;
                }
            }
        } else {
            output.push_back(arma::mat(1, 4, arma::fill::zeros));
            arma::mat& output_t(output.back());
            for (uint32 x = 0; x < n_x; x++) {
                for (uint32 y = 0; y < n_y; y++) {
                    const OnePlant& plant(plants[x][y]);
                    output_t(0,0) += static_cast<double>(plant.infectious);
                    output_t(0,1) += plant.insects.A();
                    output_t(0,2) += plant.insects.W();
                    output_t(0,3) += plant.insects.P;
                }
            }
        }

        return;
    }




    /*
     ==========================================================================*
     Iterate for one time point, and return bool for whether all plants
     are infected.
     ==========================================================================*
     */
    bool iterate() {

        alate_plants.clear();

        // Go through once, calculating and extracting alates, and updating
        // population dynamics and infectiousness
        // (Infectiousness only changes here due to plants transitioning
        // from exposed to infectious.)
        bool infectious0;
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                OnePlant& plant_xy(plants[x][y]);
                infectious0 = plant_xy.infectious;
                plant_xy.iterate(n_alates(x,y), eng);
                if (n_alates(x,y) > 0) alate_plants.push_back(XY(x,y));
                // If newly infectious, update landscape:
                if (!infectious0 && plant_xy.infectious) {
                    flight.virus[x][y] = true;
                }
            }
        }

        // Now go back through and simulate virus spread:
        spread_virus();

        // Fill `output` with current conditions:
        fill_output();

        // Lastly, go through and check whether all plants are infected:
        bool all_infected = true;
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                const OnePlant& plant_xy(plants[x][y]);
                if (!plant_xy.infectious) {
                    all_infected = false;
                    break;
                }
            }
            if (!all_infected) break;
        }

        return all_infected;

    }



public:

    // Output object:
    std::vector<arma::mat> output;


    PlantScape(const uint32& max_t_,
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
          out_by_plant(out_by_plant_),
          max_t(max_t_),
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
                    insects.set_aphids(A0(0,0), W0(0,0));
                    insects.P = P0(0,0);
                } else {
                    insects.set_aphids(A0(x,y), W0(x,y));
                    insects.P = P0(x,y);
                }
                if (!pseudo) insects.set_B(0.0);
            }
        }

        // Reserve max memory required:
        output.reserve(max_t_+1U);
        // fill starting conditions:
        fill_output();

    }


    /*
     ==========================================================================*
     Run this PlantScape:
     ==========================================================================*
     */
    void run(const bool& infect_stop,
             RcppThread::ProgressBar& prog_bar,
             const bool& show_progress) {
        bool all_infected = false;
        for (uint32 t = 0; t < max_t; t++) {
            all_infected = iterate();
            if (infect_stop && all_infected) break;
            if (show_progress) prog_bar++;
            if (t % 10 == 0) RcppThread::checkUserInterrupt();
        }
        return;
    }



};




#endif
