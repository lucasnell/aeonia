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
#include "one-plant.hpp"        // OnePlant class
#include "alate-dispersal.hpp"  // AlateFlightInfo class
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;







class PlantScape {

    AlateFlightInfo flight;
    std::vector<std::vector<OnePlant>> plants;

    // Probability that an uninoculated alate is loaded with a virus if it
    // probes an infectious plant:
    double delta_a;
    // Probability that an uninfected plant is loaded with a virus if it
    // is probed by a virus-bearing aphid:
    double delta_p;

    // for storing numbers of alates per plant:
    arma::umat n_alates;
    // for storing indices of plants that have produced >=1 alate:
    std::vector<XY> alate_plants;

    uint32 n_x;
    uint32 n_y;
    pcg32 eng;

    // How to summarize output (if at all):
    std::string summ;
    // Max time points to simulate:
    uint32 max_t;



    /*
     ==========================================================================*
     Write the current state of the PlantScape to the `output` field.
     `summ` is for how to summarize output (if at all).
     `summ = "none"` isn't recommended for long time series of many plants!

     If `summ == "none"`, `n_x * n_y` rows are output for each time point.
     The columns are...
         1) plant x
         2) plant y
         3) plant infectious with virus (0 or 1)
         4) aphid (non-winged) density
         5) alate density
         6) predator density

     If `summ == "pseudo"`, two rows are output for each time point.
     The columns are...
         1) plant pseudo (0 or 1)
         2) total plants of this type
         3) number of plants infectious with virus
         4) total aphid (non-winged) density summed across all plants of this type
         5) total alate density summed across all plants of this type
         6) total predator density summed across all plants of this type

     If `summ == "all"`, only one row is output per time point.
     The columns are...
         1) number of plants infectious with virus
         2) total aphid (non-winged) density summed across all plants
         3) total alate density summed across all plants
         4) total predator density summed across all plants
     ==========================================================================*
     */
    void fill_output() {

        if (summ == "none") {
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
        } else if (summ == "pseudo") {

            output.push_back(arma::mat(2, 6, arma::fill::zeros));
            arma::mat& output_t(output.back());
            // Set value indicating Pseudomonas present:
            output_t(1,0) = 1;
            // Fill the rest:
            uint32 k;
            for (uint32 x = 0; x < n_x; x++) {
                for (uint32 y = 0; y < n_y; y++) {
                    const OnePlant& plant(plants[x][y]);
                    k = (plant.pseudo) ? 1UL : 0UL;
                    output_t(k,1) += 1.0;
                    output_t(k,2) += static_cast<double>(plant.infectious);
                    output_t(k,3) += plant.insects.A();
                    output_t(k,4) += plant.insects.W();
                    output_t(k,5) += plant.insects.P;
                }
            }

        } else if (summ == "all") {

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

        } else {

            stop("INTERNAL ERROR: `! summ %in% c('none', 'pseudo', 'all')`");

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

        /*
         Go through once, calculating and extracting alates, and updating
         population dynamics and infectiousness
         (Infectiousness changes here due to plants transitioning
          from exposed to infectious.)
         */
        bool infectious0;
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                OnePlant& plant_xy(plants[x][y]);
                infectious0 = plant_xy.infectious;
                plant_xy.iterate(n_alates(x,y), eng);
                if (n_alates(x,y) > 0) alate_plants.push_back(XY(x,y));
                // If newly infectious, update landscape and let `flight` know
                // that samplers need to be updated:
                if (!infectious0 && plant_xy.infectious) {
                    flight.newly_infected(x, y);
                }
            }
        }

        // Now go back through and simulate virus spread:
        flight.infest(delta_a, delta_p, alate_plants, plants, n_alates, eng);

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
               const std::string& summ_,
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
          summ(summ_),
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
