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
#include "insect-pops.hpp"      // AdultWaspPop class
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;







// Class to store densities by whatever grouping is desired for output.
// These values are shared no matter what type of output summary is used.
struct OutDensities {

    std::vector<double> virus;          // virus density (0 or 1 when this summarizes 1 patch)
    std::vector<double> aphids;         // aphids (non-winged) density
    std::vector<double> alates;         // alates density
    std::vector<double> parasitized;    // parasitized aphid density
    std::vector<double> mummies;        // mummy density
    double wasps;                       // adult, female parasitoid wasp density


    // Don't include any values, but reserve vectors for length `n`.
    OutDensities(const uint32& n)
        : virus(), aphids(), alates(), parasitized(), mummies(), wasps() {
        virus.reserve(n);
        aphids.reserve(n);
        alates.reserve(n);
        parasitized.reserve(n);
        mummies.reserve(n);
    };

    // Fill values for vectors of length `n`.
    // `wasps` double is simply set to `value`
    OutDensities(const uint32& n, const double& value)
        : virus(n, value), aphids(n, value), alates(n, value),
          parasitized(n, value), mummies(n, value),
          wasps(value) {};

    // Fill a single value for all vectors and `wasps` double
    OutDensities(const double& value)
        : virus(1, value), aphids(1, value), alates(1, value),
          parasitized(1, value), mummies(1, value),
          wasps(value) {};


    // Push back all values except wasps (bc wasps are at scale of all plants):
    void push_back(const double& virus_,
                   const double& aphids_,
                   const double& alates_,
                   const double& parasitized_,
                   const double& mummies_) {
        virus.push_back(virus_);
        aphids.push_back(aphids_);
        alates.push_back(alates_);
        parasitized.push_back(parasitized_);
        mummies.push_back(mummies_);
        return;
    }


    uint32 size() const {
        return virus.size();
    }

};






class PlantScape {

    AdultWaspPop wasps;
    AlateFlightInfo flight;
    std::vector<std::vector<OnePlant>> plants;
    // Attractiveness to parasitoids:
    arma::mat wasp_attract;

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
     Note also that `summ = "time"` and `summ = "all"` result in the same
     output here, but they get summarized differently in the `sim_plantscape`
     function.

     If `summ == "none"`, `n_x * n_y` rows are output for each time point.
     The columns are...
         1) plant x
         2) plant y
         3) plant infectious with virus (0 or 1)
         4) aphid (non-winged) density
         5) alate density
         6) parasitized aphid density
         7) mummy density
         8) parasitoid density

     If `summ == "pseudo"`, two rows are output for each time point.
     The columns are...
         1) plant pseudo (0 or 1)
         2) total plants of this type
         3) number of plants infectious with virus
         4) total aphid (non-winged) density summed across all plants of this type
         5) total alate density summed across all plants of this type
         6) total parasitized aphid density summed across all plants of this type
         7) total mummy density summed across all plants of this type
         8) total parasitoid density summed across all plants of this type

     If `summ %in% c("time", "all")`, only one row is output per time point.
     The columns are...
         1) number of plants infectious with virus
         2) total aphid (non-winged) density summed across all plants
         3) total alate density summed across all plants
         4) total parasitized aphid density summed across all plants
         5) total mummy density summed across all plants
         6) total parasitoid density summed across all plants
     ==========================================================================*
     */
    void fill_output() {

        if (summ == "none") {
            // Densities (all vectors except wasps start with length 0 but
            // are reserved for length n_x * n_y,
            // wasps are reserved for length 1):
            output_dens.push_back(OutDensities(n_x * n_y));
            OutDensities& output_dens_t(output_dens.back());
            // plant x, plant y:
            output_ids.push_back(arma::umat(n_x * n_y, 2, arma::fill::none));
            arma::umat& output_ids_t(output_ids.back());

            uint32 k = 0;
            output_dens_t.wasps = wasps.Y;
            for (uint32 x = 0; x < n_x; x++) {
                for (uint32 y = 0; y < n_y; y++) {

                    const OnePlant& plant(plants[x][y]);

                    output_ids_t(k,0) = x+1U;
                    output_ids_t(k,1) = y+1U;

                    output_dens_t.push_back(static_cast<double>(plant.infectious),
                                            plant.aphids.aphids(),
                                            plant.aphids.alates(),
                                            plant.aphids.P,
                                            plant.aphids.M);
                    k++;
                }
            }

        } else if (summ == "pseudo") {

            // Densities (wasps vector starts with length 1, all others start
            // with length 2, all values (including wasps) = 0):
            output_dens.push_back(OutDensities(2, 0.0));
            OutDensities& output_dens_t(output_dens.back());
            // plant pseudo (0 or 1), # of each type:
            output_ids.push_back(arma::umat({{0, 0},
                                             {1, 0}}));
            arma::umat& output_ids_t(output_ids.back());

            output_dens_t.wasps = wasps.Y;
            uint32 k;
            for (uint32 x = 0; x < n_x; x++) {
                for (const OnePlant& plant : plants[x]) {

                    k = (plant.pseudo) ? 1UL : 0UL;

                    output_ids_t(k,1) += 1U;

                    output_dens_t.virus[k] += static_cast<double>(plant.infectious);
                    output_dens_t.aphids[k] += plant.aphids.aphids();
                    output_dens_t.alates[k] += plant.aphids.alates();
                    output_dens_t.parasitized[k] += plant.aphids.P;
                    output_dens_t.mummies[k] += plant.aphids.M;

                    k++;
                }
            }

        } else if (summ == "time" || summ == "all") {

            // Densities (all vectors start with length 1, values = 0):
            output_dens.push_back(OutDensities(0.0));
            OutDensities& output_dens_t(output_dens.back());
            // (no ids here)

            output_dens_t.wasps = wasps.Y;

            for (uint32 x = 0; x < n_x; x++) {
                for (const OnePlant& plant : plants[x]) {
                    output_dens_t.virus[0] += static_cast<double>(plant.infectious);
                    output_dens_t.aphids[0] += plant.aphids.aphids();
                    output_dens_t.alates[0] += plant.aphids.alates();
                    output_dens_t.parasitized[0] += plant.aphids.P;
                    output_dens_t.mummies[0] += plant.aphids.M;
                }
            }

        } else {

            stop("INTERNAL ERROR: `! summ %in% c('none', 'pseudo', 'time', 'all')`");

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

        // Fill adult, female parasitoid densities by patch
        // (stored in `wasps.Yi_mat`):
        wasps.fill_Yi<OnePlant>(plants, wasp_attract);

        /*
         Go through once, calculating and extracting alates, and updating
         population dynamics and infectiousness
         (Infectiousness changes here due to plants transitioning
          from exposed to infectious.)
         */
        bool infectious0;
        double new_Y = 0;  // new adult parasitoids (male and female)
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                OnePlant& plant_xy(plants[x][y]);
                infectious0 = plant_xy.infectious;
                plant_xy.iterate(wasps.Yi_mat(x,y), n_alates(x,y), new_Y, eng);
                if (n_alates(x,y) > 0) alate_plants.push_back(XY(x,y));
                // If newly infectious, update landscape and let `flight` know
                // that samplers need to be updated:
                if (!infectious0 && plant_xy.infectious) {
                    flight.newly_infected(x, y);
                }
            }
        }
        // now update adult parasitoids:
        wasps.iterate(new_Y);

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

    // Output objects:
    std::vector<OutDensities> output_dens;  // densities of organism types
    std::vector<arma::umat> output_ids;     // identifiers for each set of densities


    PlantScape(const uint32& max_t_,
               const std::string& summ_,
               const uint32& max_fly_t_,
               const arma::umat& landscape_,
               const double& radius_,
               const double& alpha_,
               const double& beta_,
               const double& epsilon_,
               const double& w_,
               const arma::mat& wasp_attract_,
               const double& delta_a_,
               const double& delta_p_,
               const uint32& total_exp_days_,
               const InsectPops& insects,
               const arma::mat& N0,
               const arma::mat& W0,
               const double& Y0,
               const std::vector<uint64>& seeds)
        : wasps(insects.wasps),
          flight(max_fly_t_, landscape_, radius_, alpha_, beta_, epsilon_, w_),
          plants(),
          wasp_attract(wasp_attract_),
          delta_a(delta_a_),
          delta_p(delta_p_),
          n_alates(arma::size(landscape_), arma::fill::none),
          alate_plants(),
          n_x(landscape_.n_rows),
          n_y(landscape_.n_cols),
          eng(),
          summ(summ_),
          max_t(max_t_),
          output_dens(),
          output_ids() {

        wasps.Y = Y0;

        alate_plants.reserve(landscape_.n_elem);

        // Make this sum to one:
        double wa_sum = arma::accu(wasp_attract);
        if (wa_sum != 1) wasp_attract /= wa_sum;

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
                                            insects.aphids));
                AphidPops& aphids(plants_x.back().aphids);
                aphids.set_aphids(N0(x,y), W0(x,y));
                if (!pseudo) aphids.set_pseudo_surv(1.0);
            }
        }

        // Reserve max memory required:
        output_dens.reserve(max_t_+1U);
        if (summ == "none" || summ == "pseudo") output_ids.reserve(max_t_+1U);
        // fill starting conditions:
        fill_output();

    }

    // Adjust density dependence across plants:
    void adjust_K(const arma::mat& Kmat) {
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                plants[x][y].aphids.set_K(Kmat(x,y));
            }
        }
        return;
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
