#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng


#include "aeonia_types.hpp"     // integer types
#include "plantscape.hpp"       // PlantScape class
#include "plantscape-outputs.hpp" // Functions to create outputs
#include "pcg.hpp"              // mt_seeds fxn
#include "util.hpp"             // thread_check fxn



using namespace Rcpp;



// Fill adult, female parasitoid densities by patch
// (stored in `Yi_mat`):
void PlantScape::fill_Yi_mat() {

    const double& zeta(wasps.zeta);
    const double& Y(wasps.Y);

    double z_tot = 0;
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            const AphidPop& aphids_xy(aphids[x][y]);
            z_mat(x,y) = aphids_xy.total();
            z_tot += z_mat(x,y);
        }
    }

    // Now go back through and calculate Y for each patch:
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            z_mat(x,y) /= z_tot;
            // Note: `wasp_attract` sums to 1 (verified inside `PlantScape`
            // constructor), and by default is the same for all patches
            Yi_mat(x,y) = Y * ((1-zeta) * wasp_attract(x,y) + zeta * z_mat(x,y));
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
bool PlantScape::iterate(arma::umat& dispersals) {

    alate_plants.clear();

    /*
     Go through once, calculating and extracting alates, and updating
     population dynamics and infectiousness
     (Infectiousness changes here due to plants transitioning
     from exposed to infectious.)
     */
    bool infectious0;
    double new_Y = 0;  // new adult parasitoids (male and female)
    double new_M, x_xy;
    arma::vec A_surv;
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {

            AphidPop& aphids_xy(aphids[x][y]);
            MummyPop& mummies_xy(mummies[x][y]);
            bool& exposed_xy(exposed[x][y]);
            bool& infectious_xy(infectious[x][y]);
            bool& pseudo_xy(pseudo[x][y]);
            uint32& exp_days_xy(exp_days[x][y]);
            const double& Yi(Yi_mat(x,y));
            arma::uvec& n_alates_xy(n_alates[x][y]);


            infectious0 = infectious_xy;

            // Set attack survival vector:
            x_xy = aphids_xy.total_unpar();
            wasps.A_mats(A_surv, aphids_xy.attack_surv, Yi, x_xy);

            // Now iterate aphids, then mummies
            new_M = aphids_xy.iterate(n_alates_xy, A_surv, eng);
            new_Y += mummies_xy.iterate(new_M); // also update new adult parasitoids

            if (exposed_xy) {
                exp_days_xy++;
                if (exp_days_xy >= total_exp_days) {
                    infectious_xy = true;
                    exposed_xy = false;
                    exp_days_xy = 0;
                }
            }

            if (arma::accu(n_alates_xy) > 0U) alate_plants.push_back(XY(x,y));

            // If newly infectious, update landscape and let `flight` know
            // that samplers need to be updated:
            if (!infectious0 && infectious_xy) {
                flight.newly_infected(x, y);
            }
        }
    }

    // now update adult parasitoids:
    wasps.iterate(new_Y, eng);

    // Now go back through and simulate virus spread:
    flight.infest(p_load_alate, p_load_plant, alate_plants,
                  exposed, infectious, exp_days, aphids,
                  n_alates, dispersals, eng);


    // Lastly, go through and check whether all plants are infected:
    bool all_infected = true;
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            const bool& infectious_xy(infectious[x][y]);
            if (!infectious_xy) {
                all_infected = false;
                break;
            }
        }
        if (!all_infected) break;
    }


    // Fill adult, female parasitoid densities by patch
    // This is done at the end of each iteration (and at PlantScape construction)
    // so that densities summed by time and separate by plant are equivalent
    // in the output objects.
    fill_Yi_mat();



    return all_infected;

}





PlantScape::PlantScape(const arma::umat& landscape_,
                       const DiseaseDispersal& disp_dis,
                       const InsectPops& insects,
                       const arma::mat& N0,
                       const arma::mat& W0,
                       const arma::mat& M0,
                       const double& Y0,
                       const arma::mat& wasp_attract_,
                       const std::vector<uint64>& seeds)
    : n_x(landscape_.n_rows),
      n_y(landscape_.n_cols),
      flight(disp_dis.max_fly_t, landscape_, disp_dis.radius,
             disp_dis.virus_attract, disp_dis.pseudo_repel,
             disp_dis.epsilon, disp_dis.w),
      wasp_attract(wasp_attract_),
      p_load_alate(disp_dis.p_load_alate),
      p_load_plant(disp_dis.p_load_plant),
      total_exp_days(disp_dis.total_exp_days),
      n_alates(vMatSize({n_x, n_y})),
      alate_plants(),
      eng(),
      z_mat(arma::size(landscape_), arma::fill::none),
      Yi_mat(arma::size(landscape_), arma::fill::none),
      wasps(insects.wasps),
      aphids(vMatSize({n_x, n_y})),
      mummies(vMatSize({n_x, n_y})),
      exposed(vMatSize({n_x, n_y})),
      infectious(vMatSize({n_x, n_y})),
      pseudo(vMatSize({n_x, n_y})),
      exp_days(vMatSize({n_x, n_y})) {

    alate_plants.reserve(n_x * n_y);

    if (wasp_attract.n_rows != n_x || wasp_attract.n_cols != n_y) {
        stop("wasp_attract not correct dimensions");
    }
    // Make this sum to one:
    double wa_sum = arma::accu(wasp_attract);
    if (wa_sum != 1) wasp_attract /= wa_sum;

    seed_pcg(eng, seeds);

    // for n_alates:
    arma::uvec tmp(insects.aphids.alates.X.n_elem - insects.aphids.adult_age,
                   arma::fill::none);

    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            aphids[x][y] = insects.aphids;
            mummies[x][y] = insects.mummies;
            exposed[x][y] = false;
            exp_days[x][y] = 0;
            infectious[x][y] = get_bit_bool(0U, landscape_(x, y));
            pseudo[x][y] = get_bit_bool(1U, landscape_(x, y));
            AphidPop& aphids_xy(aphids[x][y]);
            if (!pseudo[x][y]) aphids_xy.pseudo_surv = 1.0;
            aphids_xy.refresh_abunds(N0(x,y), W0(x,y));
            mummies[x][y].refresh_abunds(M0(x,y));
            n_alates[x][y] = tmp;
        }
    }

    // Set initial wasp density
    wasps.refresh_abunds(Y0);

    // To allow initial call to fill_output() to fill Y by plant:
    fill_Yi_mat();


    return;

}














