
#include "plantscape.hpp"


using namespace Rcpp;



// Fill adult, female parasitoid densities by patch
// (stored in `Yi_mat`):
void PlantScape::fill_Yi_mat() {

    const double& zeta(wasps.zeta);
    const double& Y(wasps.Y);

    uint32 k = 0;
    double z_tot = 0;
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            const AphidPop& aphids_xy(aphids[k]);
            z_mat.at(x,y) = aphids_xy.total();
            z_tot += z_mat.at(x,y);
            k++;
        }
    }

    // Now go back through and calculate Y for each patch:
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            z_mat.at(x,y) /= z_tot;
            // Note: `wasp_attract` sums to 1 (verified inside `PlantScape`
            // constructor), and by default is the same for all patches
            Yi_mat.at(x,y) = Y * ((1-zeta) * wasp_attract.at(x,y) + zeta * z_mat.at(x,y));
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
    bool has_alates;

    uint32 k = 0;
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {

            AphidPop& aphids_xy(aphids[k]);
            MummyPop& mummies_xy(mummies[k]);
            k++;
            uint16& exposed_xy(exposed.at(x, y));
            uint16& infectious_xy(infectious.at(x, y));
            uint32& exp_days_xy(exp_days.at(x, y));
            const double& Yi(Yi_mat.at(x, y));
            arma::uvec n_alates_xy = n_alates.slice(y).unsafe_col(x);

            infectious0 = infectious_xy;

            // Set attack survival vector:
            x_xy = aphids_xy.total_unpar();
            wasps.A_mats(A_surv, aphids_xy.attack_surv, Yi, x_xy);

            // Now iterate aphids, then mummies
            new_M = aphids_xy.iterate(has_alates, n_alates_xy, A_surv, eng);
            new_Y += mummies_xy.iterate(new_M); // also update new adult parasitoids

            if (exposed_xy) {
                exp_days_xy++;
                if (exp_days_xy >= total_exp_days) {
                    infectious_xy = true;
                    exposed_xy = false;
                    exp_days_xy = 0;
                }
            }

            if (has_alates) alate_plants.push_back({x, y});

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
            const uint16& infectious_xy(infectious.at(x, y));
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
      n_plants(n_x * n_y),
      dim_conv(n_x, n_y),
      flight(disp_dis.max_fly_t, landscape_, disp_dis.radius,
             disp_dis.virus_attract, disp_dis.pseudo_repel,
             disp_dis.epsilon, disp_dis.w),
      wasp_attract(wasp_attract_),
      p_load_alate(disp_dis.p_load_alate),
      p_load_plant(disp_dis.p_load_plant),
      total_exp_days(disp_dis.total_exp_days),
      n_alates(insects.aphids.alates.X.n_elem - insects.aphids.adult_age,
               n_x, n_y, arma::fill::none),
      alate_plants(),
      eng(),
      z_mat(arma::size(landscape_), arma::fill::none),
      Yi_mat(arma::size(landscape_), arma::fill::none),
      wasps(insects.wasps),
      aphids(),
      mummies(),
      exposed(n_x, n_y, arma::fill::none),
      infectious(n_x, n_y, arma::fill::none),
      pseudo(n_x, n_y, arma::fill::none),
      exp_days(n_x, n_y, arma::fill::none) {


    alate_plants.reserve(n_x * n_y);

    if (wasp_attract.n_rows != n_x || wasp_attract.n_cols != n_y) {
        stop("wasp_attract not correct dimensions");
    }
    // Make this sum to one:
    double wa_sum = arma::accu(wasp_attract);
    if (wa_sum != 1) wasp_attract /= wa_sum;

    seed_pcg(eng, seeds);

    aphids.reserve(n_plants);
    mummies.reserve(n_plants);

    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            aphids.emplace_back(insects.aphids);
            mummies.emplace_back(insects.mummies);
            exposed.at(x, y) = false;
            exp_days.at(x, y) = 0;
            infectious.at(x, y) = get_bit_bool(0U, landscape_(x, y));
            pseudo.at(x, y) = get_bit_bool(1U, landscape_(x, y));
            AphidPop& aphids_xy(aphids.back());
            if (!pseudo.at(x, y)) aphids_xy.pseudo_surv = 1.0;
            aphids_xy.refresh_abunds(N0(x, y), W0(x, y));
            mummies.back().refresh_abunds(M0(x, y));
        }
    }

    // Set initial wasp density
    wasps.refresh_abunds(Y0);

    // To allow initial call to fill_output() to fill Y by plant:
    fill_Yi_mat();


    return;

}














