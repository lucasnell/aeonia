/*
 This contains code for alate dispersal.
 */


#include "alate-dispersal.hpp"


using namespace Rcpp;







/*
 Small function used in AlateFlightInfo::fill_samplers to create a matrix
 of dx and dy values for all neighbors within a given radius:
 */
arma::imat make_neigh_dxdy(const double& radius) {

    // Fill `neigh_dxdy` based on radius:
    uint32 total_rows = 0;
    int32 fl_radius = std::floor(radius);
    std::vector<arma::imat> dxdy_vec;
    dxdy_vec.reserve(fl_radius * 2U + 1U);
    int32 max_dy;
    uint32 rows_x;
    double radius2 = radius * radius;
    for (int32 dx = -fl_radius; dx <= fl_radius; dx++) {
        max_dy = std::floor(std::sqrt(radius2 - static_cast<double>(dx * dx)));
        rows_x = max_dy * 2U;
        if (dx != 0) rows_x++;
        dxdy_vec.emplace_back(rows_x, 2U, arma::fill::none);
        arma::imat& dxdy_i(dxdy_vec.back());
        uint32 i = 0;
        for (int32 dy = -max_dy; dy <= max_dy; dy++) {
            if (dy == 0 && dx == 0) continue;
            dxdy_i.at(i,0) = dx;
            dxdy_i.at(i,1) = dy;
            i++;
        }
        total_rows += dxdy_i.n_rows;
    }

    arma::imat neigh_dxdy(total_rows, 2U);
    uint32 k = 0;
    for (const arma::imat& dxdy : dxdy_vec) {
        for (uint32 j = 0; j < dxdy.n_rows; j++) {
            neigh_dxdy.at(k,0) = dxdy.at(j,0);
            neigh_dxdy.at(k,1) = dxdy.at(j,1);
            k++;
        }
    }

    return neigh_dxdy;
}



/*
 Fill statuses (`virus` and `pseudo`), plus `samplers` and related objects
 (notably `land_wts` and `neighbors`).
 */
AlateFlightInfo::AlateFlightInfo(const uint32& max_fly_t_,
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
      land_wts(landscape_.n_rows * landscape_.n_cols, arma::fill::none),
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
      virus(n_x, n_y, arma::fill::none),
      pseudo(n_x, n_y, arma::fill::none) {

    // Fill `virus` and `pseudo` landscapes:
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            virus.at(x, y) = get_bit_bool(0U, landscape_.at(x, y));
            pseudo.at(x, y) = get_bit_bool(1U, landscape_.at(x, y));
        }
    }


    // dx and dy for all neighbors within `radius`:
    arma::imat neigh_dxdy = make_neigh_dxdy(radius);
    n_neigh = neigh_dxdy.n_rows; // not true for plants on edge


    // First set `land_wts` and `neighbors`:
    neighbors.reserve(n_plants);
    uint32 ki, xi, yi;
    std::vector<uint32> neighbors_k;
    neighbors_k.reserve(n_neigh);
    uint32 k = 0;
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            // --------------------------
            // Set this plant's weight:
            land_wts.at(k) = 1.0;
            if (virus.at(x, y)) land_wts.at(k) *= virus_attract;
            if (pseudo.at(x, y)) land_wts.at(k) /= pseudo_repel;
            if (land_wts.at(k) == 0) {
                std::string err_msg = "\nERROR: virus_attract is very low or ";
                err_msg += "pseudo_repel is very high, ";
                err_msg += "resulting in sampling weights equal to zero. ";
                err_msg += "This results in computational problems.";
                throw std::runtime_error(err_msg.c_str());
            }
            // --------------------------
            // Set this plant's neighbors:
            neighbors_k.clear();
            for (uint32 i = 0; i < n_neigh; i++) {
                if (i >= neigh_dxdy.n_rows) throw std::runtime_error("Huh");
                const int32& dx(neigh_dxdy.at(i,0));
                const int32& dy(neigh_dxdy.at(i,1));
                xi = x+dx;
                yi = y+dy;
                /*
                 These first two checks prevent (1) integer overflow and
                 (2) going past bounds of landscape matrix.
                 The checks need to happen in this order because if the first
                 check is true, `xi` and `yi` are massive numbers near 4e9.
                 */
                if ((dx < 0 && std::abs(dx) > x) || (dy < 0 && std::abs(dy) > y)) {
                    continue;
                } else if (xi >= n_x || yi >= n_y) {
                    continue;
                } else {
                    dim_conv.to_1d(ki, xi, yi); // set `ki`
                    neighbors_k.push_back(ki);
                }
            }
            neighbors.emplace_back(neighbors_k);
            k++;
        }
    }


    // Now go back through and populate the samplers:
    samplers.reserve(n_plants);
    std::vector<double> weights_k;
    weights_k.reserve(n_neigh);
    for (uint32 k = 0; k < n_plants; k++) {
        weights_k.clear();
        for (const uint32& nk : neighbors[k]) {
            weights_k.push_back(land_wts.at(nk));
        }
        samplers.emplace_back(weights_k);
    }

    return;

}








/*
 Sample for a single potential inoculation.
 */
inline void AlateFlightInfo::sample_inoculation(const double& p_load_alate,
                                                const double& p_load_plant,
                                                double& u,
                                                bool& has_virus,
                                                uint16& infectious,
                                                uint16& exposed,
                                                uint32& exp_days,
                                                pcg32& eng) {


    /*
     If alate is non-virus-bearing and plant is
     infectious, sample for whether alate is inoculated:
     */
    if (!has_virus && infectious) {
        u = runif_01(eng);
        if (u < p_load_alate) has_virus = true;
    }
    /*
     If alate is virus-bearing and plant is not infectious
     (or exposed), sample for whether plant is exposed.
     */
    if (has_virus && !infectious && !exposed) {
        u = runif_01(eng);
        if (u < p_load_plant) {
            exposed = true;
            infectious = false;
            exp_days = 0U;
        }
    }

    return;

}







/*
 Have all alates across the landscape fly and eventually settle.
 It also samples for whether virus spread happens, and adjusts the
 `OnePlant` objects accordingly.
 Also updates `dispersals` if it has proper size.
 */
void AlateFlightInfo::infest(const double& p_load_alate,
                             const double& p_load_plant,
                             std::vector<std::array<uint32, 2U>>& alate_plants,
                             arma::Mat<uint16>& exposed,
                             arma::Mat<uint16>& infectious,
                             arma::umat& exp_days,
                             std::vector<AphidPop>& aphids,
                             arma::ucube& n_alates,
                             arma::umat& dispersals,
                             pcg32& eng) {

    if (alate_plants.empty()) return;

    // Store bool for whether to update `dispersals`:
    bool update_any_disps = dispersals.n_rows == n_plants || dispersals.n_rows == n_x;
    // Also store whether to record where alates came from:
    bool update_all_disps = dispersals.n_rows == n_plants;


    // Check for whether samplers need reconstructed and do it if needed:
    if (any_changed) {
        std::vector<double> new_wts;
        new_wts.reserve(n_neigh);
        for (uint32 k = 0; k < n_plants; k++) {
            if (update_sampler[k]) {
                new_wts.clear();
                for (const uint32& ki : neighbors[k]) {
                    new_wts.push_back(land_wts.at(ki));
                }
                samplers[k].reconstruct(new_wts);
                update_sampler[k] = false;
            }
        }
        any_changed = false;
    }

    // -----------------------------------------------------------------
    // Now do alate flying and virus spread:
    // -----------------------------------------------------------------
    // Fast-shuffle `alate_plants` first:
    for (uint32 i = alate_plants.size(); i > 1; i--) {
        uint32 j = runif_01(eng) * i;
        std::swap(alate_plants[i-1], alate_plants[j]);
    }
    double feed_p, u;
    bool has_virus;
    uint32 k0, k_old, x_old, y_old, k_new, x_new, y_new;
    bool keep_going;
    uint32 t;

    // Now go through in this random order to sample virus spread:
    for (const std::array<uint32, 2U>& coords0 : alate_plants) {

        const uint32& x0(coords0[0]);
        const uint32& y0(coords0[1]);
        k0 = dim_conv.to_1d(x0, y0);
        const AphidPop& aphids_xy(aphids[k0]);
        arma::uvec n_alates_xy = n_alates.slice(y0).unsafe_col(x0);

        uint32 adult_age = aphids_xy.adult_age;

        for (uint32 i = 0; i < n_alates_xy.n_elem; i++) {

            uint32& n_alates_xyi(n_alates_xy.at(i));

            while (n_alates_xyi > 0) {

                k_old = k0;
                x_old = x0;
                y_old = y0;

                // If it starts on an infectious plant, then sample for whether
                // the alate is virus-bearing:
                if (infectious.at(x0, y0)) {
                    u = runif_01(eng);
                    has_virus = u < p_load_alate;
                } else has_virus = false;

                // Sample for new location (alate always leaves first plant)
                k_new = sample(k_old, eng);
                dim_conv.to_2d(x_new, y_new, k_new);
                if (update_any_disps) {
                    if (update_all_disps) {
                        dispersals(k_new, k_old)++;
                    } else dispersals(x_new, y_new)++;
                }

                // Sample for whether aphid or plant is inoculated:
                sample_inoculation(p_load_alate, p_load_plant, u, has_virus,
                                   infectious.at(x_new, y_new),
                                   exposed.at(x_new, y_new),
                                   exp_days.at(x_new, y_new),
                                   eng);

                k_old = k_new;
                x_old = x_new;
                y_old = y_new;

                keep_going = true;
                t = 1;

                while (t < max_fly_t && keep_going) {

                    // Sample for whether alate will stay to feed at this plant:
                    feed_p = w;
                    if (virus.at(x_old, y_old)) feed_p *= epsilon;
                    if (runif_01(eng) < feed_p) {
                        keep_going = false;
                        // Adding alate to winged population of plant settled on:
                        AphidPop& aphids_old(aphids[k_old]);
                        aphids_old.alates.X.at(adult_age + i) += 1;
                        break;
                    }

                    // Sample for new location
                    k_new = sample(k_old, eng);
                    dim_conv.to_2d(x_new, y_new, k_new);
                    if (update_any_disps) {
                        if (update_all_disps) {
                            dispersals(k_new, k_old)++;
                        } else dispersals(x_new, y_new)++;
                    }

                    // Sample for whether aphid or plant is inoculated:
                    sample_inoculation(p_load_alate, p_load_plant, u, has_virus,
                                       infectious.at(x_new, y_new),
                                       exposed.at(x_new, y_new),
                                       exp_days.at(x_new, y_new),
                                       eng);

                    k_old = k_new;
                    x_old = x_new;
                    y_old = y_new;
                    t++;
                }

                // If reaching max fly iterations, adding alate to winged
                // population of plant settled on:
                if (t >= max_fly_t && keep_going) {
                    AphidPop& aphids_old(aphids[k_old]);
                    aphids_old.alates.X.at(adult_age + i) += 1;
                }

                n_alates_xyi--;

            }
        }


    }


    // const std::vector<XY>& alate_plants
    // uint32 x = x0;
    // uint32 y = y0;
    // // Do first iteration where alate has to leave plant:
    // first_iterate(x, y, eng);
    // // Do remaining iterations:
    // bool keep_going = false;
    // uint32 t = 1;
    // while (t < max_fly_t && keep_going) {
    //     keep_going = this->iterate(x, y, eng);
    //     t++;
    // }
    return;
}
