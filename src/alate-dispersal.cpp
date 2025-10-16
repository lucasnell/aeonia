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
    arma::imat dxdy_i;
    double radius2 = radius * radius;
    for (int32 dx = -fl_radius; dx <= fl_radius; dx++) {
        max_dy = std::floor(std::sqrt(radius2 - static_cast<double>(dx * dx)));
        rows_x = max_dy * 2U;
        if (dx != 0) rows_x++;
        dxdy_i.set_size(rows_x, 2U);
        uint32 i = 0;
        for (int32 dy = -max_dy; dy <= max_dy; dy++) {
            if (dy == 0 && dx == 0) continue;
            dxdy_i(i,0) = dx;
            dxdy_i(i,1) = dy;
            i++;
        }
        dxdy_vec.push_back(dxdy_i);
        total_rows += dxdy_i.n_rows;
    }

    arma::imat neigh_dxdy(total_rows, 2U);
    uint32 k = 0;
    for (const arma::imat& dxdy : dxdy_vec) {
        for (uint32 j = 0; j < dxdy.n_rows; j++) {
            neigh_dxdy(k,0) = dxdy(j,0);
            neigh_dxdy(k,1) = dxdy(j,1);
            k++;
        }
    }

    return neigh_dxdy;
}



/*
 Fill statuses (`virus` and `pseudo`), plus `samplers` and related objects
 (notably `land_wts` and `neighbors`).
 */
void AlateFlightInfo::fill_status_samplers(const arma::umat& landscape_,
                                           const double& radius) {


    // Fill `virus` and `pseudo` landscapes:
    virus.reserve(n_x);
    pseudo.reserve(n_x);
    bool virus_xy, pseudo_xy;
    std::vector<bool> virus_x;
    std::vector<bool> pseudo_x;
    virus_x.reserve(n_y);
    pseudo_x.reserve(n_y);
    for (uint32 x = 0; x < n_x; x++) {
        virus_x.clear();
        pseudo_x.clear();
        for (uint32 y = 0; y < n_y; y++) {
            virus_xy = get_bit_bool(0U, landscape_(x, y));
            pseudo_xy = get_bit_bool(1U, landscape_(x, y));
            virus_x.push_back(virus_xy);
            pseudo_x.push_back(pseudo_xy);
        }
        virus.push_back(virus_x);
        pseudo.push_back(pseudo_x);
    }


    // dx and dy for all neighbors within `radius`:
    arma::imat neigh_dxdy = make_neigh_dxdy(radius);
    n_neigh = neigh_dxdy.n_rows; // not true for plants on edge


    // First set `land_wts` and `neighbors`:
    land_wts.reserve(n_plants);
    neighbors.reserve(n_plants);
    uint32 ki, xi, yi;
    std::vector<uint32> neighbors_k;
    neighbors_k.reserve(n_neigh);
    for (uint32 x = 0; x < n_x; x++) {
        for (uint32 y = 0; y < n_y; y++) {
            // --------------------------
            // Set this plant's weight:
            land_wts.push_back(1.0);
            if (virus[x][y]) land_wts.back() *= virus_attract;
            if (pseudo[x][y]) land_wts.back() /= pseudo_repel;
            if (land_wts.back() == 0) {
                std::string err_msg = "\nERROR: virus_attract is very low or ";
                err_msg += "pseudo_repel is very high, ";
                err_msg += "resulting in sampling weights equal to zero. ";
                err_msg += "This results in computational problems.";
                stop(err_msg.c_str());
            }
            // --------------------------
            // Set this plant's neighbors:
            neighbors_k.clear();
            for (uint32 i = 0; i < n_neigh; i++) {
                if (i >= neigh_dxdy.n_rows) stop("Huh");
                const int32& dx(neigh_dxdy(i,0));
                const int32& dy(neigh_dxdy(i,1));
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
            neighbors.push_back(neighbors_k);
        }
    }


    // Now go back through and populate the samplers:
    samplers.reserve(n_plants);
    std::vector<double> weights_k;
    weights_k.reserve(n_neigh);
    for (uint32 k = 0; k < n_plants; k++) {
        weights_k.clear();
        for (const uint32& nk : neighbors[k]) {
            weights_k.push_back(land_wts[nk]);
        }
        samplers.push_back(AliasSampler(weights_k));
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
                                                bool& infectious,
                                                bool& exposed,
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
 */
void AlateFlightInfo::infest(const double& p_load_alate,
                             const double& p_load_plant,
                             std::vector<XY>& alate_plants,
                             std::vector<std::vector<OnePlant>>& plants,
                             arma::umat& n_alates,
                             pcg32& eng) {

    if (alate_plants.empty()) return;

    // Check for whether samplers need reconstructed and do it if needed:
    if (any_changed) {
        std::vector<double> new_wts;
        new_wts.reserve(n_neigh);
        for (uint32 k = 0; k < n_plants; k++) {
            if (update_sampler[k]) {
                new_wts.clear();
                for (const uint32& ki : neighbors[k]) {
                    new_wts.push_back(land_wts[ki]);
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
    uint32 k_old, x_old, y_old, k_new, x_new, y_new;
    bool keep_going;
    uint32 t;

    // Now go through in this random order to sample virus spread:
    for (const XY& alate_coords : alate_plants) {

        const uint32& x0(alate_coords.x);
        const uint32& y0(alate_coords.y);
        const uint32  k0 = dim_conv.to_1d(x0, y0);
        uint32& n_alates_xy(n_alates(x0,y0));
        OnePlant& plant0(plants[x0][y0]);

        while (n_alates_xy > 0) {

            k_old = k0;
            x_old = x0;
            y_old = y0;

            // If it starts on an infectious plant, then sample for whether
            // the alate is virus-bearing:
            if (plant0.infectious) {
                u = runif_01(eng);
                has_virus = u < p_load_alate;
            } else has_virus = false;

            // Sample for new location (alate always leaves first plant)
            sample(k_new, k_old, eng);
            dim_conv.to_2d(x_new, y_new, k_new);

            // Sample for whether aphid or plant is inoculated:
            sample_inoculation(p_load_alate, p_load_plant, u, has_virus,
                               plants[x_new][y_new].infectious,
                               plants[x_new][y_new].exposed,
                               plants[x_new][y_new].exp_days,
                               eng);

            k_old = k_new;
            x_old = x_new;
            y_old = y_new;

            keep_going = true;
            t = 1;

            while (t < max_fly_t && keep_going) {

                // Sample for whether alate will stay to feed at this plant:
                feed_p = w;
                if (virus[x_old][y_old]) feed_p *= epsilon;
                if (runif_01(eng) < feed_p) {
                    keep_going = false;
                    // Adding alate to winged population of plant settled on:
                    plants[x_old][y_old].aphids.alate_adults() += 1;
                    break;
                }

                // Sample for new location
                sample(k_new, k_old, eng);
                dim_conv.to_2d(x_new, y_new, k_new);

                // Sample for whether aphid or plant is inoculated:
                sample_inoculation(p_load_alate, p_load_plant, u, has_virus,
                                   plants[x_new][y_new].infectious,
                                   plants[x_new][y_new].exposed,
                                   plants[x_new][y_new].exp_days,
                                   eng);

                k_old = k_new;
                x_old = x_new;
                y_old = y_new;
                t++;
            }

            n_alates_xy--;

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
