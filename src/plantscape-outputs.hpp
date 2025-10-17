#ifndef __AEONIA_PLANTSCAPE_OUTPUTS_H
#define __AEONIA_PLANTSCAPE_OUTPUTS_H


#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng


#include "aeonia_types.hpp"     // integer types
#include "plantscape.hpp"       // PlantScape class
#include "pcg.hpp"              // mt_seeds fxn
#include "util.hpp"             // thread_check fxn



using namespace Rcpp;




inline void fill_dens_tots(const OutDensities& out_dens,
                           std::vector<double>& virus,
                           std::vector<double>& aphids,
                           std::vector<double>& alates,
                           std::vector<double>& parasitized,
                           std::vector<double>& mummies,
                           std::vector<double>& wasps) {


    if (out_dens.size() > 1) {

        double tot_virus = 0;
        double tot_aphids = 0;
        double tot_alates = 0;
        double tot_parasitized = 0;
        double tot_mummies = 0;

        // Aphid densities:
        for (uint32 i = 0; i < out_dens.size(); i++) {

            virus.push_back(out_dens.virus[i]);
            aphids.push_back(out_dens.aphids[i]);
            alates.push_back(out_dens.alates[i]);
            parasitized.push_back(out_dens.parasitized[i]);
            mummies.push_back(out_dens.mummies[i]);
            wasps.push_back(NA_REAL);

            tot_virus += out_dens.virus[i];
            tot_aphids += out_dens.aphids[i];
            tot_alates += out_dens.alates[i];
            tot_parasitized += out_dens.parasitized[i];
            tot_mummies += out_dens.mummies[i];
        }

        // Adult parasitoids (and totals across all plants):
        virus.push_back(tot_virus);
        aphids.push_back(tot_aphids);
        alates.push_back(tot_alates);
        parasitized.push_back(tot_parasitized);
        mummies.push_back(tot_mummies);
        wasps.push_back(out_dens.wasps);

    } else if (out_dens.size() == 1) {

        virus.push_back(out_dens.virus[0]);
        aphids.push_back(out_dens.aphids[0]);
        alates.push_back(out_dens.alates[0]);
        parasitized.push_back(out_dens.parasitized[0]);
        mummies.push_back(out_dens.mummies[0]);
        wasps.push_back(out_dens.wasps);

    } else {

        stop("ERROR: empty `out_dens` inside `fill_dens_tots`");

    }

    return;


}







void ps_out_none(DataFrame& out_df,
                 const std::vector<PlantScape>& plantscapes,
                 const arma::ucube& landscapes,
                 const uint32& max_t,
                 const bool& out_pseudo) {

    uint32 n_x = landscapes.n_rows;
    uint32 n_y = landscapes.n_cols;
    uint32 n_reps = landscapes.n_slices;

    // Adding 1 to n_x*n_y below to account for separate row for adult
    // parasitoids (since they operate across all plants)
    uint32 n_rows = n_reps * (max_t + (uint32)1U) * (n_x * n_y + 1);

    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<uint32> plant_x;
    std::vector<uint32> plant_y;
    std::vector<int> pseudo;
    std::vector<double> virus;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> parasitized;
    std::vector<double> mummies;
    std::vector<double> wasps;

    rep.reserve(n_rows);
    time.reserve(n_rows);
    plant_x.reserve(n_rows);
    plant_y.reserve(n_rows);
    if (out_pseudo) pseudo.reserve(n_rows);
    virus.reserve(n_rows);
    aphids.reserve(n_rows);
    alates.reserve(n_rows);
    parasitized.reserve(n_rows);
    mummies.reserve(n_rows);
    wasps.reserve(n_rows);

    for (uint32 r = 0; r < n_reps; r++) {
        for (uint32 t = 0; t < plantscapes[r].output_dens.size(); t++) {

            const OutDensities& out_dens(plantscapes[r].output_dens[t]);
            const arma::umat& out_ids(plantscapes[r].output_ids[t]);

            // Fill aphids and wasps:
            fill_dens_tots(out_dens, virus, aphids, alates,
                           parasitized, mummies, wasps);

            // Fill rep, time, x, y, and optionally pseudo:
            for (uint32 i = 0; i < out_dens.size(); i++) {
                rep.push_back(r+1);
                time.push_back(t+1);
                plant_x.push_back(out_ids(i,0));
                plant_y.push_back(out_ids(i,1));
                if (out_pseudo) {
                    const uint32& lxy = landscapes(static_cast<uint32>(out_ids(i,0)-1),
                                                   static_cast<uint32>(out_ids(i,1)-1),
                                                   r);
                    pseudo.push_back(get_bit_int(1U, lxy));
                }
            }
            // Above IDs for adult parasitoids (and totals across all plants):
            rep.push_back(r+1);
            time.push_back(t+1);
            plant_x.push_back(0);
            plant_y.push_back(0);
            if (out_pseudo) pseudo.push_back(NA_LOGICAL);
        }
    }

    if (out_pseudo) {
        out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                                   _["x"] = plant_x, _["y"] = plant_y,
                                   _["pseudo"] = LogicalVector(pseudo.begin(), pseudo.end()),
                                   _["virus"] = virus, _["aphids"] = aphids,
                                   _["alates"] = alates, _["parasitized"] = parasitized,
                                   _["mummies"] = mummies, _["wasps"] = wasps);
    } else {
        out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                                   _["x"] = plant_x, _["y"] = plant_y,
                                   _["virus"] = virus, _["aphids"] = aphids,
                                   _["alates"] = alates, _["parasitized"] = parasitized,
                                   _["mummies"] = mummies, _["wasps"] = wasps);
    }
    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return;

}






void ps_out_pseudo(DataFrame& out_df,
                   const std::vector<PlantScape>& plantscapes,
                   const arma::ucube& landscapes,
                   const uint32& max_t) {

    uint32 n_reps = landscapes.n_slices;
    // Using 3 instead of 2 below because I need an extra row for adult
    // parasitoids that operate across all plants
    uint32 n_rows = n_reps * (max_t + (uint32)1U) * (uint32)3U;

    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<int> pseudo;
    std::vector<uint32> np;
    std::vector<double> virus;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> parasitized;
    std::vector<double> mummies;
    std::vector<double> wasps;

    rep.reserve(n_rows);
    time.reserve(n_rows);
    pseudo.reserve(n_rows);
    np.reserve(n_rows);
    virus.reserve(n_rows);
    aphids.reserve(n_rows);
    alates.reserve(n_rows);
    parasitized.reserve(n_rows);
    mummies.reserve(n_rows);
    wasps.reserve(n_rows);

    for (uint32 r = 0; r < n_reps; r++) {
        for (uint32 t = 0; t < plantscapes[r].output_dens.size(); t++) {

            const OutDensities& out_dens(plantscapes[r].output_dens[t]);
            const arma::umat& out_ids(plantscapes[r].output_ids[t]);

            // Fill aphids and wasps:
            fill_dens_tots(out_dens, virus, aphids, alates,
                           parasitized, mummies, wasps);

            // Fill rep, time, pseudo, and # plants:
            for (uint32 i = 0; i < out_dens.size(); i++) {
                rep.push_back(r+1);
                time.push_back(t+1);
                pseudo.push_back(out_ids(i,0));
                np.push_back(out_ids(i,1));
            }

            // Above IDs for adult parasitoids (and totals across all plants):
            rep.push_back(r+1);
            time.push_back(t+1);
            pseudo.push_back(NA_LOGICAL);
            np.push_back(arma::accu(out_ids.col(1)));

        }
    }

    out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                               _["pseudo"] = LogicalVector(pseudo.begin(), pseudo.end()),
                               _["n"] = np,
                               _["virus"] = virus, _["aphids"] = aphids,
                               _["alates"] = alates, _["parasitized"] = parasitized,
                               _["mummies"] = mummies, _["wasps"] = wasps);
    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return;

}




void ps_out_time(DataFrame& out_df,
                 const std::vector<PlantScape>& plantscapes,
                 const arma::ucube& landscapes,
                 const uint32& max_t) {

    uint32 n_reps = landscapes.n_slices;
    uint32 n_rows = n_reps * (max_t + (uint32)1U);

    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<double> virus;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> parasitized;
    std::vector<double> mummies;
    std::vector<double> wasps;

    rep.reserve(n_rows);
    time.reserve(n_rows);
    virus.reserve(n_rows);
    aphids.reserve(n_rows);
    alates.reserve(n_rows);
    parasitized.reserve(n_rows);
    mummies.reserve(n_rows);
    wasps.reserve(n_rows);

    for (uint32 r = 0; r < n_reps; r++) {
        for (uint32 t = 0; t < plantscapes[r].output_dens.size(); t++) {

            const OutDensities& out_dens(plantscapes[r].output_dens[t]);
            // < No IDs for this summary type! >

            // Fill aphids and wasps:
            fill_dens_tots(out_dens, virus, aphids, alates,
                           parasitized, mummies, wasps);
            // Fill rep and time:
            rep.push_back(r+1);
            time.push_back(t+1);

        }
    }

    out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                               _["virus"] = virus, _["aphids"] = aphids,
                               _["alates"] = alates, _["parasitized"] = parasitized,
                               _["mummies"] = mummies, _["wasps"] = wasps);
    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return;

}



void ps_out_all(DataFrame& out_df,
                const std::vector<PlantScape>& plantscapes,
                const arma::ucube& landscapes,
                const uint32& max_t,
                const uint32& infect_time_n) {

    uint32 n_reps = landscapes.n_slices;
    uint32 n_rows = n_reps;

    arma::uvec rep(n_rows, arma::fill::none);
    arma::vec p_alates(n_rows, arma::fill::zeros);
    arma::vec log_aphids(n_rows, arma::fill::zeros);
    arma::vec aphids(n_rows, arma::fill::zeros);
    arma::vec log_alates(n_rows, arma::fill::zeros);
    arma::vec alates(n_rows, arma::fill::zeros);
    arma::vec log_parasitized(n_rows, arma::fill::zeros);
    arma::vec parasitized(n_rows, arma::fill::zeros);
    arma::vec log_mummies(n_rows, arma::fill::zeros);
    arma::vec mummies(n_rows, arma::fill::zeros);
    arma::vec log_wasps(n_rows, arma::fill::zeros);
    arma::vec wasps(n_rows, arma::fill::zeros);
    arma::vec infect_time(n_rows, arma::fill::none);
    arma::uvec outbreak_size(n_rows, arma::fill::zeros);

    for (uint32 r = 0; r < n_reps; r++) {

        rep(r) = r+1;
        infect_time(r) = plantscapes[r].output_dens.size() + 1;

        for (uint32 t = 0; t < plantscapes[r].output_dens.size(); t++) {

            const OutDensities& out_dens(plantscapes[r].output_dens[t]);
            // < No IDs for this summary type! >

            const double& virus_rt(out_dens.virus[0]);
            const double& aphids_rt(out_dens.aphids[0]);
            const double& alates_rt(out_dens.alates[0]);
            const double& parasitized_rt(out_dens.parasitized[0]);
            const double& mummies_rt(out_dens.mummies[0]);
            const double& wasps_rt(out_dens.wasps);

            p_alates(r) += alates_rt / (alates_rt + aphids_rt);
            log_aphids(r) += std::log(aphids_rt + 1);
            aphids(r) += aphids_rt;
            log_alates(r) += std::log(alates_rt + 1);
            alates(r) += alates_rt;
            log_parasitized(r) += std::log(parasitized_rt + 1);
            parasitized(r) += parasitized_rt;
            log_mummies(r) += std::log(mummies_rt + 1);
            mummies(r) += mummies_rt;
            log_wasps(r) += std::log(wasps_rt + 1);
            wasps(r) += wasps_rt;
            if (t < infect_time(r) && (uint32)virus_rt >= infect_time_n) {
                infect_time(r) = t;
            }
            if (virus_rt > outbreak_size(r)) outbreak_size(r) = virus_rt;

        }

        if (infect_time(r) > plantscapes[r].output_dens.size()) {
            infect_time(r) = arma::datum::nan;
        }

        // Convert from sums to means:
        double n_dbl = static_cast<double>(plantscapes[r].output_dens.size());
        p_alates(r) /= n_dbl;
        log_aphids(r) /= n_dbl;
        aphids(r) /= n_dbl;
        log_alates(r) /= n_dbl;
        alates(r) /= n_dbl;
        log_parasitized(r) /= n_dbl;
        parasitized(r) /= n_dbl;
        log_mummies(r) /= n_dbl;
        mummies(r) /= n_dbl;
        log_wasps(r) /= n_dbl;
        wasps(r) /= n_dbl;
    }

    out_df = DataFrame::create(_["rep"] = rep,
                               _["p_alates"] = p_alates,
                               _["log_aphids"] = log_aphids,
                               _["aphids"] = aphids,
                               _["log_alates"] = log_alates,
                               _["alates"] = alates,
                               _["log_parasitized"] = log_parasitized,
                               _["parasitized"] = parasitized,
                               _["log_mummies"] = log_mummies,
                               _["mummies"] = mummies,
                               _["log_wasps"] = log_wasps,
                               _["wasps"] = wasps,
                               _["infect_time"] = infect_time,
                               _["outbreak_size"] = outbreak_size);
    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return;

}








#endif

