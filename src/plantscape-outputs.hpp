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


void ps_out_none(DataFrame& out_df,
                 const std::vector<PlantScape>& plantscapes,
                 const arma::ucube& landscapes,
                 const uint32& max_t,
                 const bool& out_pseudo) {

    uint32 n_x = landscapes.n_rows;
    uint32 n_y = landscapes.n_cols;
    uint32 n_reps = landscapes.n_slices;

    uint32 n_rows = n_reps * (max_t + (uint32)1U) * n_x * n_y;

    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<double> virus;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> parasitized;
    std::vector<double> mummies;
    std::vector<double> wasps;
    std::vector<uint32> plant_x;
    std::vector<uint32> plant_y;
    std::vector<bool> pseudo;

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
        for (uint32 t = 0; t < plantscapes[r].output.size(); t++) {
            const arma::mat& rep_out(plantscapes[r].output[t]);
            for (uint32 i = 0; i < rep_out.n_rows; i++) {
                rep.push_back(r+1);
                time.push_back(t+1);
                plant_x.push_back(rep_out(i,0));
                plant_y.push_back(rep_out(i,1));
                if (out_pseudo) {
                    const uint32& lxy = landscapes(static_cast<uint32>(rep_out(i,0)-1),
                                                   static_cast<uint32>(rep_out(i,1)-1),
                                                   r);
                    pseudo.push_back(get_bit_bool(1U, lxy));
                }
                virus.push_back(rep_out(i,2));
                aphids.push_back(rep_out(i,3));
                alates.push_back(rep_out(i,4));
                parasitized.push_back(rep_out(i,5));
                mummies.push_back(rep_out(i,6));
                wasps.push_back(rep_out(i,7));
            }
        }
    }

    if (out_pseudo) {
        out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                                   _["x"] = plant_x, _["y"] = plant_y,
                                   _["pseudo"] = pseudo,
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
    uint32 n_rows = n_reps * (max_t + (uint32)1U) * (uint32)2U;

    std::vector<uint32> rep;
    std::vector<uint32> time;
    std::vector<double> virus;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> parasitized;
    std::vector<double> mummies;
    std::vector<double> wasps;
    std::vector<bool> pseudo;
    std::vector<uint32> np;

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
        for (uint32 t = 0; t < plantscapes[r].output.size(); t++) {
            const arma::mat& rep_out(plantscapes[r].output[t]);
            for (uint32 i = 0; i < rep_out.n_rows; i++) {
                rep.push_back(r+1);
                time.push_back(t+1);
                pseudo.push_back(rep_out(i,0));
                np.push_back(rep_out(i,1));
                virus.push_back(rep_out(i,2));
                aphids.push_back(rep_out(i,3));
                alates.push_back(rep_out(i,4));
                parasitized.push_back(rep_out(i,5));
                mummies.push_back(rep_out(i,6));
                wasps.push_back(rep_out(i,7));
            }
        }
    }

    out_df = DataFrame::create(_["rep"] = rep, _["time"] = time,
                               _["pseudo"] = pseudo, _["n"] = np,
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
        for (uint32 t = 0; t < plantscapes[r].output.size(); t++) {
            const arma::mat& rep_out(plantscapes[r].output[t]);
            rep.push_back(r+1);
            time.push_back(t+1);
            virus.push_back(rep_out(0,0));
            aphids.push_back(rep_out(0,1));
            alates.push_back(rep_out(0,2));
            parasitized.push_back(rep_out(0,3));
            mummies.push_back(rep_out(0,4));
            wasps.push_back(rep_out(0,5));
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
        infect_time(r) = plantscapes[r].output.size() + 1;

        for (uint32 t = 0; t < plantscapes[r].output.size(); t++) {

            const arma::mat& rep_out(plantscapes[r].output[t]);
            const double& virus_rt(rep_out(0,0));
            const double& aphids_rt(rep_out(0,1));
            const double& alates_rt(rep_out(0,2));
            const double& parasitized_rt(rep_out(0,3));
            const double& mummies_rt(rep_out(0,4));
            const double& wasps_rt(rep_out(0,5));

            p_alates(r) += alates_rt / (alates_rt + aphids_rt);
            log_aphids(r) += std::log(alates_rt + 1);
            aphids(r) += alates_rt;
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

        if (infect_time(r) > plantscapes[r].output.size()) {
            infect_time(r) = arma::datum::nan;
        }

        // Convert from sums to means:
        double n_dbl = static_cast<double>(plantscapes[r].output.size());
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

