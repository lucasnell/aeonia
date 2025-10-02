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

    std::vector<uint32> rep;
    std::vector<double> p_alates;
    std::vector<double> log_aphids;
    std::vector<double> aphids;
    std::vector<double> log_alates;
    std::vector<double> alates;
    std::vector<double> log_parasitized;
    std::vector<double> parasitized;
    std::vector<double> log_mummies;
    std::vector<double> mummies;
    std::vector<double> log_wasps;
    std::vector<double> wasps;
    std::vector<uint32> infect_time;
    std::vector<uint32> outbreak_size;

    rep.reserve(n_rows);
    p_alates.reserve(n_rows);
    log_aphids.reserve(n_rows);
    aphids.reserve(n_rows);
    log_alates.reserve(n_rows);
    alates.reserve(n_rows);
    log_parasitized.reserve(n_rows);
    parasitized.reserve(n_rows);
    log_mummies.reserve(n_rows);
    mummies.reserve(n_rows);
    log_wasps.reserve(n_rows);
    wasps.reserve(n_rows);
    infect_time.reserve(n_rows);
    outbreak_size.reserve(n_rows);

    for (uint32 r = 0; r < n_reps; r++) {
        rep.push_back(r+1);
        p_alates.push_back(0);
        log_aphids.push_back(0);
        aphids.push_back(0);
        log_alates.push_back(0);
        alates.push_back(0);
        log_parasitized.push_back(0);
        parasitized.push_back(0);
        log_mummies.push_back(0);
        mummies.push_back(0);
        log_wasps.push_back(0);
        wasps.push_back(0);
        infect_time.push_back(plantscapes[r].output.size() + 1);
        outbreak_size.push_back(0);


        for (uint32 t = 0; t < plantscapes[r].output.size(); t++) {

            const arma::mat& rep_out(plantscapes[r].output[t]);
            const double& virus_rt(rep_out(0,0));
            const double& aphids_rt(rep_out(0,1));
            const double& alates_rt(rep_out(0,2));
            const double& parasitized_rt(rep_out(0,3));
            const double& mummies_rt(rep_out(0,4));
            const double& wasps_rt(rep_out(0,5));

            p_alates.back() += alates_rt / (alates_rt + aphids_rt);
            log_aphids.back() += std::log(alates_rt + 1);
            aphids.back() += alates_rt;
            log_alates.back() += std::log(alates_rt + 1);
            alates.back() += alates_rt;
            log_parasitized.back() += std::log(parasitized_rt + 1);
            parasitized.back() += parasitized_rt;
            log_mummies.back() += std::log(mummies_rt + 1);
            mummies.back() += mummies_rt;
            log_wasps.back() += std::log(wasps_rt + 1);
            wasps.back() += wasps_rt;
            if (t < infect_time.back() && (uint32)virus_rt >= infect_time_n) {
                infect_time.back() = t;
            }
            if (virus_rt > outbreak_size.back()) outbreak_size.back() = virus_rt;

        }

        // Convert from sums to means:
        double n_dbl = static_cast<double>(plantscapes[r].output.size());
        p_alates.back() /= n_dbl;
        log_aphids.back() /= n_dbl;
        aphids.back() /= n_dbl;
        log_alates.back() /= n_dbl;
        alates.back() /= n_dbl;
        log_parasitized.back() /= n_dbl;
        parasitized.back() /= n_dbl;
        log_mummies.back() /= n_dbl;
        mummies.back() /= n_dbl;
        log_wasps.back() /= n_dbl;
        wasps.back() /= n_dbl;
    }

    out_df = DataFrame::create(_["rep"] = rep, _["p_alates"] = p_alates,
                               _["log_aphids"] = log_aphids, _["aphids"] = aphids,
                               _["log_alates"] = log_alates, _["alates"] = alates,
                               _["log_parasitized"] = log_parasitized,
                               _["parasitized"] = parasitized,
                               _["log_mummies"] = log_mummies, _["mummies"] = mummies,
                               _["log_wasps"] = log_wasps, _["wasps"] = wasps,
                               _["infect_time"] = infect_time,
                               _["outbreak_size"] = outbreak_size);
    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return;

}








#endif

