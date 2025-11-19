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



/*
 This creates a character vector for output column names.
 */
CharacterVector col_namer(const PlantScape& plantscape,
                          const bool& out_pseudo,
                          const bool& out_stages) {

    const std::string summ = plantscape.summary();

    CharacterVector out = {"rep"};
    if (summ != "all") out.push_back("time");
    if (summ == "none") {
        out.push_back("x");
        out.push_back("y");
        if (out_pseudo) out.push_back("pseudo");
    }
    if (summ == "pseudo") {
        out.push_back("pseudo");
        out.push_back("n");
    }
    if (summ != "all") {
        out.push_back("virus");
        if (out_stages) {
            out.push_back("aphids_juv");
            out.push_back("aphids_adu");
            out.push_back("alates_juv");
            out.push_back("alates_adu");
        } else {
            out.push_back("aphids");
            out.push_back("alates");
        }
        out.push_back("parasitized");
        out.push_back("mummies");
        out.push_back("wasps");
    } else {
        out.push_back("p_alates");
        if (out_stages) {
            out.push_back("log_aphids_juv");
            out.push_back("log_aphids_adu");
            out.push_back("aphids_juv");
            out.push_back("aphids_adu");
            out.push_back("log_alates_juv");
            out.push_back("log_alates_adu");
            out.push_back("alates_juv");
            out.push_back("alates_adu");
        } else {
            out.push_back("log_aphids");
            out.push_back("aphids");
            out.push_back("log_alates");
            out.push_back("alates");
        }
        out.push_back("log_parasitized");
        out.push_back("parasitized");
        out.push_back("log_mummies");
        out.push_back("mummies");
        out.push_back("log_wasps");
        out.push_back("wasps");
        out.push_back("aphid_gone_n");
        out.push_back("wasp_gone_n");
        out.push_back("infect_time");
        out.push_back("outbreak_size");
    }

    return out;
}


// Convert from list of rows to a dataframe:
void list_to_data_frame(DataFrame& out_df,
                        const List& tmp_list,
                        const CharacterVector& col_names) {

    // Convert to matrix, then data.frame, then tibble (I wish I knew a
    // simpler way!):
    // (NOTE: Do not try to convert directly from matrix to tibble.
    //  It turns many of the values to zeros!)
    Function do_call("do.call");
    Function rbind("rbind");
    Function as_data_frame("as.data.frame");
    Function as_tibble("as_tibble");
    NumericMatrix tmp_matrix = do_call(rbind, tmp_list);
    DataFrame tmp_dataframe = as_data_frame(tmp_matrix);
    tmp_dataframe.names() = col_names; // <-- adds column names
    out_df = as_tibble(tmp_dataframe);

    // If the list was transposed (every item is a column), you could simply
    // do this:
    // (note: in this case, you could use `std::vector<std::vector<double>>`
    //  instead of `List`)
    // Function as_tibble("as_tibble");
    // out_df = as_tibble(tmp_list, Named(".name_repair") = "unique_quiet");
    // out_df.names() = col_names; // <-- adds column names

    return;
}




void ps_out_none_pseudo(DataFrame& out_df,
                        const std::vector<PlantScape>& plantscapes,
                        const arma::ucube& landscapes,
                        const uint32& max_t,
                        const bool& out_pseudo,
                        const bool& out_stages) {

    bool summ_none = plantscapes[0].summary() == "none";

    uint32 n_reps = landscapes.n_slices;
    uint32 n_rows;
    uint32 n_cols = 10;
    if (out_stages) n_cols += 2;

    if (summ_none) {
        uint32 n_x = landscapes.n_rows;
        uint32 n_y = landscapes.n_cols;
        // Adding 1 to n_x*n_y below to account for separate row for adult
        // parasitoids (since they operate across all plants)
        n_rows = n_reps * (max_t + (uint32)1U) * (n_x * n_y + 1);
        if (out_pseudo) n_cols++;
    } else {
        // Using 3 instead of 2 below because I need an extra row for adult
        // parasitoids that operate across all plants
        n_rows = n_reps * (max_t + (uint32)1U) * (uint32)3U;
    }

    CharacterVector col_names = col_namer(plantscapes[0], out_pseudo, out_stages);
    if (col_names.size() != n_cols) {
        for (int i = 0; i < col_names.size(); i++) {
            Rcout << col_names[i] << ", ";
        }
        Rcout << std::endl;
        Rcout << "n_cols = " << std::to_string(n_cols) << std::endl;
        stop("INTERNAL ERROR: col_namer failure");
    }


    List tmp_list(n_rows);
    std::vector<double> tmp_row;
    tmp_row.reserve(n_cols);

    double tot_virus, tot_parasitized, tot_mummies;
    arma::vec tot_aphids(((out_stages) ? 4U : 2U), arma::fill::none);

    uint32 k = 0;
    for (uint32 r = 0; r < n_reps; r++) {
        for (uint32 t = 0; t < plantscapes[r].output_dens.size(); t++) {

            const OutDensities& out_dens(plantscapes[r].output_dens[t]);
            const arma::umat& out_ids(plantscapes[r].output_ids[t]);

            tot_virus = 0;
            tot_aphids.zeros();
            tot_parasitized = 0;
            tot_mummies = 0;

            // Aphid densities:
            for (uint32 i = 0; i < out_dens.size(); i++) {

                tmp_row.clear();

                tmp_row.push_back(r+1);             // rep
                tmp_row.push_back(t+1);             // time
                tmp_row.push_back(out_ids(i,0));    // x / pseudo
                tmp_row.push_back(out_ids(i,1));    // y / n

                if (out_pseudo) {
                    const uint32& lxy = landscapes(static_cast<uint32>(out_ids(i,0)-1),
                                                   static_cast<uint32>(out_ids(i,1)-1),
                                                   r);
                    tmp_row.push_back(get_bit_int(1U, lxy));    // pseudo
                }

                tmp_row.push_back(out_dens.virus[i]);
                tot_virus += out_dens.virus[i];

                const std::array<double,4>& aphids(out_dens.aphids[i]);
                if (out_stages) {
                    // aphids_juv, aphids_adu, alates_juv, alates_adu
                    for (uint32 j = 0; j < 4; j++) {
                        tmp_row.push_back(aphids[j]);
                        tot_aphids(j) += aphids[j];
                    }
                } else {
                    tmp_row.push_back(aphids[0] + aphids[1]);  // aphids
                    tot_aphids(0) += tmp_row.back();
                    tmp_row.push_back(aphids[2] + aphids[3]);  // alates
                    tot_aphids(1) += tmp_row.back();
                }

                tmp_row.push_back(out_dens.parasitized[i]);
                tmp_row.push_back(out_dens.mummies[i]);
                tmp_row.push_back(NA_REAL);

                tot_parasitized += out_dens.parasitized[i];
                tot_mummies += out_dens.mummies[i];

                tmp_list[k] = tmp_row;
                k++;

            }

            // Adult parasitoids (and totals across all plants):
            tmp_row.clear();
            tmp_row.push_back(r+1);         // rep
            tmp_row.push_back(t+1);         // time
            tmp_row.push_back(NA_REAL);     // x / pseudo
            if (summ_none) {
                tmp_row.push_back(NA_REAL);     // y
            } else tmp_row.push_back(arma::accu(out_ids.col(1)));     // n
            if (out_pseudo) tmp_row.push_back(NA_REAL); // pseudo
            tmp_row.push_back(tot_virus); // virus
            for (double& a : tot_aphids) tmp_row.push_back(a); // aphids/alates
            tmp_row.push_back(tot_parasitized);
            tmp_row.push_back(tot_mummies);
            tmp_row.push_back(out_dens.wasps);

            tmp_list[k] = tmp_row;
            k++;

        }
    }

    list_to_data_frame(out_df, tmp_list, col_names);

    // Adjust some columns to integers:
    std::vector<std::string> int_cols;
    if (summ_none) {
        int_cols = {"rep", "time", "x", "y"};
        if (out_pseudo) int_cols.push_back("pseudo");
    } else {
        int_cols = {"rep", "time", "n"};
        // pseudo is logical:
        out_df["pseudo"] = as<LogicalVector>(out_df["pseudo"]);
    }
    for (std::string& s : int_cols) out_df[s] = as<IntegerVector>(out_df[s]);


    return;

}








void ps_out_time(DataFrame& out_df,
                 const std::vector<PlantScape>& plantscapes,
                 const arma::ucube& landscapes,
                 const uint32& max_t,
                 const bool& out_stages) {

    uint32 n_reps = landscapes.n_slices;
    uint32 n_rows = n_reps * (max_t + (uint32)1U);

    uint32 n_cols = 8;
    if (out_stages) n_cols += 2;

    CharacterVector col_names = col_namer(plantscapes[0], false, out_stages);
    if (col_names.size() != n_cols) {
        for (int i = 0; i < col_names.size(); i++) {
            Rcout << col_names[i] << ", ";
        }
        Rcout << std::endl;
        Rcout << "n_cols = " << std::to_string(n_cols) << std::endl;
        stop("INTERNAL ERROR: col_namer failure");
    }

    List tmp_list(n_rows);
    std::vector<double> tmp_row;
    tmp_row.reserve(n_cols);

    uint32 k = 0;
    for (uint32 r = 0; r < n_reps; r++) {
        for (uint32 t = 0; t < plantscapes[r].output_dens.size(); t++) {

            const OutDensities& out_dens(plantscapes[r].output_dens[t]);
            // < No IDs for this summary type! >

            tmp_row.clear();

            tmp_row.push_back(r+1);             // rep
            tmp_row.push_back(t+1);             // time

            tmp_row.push_back(out_dens.virus[0]); // virus
            const std::array<double,4>& aphids(out_dens.aphids[0]);
            if (out_stages) {
                // aphids_juv, aphids_adu, alates_juv, alates_adu
                for (uint32 j = 0; j < 4; j++) {
                    tmp_row.push_back(aphids[j]);
                }
            } else {
                tmp_row.push_back(aphids[0] + aphids[1]);  // aphids
                tmp_row.push_back(aphids[2] + aphids[3]);  // alates
            }
            tmp_row.push_back(out_dens.parasitized[0]);     // parasitized
            tmp_row.push_back(out_dens.mummies[0]);         // mummies
            tmp_row.push_back(out_dens.wasps);              // wasps

            tmp_list[k] = tmp_row;
            k++;

        }
    }

    list_to_data_frame(out_df, tmp_list, col_names);
    std::vector<std::string> int_cols = {"rep", "time"};
    for (std::string& s : int_cols) out_df[s] = as<IntegerVector>(out_df[s]);


    return;

}



void ps_out_all(DataFrame& out_df,
                const std::vector<PlantScape>& plantscapes,
                const arma::ucube& landscapes,
                const uint32& max_t,
                const bool& out_stages,
                const uint32& infect_time_n,
                const double& aphid_gone_thresh,
                const double& wasp_gone_thresh) {

    uint32 n_reps = landscapes.n_slices;
    uint32 n_rows = n_reps;

    uint32 n_cols = 16;
    if (out_stages) n_cols += 4;

    CharacterVector col_names = col_namer(plantscapes[0], false, out_stages);
    if (col_names.size() != n_cols) {
        for (int i = 0; i < col_names.size(); i++) {
            Rcout << col_names[i] << ", ";
        }
        Rcout << std::endl;
        Rcout << "n_cols = " << std::to_string(n_cols) << std::endl;
        stop("INTERNAL ERROR: col_namer failure");
    }

    List tmp_list(n_rows);
    std::vector<double> tmp_row;
    tmp_row.reserve(n_cols);

    uint32 infect_idx = n_cols - 2U;
    uint32 outbreak_idx = n_cols - 1U;

    for (uint32 r = 0; r < n_reps; r++) {

        tmp_row.clear();
        tmp_row.push_back(r+1); // rep
        for (uint32 j = 1; j < n_cols; j++) tmp_row.push_back(0.0);

        double& infect_time(tmp_row[infect_idx]);
        double& outbreak_size(tmp_row[outbreak_idx]);

        infect_time = plantscapes[r].output_dens.size() + 1;

        for (uint32 t = 0; t < plantscapes[r].output_dens.size(); t++) {

            const OutDensities& out_dens(plantscapes[r].output_dens[t]);
            // < No IDs for this summary type! >

            const double& virus(out_dens.virus[0]);
            const std::array<double,4>& aphids(out_dens.aphids[0]);
            const double& parasitized(out_dens.parasitized[0]);
            const double& mummies(out_dens.mummies[0]);
            const double& wasps(out_dens.wasps);

            uint32 j = 1;
            double total_aphids = std::accumulate(aphids.begin(), aphids.end(), 0);
            tmp_row[j] += (aphids[2] + aphids[3]) / (total_aphids);  // p_alates
            j++;
            if (out_stages) {
                // aphids_juv, aphids_adu, alates_juv, alates_adu
                uint32 l2;
                for (uint32 l = 0; l < 4; l++) {
                    l2 = j + 2U * l;
                    tmp_row[l2] += std::log(aphids[l] + 1);
                    tmp_row[l2+1U] += aphids[l];
                }
                j += 8;
            } else {
                tmp_row[j+0] += std::log(aphids[0] + aphids[1] + 1);  // log_aphids
                tmp_row[j+1] += (aphids[0] + aphids[1]);  // aphids
                tmp_row[j+2] += std::log(aphids[2] + aphids[3] + 1);  // log_alates
                tmp_row[j+3] += (aphids[2] + aphids[3]);  // alates
                j += 4;
            }

            tmp_row[j+0] += std::log(parasitized + 1);   // log_parasitized
            tmp_row[j+1] += parasitized;                 // parasitized
            tmp_row[j+2] += std::log(mummies + 1);       // log_mummies
            tmp_row[j+3] += mummies;                     // mummies
            tmp_row[j+4] += std::log(wasps + 1);         // log_wasps
            tmp_row[j+5] += wasps;                       // wasps
            if ((total_aphids + parasitized) < aphid_gone_thresh) {
                tmp_row[j+6] += 1;
            }
            if (wasps < wasp_gone_thresh) tmp_row[j+7] += 1;

            if (t < infect_time && (uint32)virus >= infect_time_n) {
                infect_time = t;
            }
            if (virus > outbreak_size) outbreak_size = virus;

        }

        if (infect_time > plantscapes[r].output_dens.size()) {
            infect_time = NA_REAL;
        }

        // Convert from sums to means:
        double n_dbl = static_cast<double>(plantscapes[r].output_dens.size());
        for (uint32 j = 1; j < infect_idx; j++) tmp_row[j] /= n_dbl;

        tmp_list[r] = tmp_row;

    }



    list_to_data_frame(out_df, tmp_list, col_names);

    std::vector<std::string> int_cols = {"rep", "outbreak_size"};
    for (std::string& s : int_cols) out_df[s] = as<IntegerVector>(out_df[s]);


    return;

}








#endif
