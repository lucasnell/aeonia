

#include "plantscape-outputs.hpp"



using namespace Rcpp;



/*
 This creates a character vector for output column names.
 This does the heavy lifting for the next two functions.
 */
CharacterVector col_namer__(const std::string& summ,
                            const bool& out_pseudo,
                            const bool& out_attack_surv,
                            const bool& out_stages) {

    CharacterVector out = {"rep"};
    if (summ != "all") out.push_back("time");
    if (summ == "none") {
        out.push_back("x");
        out.push_back("y");
        if (out_pseudo) out.push_back("pseudo");
        if (out_attack_surv) out.push_back("attack_surv");
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
        out.push_back("n_infected");
    }

    return out;
}


/*
 Extract summary type based on a `ScapeSimmer` object, create column names, then
 verify that the output is the same length as `n_cols`
 */
CharacterVector col_namer_cpp(const ScapeSimmer& simmer,
                            const bool& out_pseudo,
                            const bool& out_attack_surv,
                            const bool& out_stages,
                            const uint32& n_cols) {

    const std::string& summ(simmer.summ);

    CharacterVector out = col_namer__(summ, out_pseudo, out_attack_surv, out_stages);

    if (out.size() != n_cols) {
        if (out.size() == 0) {
            Rcout << "col_names is empty !!" << std::endl;
        } else {
            Rcout << "col_names: " << out[0];
            for (int i = 1; i < out.size(); i++) {
                Rcout << ", " << out[i];
            }
            Rcout << std::endl;
        }
        Rcout << "n_cols = " << n_cols << std::endl;
        stop("INTERNAL ERROR: col_namer failure");
    }

    return out;

}


//' Just create column names, for use in R.
//'
//' @inheritParams sim_plantscape
//'
//' @export
//'
//[[Rcpp::export]]
CharacterVector col_namer(const std::string& summ,
                          const bool& out_pseudo,
                          const bool& out_attack_surv,
                          const bool& out_stages) {

    if (summ != "none" && summ != "time" && summ != "all") {
        stop("`summ` should be 'none', 'time', or 'all'");
    }

    CharacterVector out = col_namer__(summ, out_pseudo, out_attack_surv, out_stages);

    return out;

}




// Convert from list of columns to a dataframe:
void list_to_data_frame(DataFrame& out_df,
                        const std::vector<std::vector<double>>& tmp_list,
                        const CharacterVector& col_names) {

    Environment pkg = Environment::namespace_env("tibble");
    Function as_tibble = pkg["as_tibble"];
    out_df = as_tibble(tmp_list, Named(".name_repair") = "unique_quiet");
    out_df.names() = col_names; // <-- adds column names

    return;
}




std::vector<arma::span> make_spans(const bool& out_stages,
                                   const AphidPop& aphid_pop) {

    // Needed for calculating aphid by age and morphological stage:
    uint32 n_age_stages = aphid_pop.n_age_stages();
    const uint32& adult_age(aphid_pop.adult_age);
    if (aphid_pop.n_stages() != (2U * n_age_stages)) {
        stop("Alates and apterous aphids must have the same number of age stages");
    }

    std::vector<arma::span> spans;
    if (out_stages) {
        // Used if stages are output:
        spans.reserve(4);
        spans.push_back(arma::span(0U, adult_age - 1U));
        spans.push_back(arma::span(adult_age, n_age_stages - 1U));
        spans.push_back(arma::span(n_age_stages, n_age_stages + adult_age - 1U));
        spans.push_back(arma::span(n_age_stages + adult_age, n_age_stages * 2U - 1U));
    } else {
        // Used if they're not:
        spans.reserve(2);
        spans.push_back(arma::span(0U, n_age_stages - 1U));
        spans.push_back(arma::span(n_age_stages, n_age_stages * 2U - 1U));
    }

    return spans;

}



// aphids, alates OR
// aphids_juv, aphids_adu, alates_juv, alates_adu
// (see construction of `spans` above)
inline void push_aphids(const arma::vec& aphids,
                       const std::vector<arma::span>& spans,
                       std::vector<std::vector<double>>& tmp_list,
                       uint32& k,
                       arma::vec& tot_aphids) {
    for (uint32 si = 0; si < spans.size(); si++) {
        const arma::span& sp(spans[si]);
        tmp_list[k].push_back(arma::accu(aphids(sp)));
        tot_aphids(si) += tmp_list[k].back();
        k++;
    }
    return;
}
// Same as above, but without adding to `tot_aphids`
inline void push_aphids(const arma::vec& aphids,
                       const std::vector<arma::span>& spans,
                       std::vector<std::vector<double>>& tmp_list,
                       uint32& k) {
    for (const arma::span& sp : spans) {
        tmp_list[k].push_back(arma::accu(aphids(sp)));
        k++;
    }
    return;
}
// Instead of pushing to back, it adds to each. Used in `ps_out_all`.
inline void add_aphids(const arma::vec& aphids,
                       const std::vector<arma::span>& spans,
                       std::vector<std::vector<double>>& tmp_list,
                       uint32& k) {

    double ap;
    for (const arma::span& sp : spans) {
        ap = arma::accu(aphids(sp));
        tmp_list[k].back() += std::log(ap + 1);
        k++;
        tmp_list[k].back() += ap;
        k++;
    }

    return;
}










void ps_out_none(DataFrame& out_df,
                 const std::vector<ScapeSimmer>& simmers,
                 const arma::ucube& landscapes,
                 const uint32& max_t,
                 const bool& out_pseudo,
                 const bool& out_attack_surv,
                 const bool& out_stages) {

    const AphidPop& aphid_pop(simmers[0].scape.aphids[0,0]);
    const WaspPop& wasp_pop(simmers[0].scape.wasps);

    uint32 n_reps = landscapes.n_slices;
    uint32 n_rows;
    uint32 n_cols = 10;
    if (out_stages) n_cols += 2;
    if (out_pseudo) n_cols++;
    if (out_attack_surv) n_cols++;

    uint32 n_x = landscapes.n_rows;
    uint32 n_y = landscapes.n_cols;
    // Adding 1 to n_x*n_y below to account for separate row for adult
    // parasitoids (since they operate across all plants)
    n_rows = n_reps * (max_t + (uint32)1U) * (n_x * n_y + 1);

    CharacterVector col_names = col_namer_cpp(simmers[0], out_pseudo,
                                              out_attack_surv, out_stages,
                                              n_cols);

    std::vector<std::vector<double>> tmp_list(n_cols);
    for (std::vector<double>& x : tmp_list) x.reserve(n_rows);


    double tot_virus, tot_parasitized, tot_mummies;
    arma::vec tot_aphids(((out_stages) ? 4U : 2U), arma::fill::none);

    // Used if attack survivals are output:
    arma::vec X, A_surv, A_surv_apt;
    double x, A;
    const arma::vec& attack_surv(aphid_pop.attack_surv);

    // make spans for summing by desired stages:
    std::vector<arma::span> spans = make_spans(out_stages, aphid_pop);

    uint32 k;
    for (uint32 r = 0; r < n_reps; r++) {
        for (uint32 t = 0; t < simmers[r].output_dens.size(); t++) {

            const OutDensities& out_dens(simmers[r].output_dens[t]);
            const arma::umat& out_ids(simmers[r].output_ids[t]);

            tot_virus = 0;
            tot_aphids.zeros();
            tot_parasitized = 0;
            tot_mummies = 0;

            // Aphid densities:
            for (uint32 i = 0; i < out_dens.size(); i++) {

                tmp_list[0].push_back(r+1);             // rep
                tmp_list[1].push_back(t+1);             // time
                tmp_list[2].push_back(out_ids.at(i,0));    // x
                tmp_list[3].push_back(out_ids.at(i,1));    // y

                k = 4;

                if (out_pseudo) {
                    const uint32& lxy = landscapes(out_ids.at(i,0)-(uint32)1,
                                                   out_ids.at(i,1)-(uint32)1,
                                                   r);
                    tmp_list[k].push_back(get_bit_int(1U, lxy));    // pseudo
                    k++;
                }
                if (out_attack_surv) {
                    const double& Yi(out_dens.wasps[i]);
                    X = out_dens.aphids[i];
                    x = arma::accu(X);
                    wasp_pop.A_mats(A_surv_apt, attack_surv, Yi, x); // set A_surv_apt
                    // A_surv_apt is for apterous only. Doubling vector to
                    // include alates while also making adult alates not
                    // able to be parasitized:
                    if (A_surv.n_elem != (A_surv_apt.n_elem * 2U))
                        A_surv.set_size(A_surv_apt.n_elem * 2U);
                    for (uint32 i = 0; i < A_surv_apt.n_elem; i++) {
                        A_surv.at(i) = A_surv_apt.at(i);
                        if (i < aphid_pop.adult_age) {
                            A_surv.at(i+A_surv_apt.n_elem) = A_surv_apt.at(i);
                        } else A_surv.at(i+A_surv_apt.n_elem) = 1;
                    }
                    X /= x;
                    A = arma::accu(A_surv % X); // weighted mean
                    tmp_list[k].push_back(A);    // attack_surv
                    k++;
                }

                tmp_list[k].push_back(out_dens.virus[i]);
                tot_virus += out_dens.virus[i];
                k++;

                push_aphids(out_dens.aphids[i], spans, tmp_list, k, tot_aphids);

                tmp_list[k+0].push_back(out_dens.parasitized[i]);
                tmp_list[k+1].push_back(out_dens.mummies[i]);
                tmp_list[k+2].push_back(out_dens.wasps[i]);

                tot_parasitized += out_dens.parasitized[i];
                tot_mummies += out_dens.mummies[i];

            }



            // Adult parasitoids (and totals across all plants):
            tmp_list[0].push_back(r+1);         // rep
            tmp_list[1].push_back(t+1);         // time
            tmp_list[2].push_back(NA_REAL);     // x
            tmp_list[3].push_back(NA_REAL);     // y
            k = 4;
            if (out_pseudo) {
                tmp_list[k].push_back(NA_REAL); // pseudo
                k++;
            }
            if (out_attack_surv) {
                tmp_list[k].push_back(NA_REAL); // attack_surv
                k++;
            }
            tmp_list[k].push_back(tot_virus); // virus
            k++;
            for (double& a : tot_aphids) {
                tmp_list[k].push_back(a); // aphids/alates
                k++;
            }
            tmp_list[k+0].push_back(tot_parasitized);
            tmp_list[k+1].push_back(tot_mummies);
            tmp_list[k+2].push_back(out_dens.tot_wasps);

        }
    }

    list_to_data_frame(out_df, tmp_list, col_names);

    // Adjust some columns to integers:
    std::vector<std::string> int_cols;
    int_cols = {"rep", "time", "x", "y"};
    if (out_pseudo) int_cols.push_back("pseudo");
    for (std::string& s : int_cols) out_df[s] = as<IntegerVector>(out_df[s]);


    return;

}





/*
 Note that `n_actual_rows` here should be the actual number of items in the
 output dataframe, not `n_rows` from the `ps_out_time` and `ps_out_all`
 functions since those are the max possible rows.
 There can be a mismatch between these values when simulations
 stop once all plants are infected.
 */
List make_disp_col(const std::vector<ScapeSimmer>& simmers,
                   const uint32& n_actual_rows,
                   const uint32& n_plants) {

    List disp_col(n_actual_rows);
    // Make row and column names:
    CharacterVector mat_names(n_plants);
    uint32 x,y;
    for (uint32 k = 0; k < n_plants; k++) {
        simmers[0].to_2d(x, y, k);
        mat_names[k] = std::to_string(x+1) + "_" + std::to_string(y+1);
    }

    uint32 k = 0;
    for (uint32 r = 0; r < simmers.size(); r++) {
        for (const arma::umat& disps : simmers[r].dispersals) {
            if (disps.n_rows != n_plants || disps.n_cols != n_plants) {
                stop("INTERNAL ERROR: inconsistent plantscape dispersal objects");
            }
            IntegerMatrix m = wrap(disps);
            rownames(m) = mat_names;
            colnames(m) = mat_names;
            if (k >= n_actual_rows) stop("k >= n_actual_rows");
            disp_col[k] = m;
            k++;
        }
    }

    return disp_col;
}




void ps_out_time(DataFrame& out_df,
                 const std::vector<ScapeSimmer>& simmers,
                 const arma::ucube& landscapes,
                 const uint32& max_t,
                 const bool& out_stages) {

    const AphidPop& aphid_pop(simmers[0].scape.aphids[0,0]);

    uint32 n_reps = landscapes.n_slices;
    uint32 n_rows = n_reps * (max_t + (uint32)1U);

    uint32 n_cols = 8;
    if (out_stages) n_cols += 2;

    CharacterVector col_names = col_namer_cpp(simmers[0], false, false,
                                              out_stages, n_cols);

    // Check if we need to also add dispersal events:
    uint32 n_plants = landscapes.n_rows * landscapes.n_cols;
    bool out_dispersals = simmers[0].dispersals[0].n_rows == n_plants;
    if (out_dispersals) {
        col_names.push_back("disps");
        n_cols++;
    }

    std::vector<std::vector<double>> tmp_list(n_cols);
    for (std::vector<double>& x : tmp_list) x.reserve(n_rows);

    // make spans for summing by desired stages:
    std::vector<arma::span> spans = make_spans(out_stages, aphid_pop);

    uint32 k;
    for (uint32 r = 0; r < n_reps; r++) {
        for (uint32 t = 0; t < simmers[r].output_dens.size(); t++) {

            const OutDensities& out_dens(simmers[r].output_dens[t]);
            // < No IDs for this summary type! >

            tmp_list[0].push_back(r+1);             // rep
            tmp_list[1].push_back(t+1);             // time
            tmp_list[2].push_back(out_dens.virus[0]); // virus

            k = 3;

            push_aphids(out_dens.aphids[0], spans, tmp_list, k);

            tmp_list[k+0].push_back(out_dens.parasitized[0]);   // parasitized
            tmp_list[k+1].push_back(out_dens.mummies[0]);       // mummies
            tmp_list[k+2].push_back(out_dens.tot_wasps);        // wasps

            // Just filler for disps column:
            if (out_dispersals) tmp_list.back().push_back(0);

        }
    }

    list_to_data_frame(out_df, tmp_list, col_names);
    std::vector<std::string> int_cols = {"rep", "time"};
    for (std::string& s : int_cols) out_df[s] = as<IntegerVector>(out_df[s]);

    // Now add dispersal events list column:
    if (out_dispersals) {
        List disp_col = make_disp_col(simmers, tmp_list.back().size(), n_plants);
        out_df["disps"] = disp_col;
    }


    return;

}



void ps_out_all(DataFrame& out_df,
                const std::vector<ScapeSimmer>& simmers,
                const arma::ucube& landscapes,
                const uint32& max_t,
                const bool& out_stages,
                const uint32& infect_time_n,
                const double& aphid_gone_thresh,
                const double& wasp_gone_thresh) {

    const AphidPop& aphid_pop(simmers[0].scape.aphids[0,0]);

    uint32 n_reps = landscapes.n_slices;
    uint32 n_rows = n_reps;

    uint32 n_cols = 16;
    if (out_stages) n_cols += 4;

    CharacterVector col_names = col_namer_cpp(simmers[0], false, false,
                                              out_stages, n_cols);

    // Check if we need to also add dispersal events:
    uint32 n_plants = landscapes.n_rows * landscapes.n_cols;
    bool out_dispersals = simmers[0].dispersals[0].n_rows == n_plants;
    if (out_dispersals) {
        col_names.push_back("disps");
        n_cols++;
    }

    std::vector<std::vector<double>> tmp_list(n_cols);
    for (std::vector<double>& x : tmp_list) x.reserve(n_rows);


    uint32 infect_idx = n_cols - 2U;
    uint32 infected_idx = n_cols - 1U;
    if (out_dispersals) {
        infect_idx--;
        infected_idx--;
    }

    // make spans for summing by desired stages:
    std::vector<arma::span> spans = make_spans(out_stages, aphid_pop);
    // make spans for apterous vs alate (used for `p_alate` column)
    std::vector<arma::span> ava_spans = make_spans(false, aphid_pop);
    double alate_rt, apterous_rt, total_aphids_rt;

    uint32 k;

    for (uint32 r = 0; r < n_reps; r++) {

        tmp_list[0].push_back(r+1); // rep
        for (uint32 j = 1; j < n_cols; j++) tmp_list[j].push_back(0.0);

        double& infect_time(tmp_list[infect_idx].back());
        double& n_infected(tmp_list[infected_idx].back());

        infect_time = simmers[r].output_dens.size() + 1;

        for (uint32 t = 0; t < simmers[r].output_dens.size(); t++) {

            const OutDensities& out_dens(simmers[r].output_dens[t]);
            // < No IDs for this summary type! >

            const double& virus(out_dens.virus[0]);
            const arma::vec& aphids(out_dens.aphids[0]);
            const double& parasitized(out_dens.parasitized[0]);
            const double& mummies(out_dens.mummies[0]);
            const double& wasps(out_dens.tot_wasps);

            k = 1;
            apterous_rt = arma::accu(aphids(ava_spans[0]));
            alate_rt = arma::accu(aphids(ava_spans[1]));
            total_aphids_rt = apterous_rt + alate_rt;

            tmp_list[k].back() += (alate_rt / total_aphids_rt);  // p_alates
            k++;

            add_aphids(aphids, spans, tmp_list, k);

            tmp_list[k+0].back() += std::log(parasitized + 1);   // log_parasitized
            tmp_list[k+1].back() += parasitized;                 // parasitized
            tmp_list[k+2].back() += std::log(mummies + 1);       // log_mummies
            tmp_list[k+3].back() += mummies;                     // mummies
            tmp_list[k+4].back() += std::log(wasps + 1);         // log_wasps
            tmp_list[k+5].back() += wasps;                       // wasps
            if ((total_aphids_rt + parasitized) < aphid_gone_thresh) {
                tmp_list[k+6].back() += 1;
            }
            if (wasps < wasp_gone_thresh) tmp_list[k+7].back() += 1;

            if (t < infect_time && (uint32)virus >= infect_time_n) {
                infect_time = t;
            }
            if (virus > n_infected) n_infected = virus;

        }

        if (infect_time > simmers[r].output_dens.size()) {
            infect_time = NA_REAL;
        }

        // Convert from sums to means:
        double n_dbl = static_cast<double>(simmers[r].output_dens.size());
        for (uint32 j = 1; j < infect_idx; j++) tmp_list[j].back() /= n_dbl;

    }


    list_to_data_frame(out_df, tmp_list, col_names);

    std::vector<std::string> int_cols = {"rep", "n_infected"};
    for (std::string& s : int_cols) out_df[s] = as<IntegerVector>(out_df[s]);


    // Now add dispersal events list column:
    if (out_dispersals) {
        List disp_col = make_disp_col(simmers, tmp_list.back().size(), n_plants);
        out_df["disps"] = disp_col;
    }


    return;

}



