#ifndef __AEONIA_PLANTSCAPE_OBSERVER_H
#define __AEONIA_PLANTSCAPE_OBSERVER_H

#include "aeonia_types.hpp"     // integer types


#include <RcppArmadillo.h>
#include <RcppThread.h>         // multithreading

#include "plantscape.hpp"       // PlantScape class
#include "convert-dims.hpp"     // DimensionConverter class



// Class to store densities by whatever grouping is desired for output.
// These values are shared no matter what type of output summary is used.
struct OutDensities {

    std::vector<double> virus;          // virus density (0 or 1 when this summarizes 1 patch)
    std::vector<arma::vec> aphids;      // aphids (all stages) density
    std::vector<double> parasitized;    // parasitized aphid density
    std::vector<double> mummies;        // mummy density
    std::vector<double> wasps;          // adult, female parasitoid wasp density
    double tot_wasps;


    // Don't include any values, later reserve using that method.
    OutDensities() : virus(), aphids(), parasitized(), mummies(), wasps(), tot_wasps(0) {};

    // Fill zeros for vectors of length `n`.
    OutDensities(const uint32& n, const uint32& n_stages)
        : virus(n, 0.0),
          aphids(n, arma::vec(n_stages, arma::fill::zeros)),
          parasitized(n, 0.0),
          mummies(n, 0.0),
          wasps(n, 0.0),
          tot_wasps(0) {};


    void reserve(const uint32& n) {
        virus.reserve(n);
        aphids.reserve(n);
        parasitized.reserve(n);
        mummies.reserve(n);
        wasps.reserve(n);
        return;
    }

    // Push back all values:
    void push_back(const double& virus_,
                   const arma::vec& aphids_,
                   const double& parasitized_,
                   const double& mummies_,
                   const double& wasps_) {
        virus.push_back(virus_);
        aphids.push_back(aphids_);
        parasitized.push_back(parasitized_);
        mummies.push_back(mummies_);
        wasps.push_back(wasps_);
        // tot_wasps += wasps_;
        return;
    }


    uint32 size() const {
        return virus.size();
    }

    void add_to(const uint32& k,
                const double& virus_,
                const arma::vec& aphids_,
                const double& parasitized_,
                const double& mummies_,
                const double& wasps_) {
        virus[k] += virus_;
        aphids[k] += aphids_;
        parasitized[k] += parasitized_;
        mummies[k] += mummies_;
        wasps[k] += wasps_;
        // tot_wasps += wasps_;
    }

};









struct ScapeSimmer {

    // Plantscape object to track:
    PlantScape scape;

    // Output objects:
    std::vector<OutDensities> output_dens;  // densities of organism types
    std::vector<arma::umat> output_ids;     // identifiers for each set of densities
    std::vector<arma::umat> dispersals;     // dispersal events

    // How to summarize output (if at all):
    std::string summ;
    // Max time points to simulate:
    uint32 max_t;
    // landscape x and y dimensions:
    uint32 n_x;
    uint32 n_y;

    // Whether to output dispersals:
    bool out_dispersals;
    // iterator for current dispersal matrix:
    std::vector<arma::umat>::iterator disp_iter;

    // Write the current state of the PlantScape to the `output` field.
    void fill_output();



    inline ScapeSimmer(const arma::umat& landscape_,
                       const DiseaseDispersal& disp_dis,
                       const InsectPops& insects,
                       const arma::mat& N0,
                       const arma::mat& W0,
                       const arma::mat& M0,
                       const double& Y0,
                       const arma::mat& wasp_attract_,
                       const std::vector<uint64>& seeds,
                       const std::string& summ_,
                       const uint32& max_t_,
                       const bool& out_dispersals_)
        : scape(landscape_, disp_dis, insects, N0, W0, M0, Y0, wasp_attract_, seeds),
          output_dens(),
          output_ids(),
          dispersals(),
          summ(summ_),
          max_t(max_t_),
          n_x(landscape_.n_rows),
          n_y(landscape_.n_cols),
          out_dispersals(out_dispersals_),
          disp_iter(),
          dim_conv(n_x, n_y) {

        // Reserve max memory required:
        output_dens.reserve(max_t+1U);
        if (summ == "none") output_ids.reserve(max_t+1U);
        // fill starting conditions:
        fill_output();

        // Optionally reserve `dispersals`:
        if (summ != "time" && summ != "all") out_dispersals = false;
        if (out_dispersals) {
            if (summ == "time") {
                dispersals.reserve(max_t + 1);
                for (uint32 i = 0; i < (max_t + 1); i++) {
                    dispersals.emplace_back(n_x*n_y, n_x*n_y, arma::fill::zeros);
                }
            } else {
                dispersals = std::vector<arma::umat>(1, arma::umat(
                    n_x * n_y, n_x * n_y, arma::fill::zeros));
            }
        } else dispersals.push_back(arma::umat());
        disp_iter = dispersals.begin();
        if (out_dispersals && summ == "time") disp_iter++;
    }


    // Run all time steps
    void run(const bool& infect_stop,
             RcppThread::ProgressBar& prog_bar,
             const bool& show_progress);

    // Used in make_disp_col
    inline void to_2d(uint32& x, uint32& y, const uint32& k) const {
        dim_conv.to_2d(x, y, k);
        return;
    }


private:

    DimensionConverter dim_conv;


};





#endif

