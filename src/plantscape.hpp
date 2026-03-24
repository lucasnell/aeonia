#ifndef __AEONIA_PLANTSCAPE_H
#define __AEONIA_PLANTSCAPE_H

/*
 This contains code for landscapes of plants, including
 plant disease dynamics and *Pseudomonas* presence.
 */
#include "aeonia_types.hpp"     // integer types


#include <RcppArmadillo.h>
#include <vector>
#include <mdspan>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "convert-dims.hpp"     // XY and get_bit_bool
#include "aphids.hpp"           // AphidPop class
#include "wasps.hpp"            // WaspPop, MummyPop classes
#include "insects.hpp"          // InsectPops class
#include "alate-dispersal.hpp"  // AlateFlightInfo class
#include "pcg.hpp"              // runif_01 fxn






using namespace Rcpp;

struct ScapeSimmer; // required to declare friendship



// Simple class to store many dispersal and disease inputs to PlantScape
struct DiseaseDispersal {

    double radius;
    uint32 max_fly_t;
    double virus_attract;
    double pseudo_repel;
    double epsilon;
    double w;
    double p_load_alate;
    double p_load_plant;
    uint32 total_exp_days;


    DiseaseDispersal(const double& radius_,
                     const uint32& max_fly_t_,
                     const double& virus_attract_,
                     const double& pseudo_repel_,
                     const double& epsilon_,
                     const double& w_,
                     const double& p_load_alate_,
                     const double& p_load_plant_,
                     const uint32& total_exp_days_)
        : radius(radius_),
          max_fly_t(max_fly_t_),
          virus_attract(virus_attract_),
          pseudo_repel(pseudo_repel_),
          epsilon(epsilon_),
          w(w_),
          p_load_alate(p_load_alate_),
          p_load_plant(p_load_plant_),
          total_exp_days(total_exp_days_) {}


};














class PlantScape {

    friend struct ScapeSimmer;

protected:

    uint32 n_x;
    uint32 n_y;
    uint32 n_plants;

    DimensionConverter dim_conv;

    AlateFlightInfo flight;

    // Attractiveness to parasitoids:
    arma::mat wasp_attract;

    // Probability that an uninoculated alate is loaded with a virus if it
    // probes an infectious plant:
    double p_load_alate;
    // Probability that an uninfected plant is loaded with a virus if it
    // is probed by a virus-bearing aphid:
    double p_load_plant;

    // days since exposure required to transition to infected
    uint32 total_exp_days;

    // for storing numbers of alates per plant (rows = stages, cols/slices = x/y):
    arma::ucube n_alates;
    // for storing x,y indices of plants that have produced >=1 alate:
    std::vector<std::array<uint32, 2U>> alate_plants;


    pcg32 eng;

    // relative aphid abundance and adult, female parasitoid densities by patch
    arma::mat z_mat;
    arma::mat Yi_mat;


    // Fill adult, female parasitoid densities by patch
    // (stored in `Yi_mat`):
    void fill_Yi_mat();

    // Store raw data for insect pop mdspans below:
    std::vector<AphidPop> aphids_;
    std::vector<MummyPop> mummies_;

public:

    // Insect populations:
    WaspPop wasps;
    Span2D<AphidPop> aphids;
    Span2D<MummyPop> mummies;

    // Plant infection info:
    arma::Mat<uint16> exposed;    // virus: exposed but not yet infectious?
    arma::Mat<uint16> infectious; // virus: infectious?
    arma::Mat<uint16> pseudo;     // contains Pseudomonas?
    arma::umat exp_days;   // days since exposure (ignored if not exposed)



    PlantScape(const arma::umat& landscape_,
               const DiseaseDispersal& disp_dis,
               const InsectPops& insects,
               const arma::mat& N0,
               const arma::mat& W0,
               const arma::mat& M0,
               const double& Y0,
               const arma::mat& wasp_attract_,
               const std::vector<uint64>& seeds);




    /*
     Iterate for one time point, and return bool for whether all plants
     are infected.
     */
    bool iterate(arma::umat& dispersals);



};





#endif
