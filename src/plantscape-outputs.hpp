#ifndef __AEONIA_PLANTSCAPE_OUTPUTS_H
#define __AEONIA_PLANTSCAPE_OUTPUTS_H


#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng


#include "aeonia_types.hpp"         // integer types
#include "insects.hpp"              // InsectPops class
#include "plantscape.hpp"           // PlantScape class
#include "plantscape-simmer.hpp"    // ScapeSimmer class
#include "pcg.hpp"                  // mt_seeds fxn
#include "util.hpp"                 // thread_check fxn



using namespace Rcpp;





void ps_out_none(DataFrame& out_df,
                 const std::vector<ScapeSimmer>& simmers,
                 const arma::ucube& landscapes,
                 const uint32& max_t,
                 const bool& out_pseudo,
                 const bool& out_attack_surv,
                 const bool& out_stages);




void ps_out_time(DataFrame& out_df,
                 const std::vector<ScapeSimmer>& simmers,
                 const arma::ucube& landscapes,
                 const uint32& max_t,
                 const bool& out_stages);



void ps_out_all(DataFrame& out_df,
                const std::vector<ScapeSimmer>& simmers,
                const arma::ucube& landscapes,
                const uint32& max_t,
                const bool& out_stages,
                const uint32& infect_time_n,
                const double& aphid_gone_thresh,
                const double& wasp_gone_thresh);








#endif
