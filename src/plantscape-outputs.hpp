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



/*
 This creates a character vector for output column names.
 This does the heavy lifting for the next two functions.
 */
CharacterVector col_namer__(const std::string& summ,
                            const bool& out_pseudo,
                            const bool& out_attack_surv,
                            const bool& out_stages);

/*
 Extract summary type based on a `ScapeSimmer` object, create column names, then
 verify that the output is the same length as `n_cols`
 */
CharacterVector col_namer_cpp(const ScapeSimmer& simmer,
                            const bool& out_pseudo,
                            const bool& out_attack_surv,
                            const bool& out_stages,
                            const uint32& n_cols);





// Convert from list of columns to a dataframe:
void list_to_data_frame(DataFrame& out_df,
                        const std::vector<std::vector<double>>& tmp_list,
                        const CharacterVector& col_names);




void ps_out_none(DataFrame& out_df,
                 const std::vector<ScapeSimmer>& simmers,
                 const arma::ucube& landscapes,
                 const uint32& max_t,
                 const bool& out_pseudo,
                 const bool& out_stages);





/*
 Note that `n_actual_rows` here should be the actual number of items in the
 output dataframe, not `n_rows` from the `ps_out_time` and `ps_out_all`
 functions since those are the max possible rows.
 There can be a mismatch between these values when simulations
 stop once all plants are infected.
 */
List make_disp_col(const std::vector<ScapeSimmer>& simmers,
                   const uint32& n_actual_rows,
                   const uint32& n_plants);




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
