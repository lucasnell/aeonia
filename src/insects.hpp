# ifndef __AEONIA_INSECT_POPULATIONS_H
# define __AEONIA_INSECT_POPULATIONS_H


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <cmath>                // log, exp
#include <random>               // normal distribution
#include <cstdint>              // integer types
#include <algorithm>            // find
#include <pcg/pcg_random.hpp>   // pcg prng
#include "aeonia_types.hpp"  // integer types
#include "aphids.hpp"  // aphid classes
#include "wasps.hpp"  // wasp classes



using namespace Rcpp;




struct InsectPops {

    AphidPop aphids;
    WaspPop wasps;
    MummyPop mummies;

    InsectPops(const AphidPop& aphids_,
               const WaspPop& wasps_,
               const MummyPop& mummies_)
        : aphids(aphids_), wasps(wasps_), mummies(mummies_) {};

};








#endif
