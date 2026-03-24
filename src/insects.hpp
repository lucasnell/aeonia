# ifndef __AEONIA_INSECT_POPULATIONS_H
# define __AEONIA_INSECT_POPULATIONS_H

#include "aeonia_types.hpp"  // integer types


#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include <cmath>                // log, exp
#include <random>               // normal distribution
#include <cstdint>              // integer types
#include <algorithm>            // find
#include <pcg/pcg_random.hpp>   // pcg prng

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

    // Iterate aphid and wasp populations. For use in `test_insect_pops`
    void iterate(pcg32& eng) {

        arma::vec A_surv;

        // Set attack survival vector:
        double x = aphids.total_unpar();
        wasps.A_mats(A_surv, aphids.attack_surv, wasps.Y, x);

        // Now iterate aphids, then mummies
        double new_M = aphids.iterate(A_surv, eng);
        double new_Y = mummies.iterate(new_M);

        wasps.iterate(new_Y, eng);

        return;
    }

};








#endif
