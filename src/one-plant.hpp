#ifndef __AEONIA_ONE_PLANT_H
#define __AEONIA_ONE_PLANT_H

/*
 This contains code for one plant.
 */

#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng

#include "aeonia_types.hpp"     // integer types
#include "insect-pops.hpp"      // AphidPops class
#include "pcg.hpp"              // runif_01 fxn



using namespace Rcpp;





struct OnePlant {

    // virus: exposed but not yet infectious?
    bool exposed;
    // virus: infectious?
    bool infectious;
    // contains Pseudomonas?
    bool pseudo;
    // days since exposure (ignored if not exposed):
    uint32 exp_days;

    // Aphid populations:
    AphidPops aphids;

    OnePlant(const bool& infectious_,
             const bool& pseudo_,
             const uint32& total_exp_days_,
             const AphidPops& aphids_)
        : exposed(false),
          infectious(infectious_),
          pseudo(pseudo_),
          exp_days(0),
          aphids(aphids_),
          total_exp_days(total_exp_days_) {
        if (!pseudo) aphids.set_pseudo_surv(1.0);
    }

    // iterate and (1) set number of alates moving from this plant and
    // (2) add to the number of new parasitoids (both male and female)
    void iterate(const double& Yi,
                 uint32& n_alates,
                 double& new_Ys,
                 pcg32& eng) {

        aphids.iterate(Yi, n_alates, new_Ys, eng);

        if (exposed) {
            exp_days++;
            if (exp_days >= total_exp_days) {
                infectious = true;
                exposed = false;
                exp_days = 0;
            }
        }
        return;
    }



private:

    uint32 total_exp_days; // days since exposure required to transition to infected

};



#endif
