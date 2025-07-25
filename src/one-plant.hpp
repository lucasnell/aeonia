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
#include "insect-pops.hpp"      // InsectPops class
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

    // Insect populations:
    InsectPops insects;

    OnePlant(const bool& infectious_,
             const bool& pseudo_,
             const uint32& total_exp_days_,
             const InsectPops& insects_)
        : exposed(false),
          infectious(infectious_),
          pseudo(pseudo_),
          exp_days(0),
          insects(insects_),
          total_exp_days(total_exp_days_) {
        if (!pseudo) insects.set_B(0.0);
    }

    // iterate and set number of alates moving from this plant:
    void iterate(uint32& n_alates, double& wasp_disp_pool, pcg32& eng) {
        insects.iterate(n_alates, wasp_disp_pool, eng);
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
