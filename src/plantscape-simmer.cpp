#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng
#include <RcppThread.h>         // multithreading


#include "aeonia_types.hpp"         // integer types
#include "plantscape.hpp"           // PlantScape class
#include "plantscape-simmer.hpp"    // Functions to simulate plantscapes
#include "plantscape-outputs.hpp"   // Functions to create outputs
#include "pcg.hpp"                  // mt_seeds fxn
#include "util.hpp"                 // thread_check fxn



using namespace Rcpp;






void ScapeSimmer::run(const bool& infect_stop,
                      RcppThread::ProgressBar& prog_bar,
                      const bool& show_progress) {
    bool all_infected = false;
    for (uint32 t = 0; t < max_t; t++) {
        all_infected = scape.iterate(*disp_iter);
        // Move to next dispersal matrix if necessary:
        if (out_dispersals && summ == "time") disp_iter++;
        // Fill `output` with current conditions:
        fill_output();
        if (infect_stop && all_infected) {
            // Shorten dispersals if stopping early bc of full infection:
            if (out_dispersals && summ == "time" && disp_iter != dispersals.end()) {
                // Note: no need to add one to `curr_size` bc line after
                // `scape.iterate()` already iterates `disp_iter`
                size_t curr_size = disp_iter - dispersals.begin();
                dispersals.resize(curr_size);
            }
            break;
        }
        if (show_progress) prog_bar++;
        if (t % 10 == 0) RcppThread::checkUserInterrupt();
    }
    return;
}





/*
 ==========================================================================*
 Write the current state of the PlantScape to the `output` field.
 `summ` is for how to summarize output (if at all).
 `summ = "none"` isn't recommended for long time series of many plants!
 Note also that `summ = "time"` and `summ = "all"` result in the same
 output here, but they get summarized differently in the `sim_plantscape`
 function.
 Note lastly that for each of the column lists below,
 if `out_stages = TRUE`, there are two columns each for aphids and alates,
 one for juveniles and one for adults.

 If `summ == "none"`, `n_x * n_y` rows are output for each time point.
 The columns are...
 1) plant x
 2) plant y
 3) plant infectious with virus (0 or 1)
 4) aphid (non-winged) density
 5) alate density
 6) parasitized aphid density
 7) mummy density
 8) parasitoid density

 If `summ %in% c("time", "all")`, only one row is output per time point.
 The columns are...
 1) number of plants infectious with virus
 2) total aphid (non-winged) density summed across all plants
 3) total alate density summed across all plants
 4) total parasitized aphid density summed across all plants
 5) total mummy density summed across all plants
 6) total parasitoid density summed across all plants
 ==========================================================================*
 */
void ScapeSimmer::fill_output() {

    if (summ == "none") {
        // Densities (all vectors except wasps start with length 0 but
        // are reserved for length n_x * n_y,
        // wasps are a double, not a vector):
        output_dens.push_back(OutDensities());
        OutDensities& output_dens_t(output_dens.back());
        output_dens_t.reserve(n_x * n_y);
        // plant x, plant y:
        output_ids.push_back(arma::umat(n_x * n_y, 2, arma::fill::none));
        arma::umat& output_ids_t(output_ids.back());

        uint32 k = 0;
        output_dens_t.tot_wasps = scape.wasps.Y;
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {

                output_ids_t(k,0) = x+1U;
                output_ids_t(k,1) = y+1U;

                output_dens_t.push_back(static_cast<double>(scape.infectious[x][y]),
                                        scape.aphids[x][y].unparas_X(),
                                        scape.aphids[x][y].paras.total(),
                                        scape.mummies[x][y].total(),
                                        scape.Yi_mat(x,y));
                k++;
            }
        }

    } else if (summ == "time" || summ == "all") {

        // Densities (all vectors start with length 1, values = 0):
        output_dens.push_back(OutDensities(1, scape.aphids[0][0].n_stages()));
        OutDensities& output_dens_t(output_dens.back());
        // (no ids here)

        output_dens_t.tot_wasps = scape.wasps.Y;

        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                output_dens_t.add_to(0,
                                     static_cast<double>(scape.infectious[x][y]),
                                     scape.aphids[x][y].unparas_X(),
                                     scape.aphids[x][y].paras.total(),
                                     scape.mummies[x][y].total(),
                                     0.0);
            }
        }

    } else {

        throw std::runtime_error("INTERNAL ERROR: `! summ %in% c('none', 'time', 'all')`");

    }

    return;
}






//' Create a pointer object in which to store disease (and dispersal) info.
//'
//' This pointer is used as an argument to [sim_plantscape()].
//'
//' @details # Radius
//' From "The Role of Aphid Behaviour in the Epidemiology of Potato Virus Y:
//' a Simulation Study" by Thomas Nemecek (1993; p. 72), dispersal distances
//' follow a Weibull distribution with shape = 0.6569 and scale = 9.613.
//'
//' For use in larger landscapes, the object `pop_info$radius` is the median
//' of this distribution.
//' See `raw-data/pop_info.R` for the code used to generate it.
//'
//' For smaller landscapes (e.g., 3x3), I use `radius = 1`.
//'
//' @param n_x Single integer indicating x dimension of the landscape.
//'     Must be at least `2` and less than `1e6`.
//' @param n_y Single integer indicating y dimension of the landscape.
//'     Must be at least `2` and less than `1e6`.
//' @param radius Max distance that alates will travel between plants.
//'     Must be >= 1.
//'     See "Radius" section below for details.
//' @param virus_attract Effect of virus infection on alate alighting.
//'     Sampling weights for virus-infectious plants is `virus_attract`,
//'     compared to empty (no virus or *Pseudomonas*) plants whose weight is 1.
//'     Thus, when `virus_attract > 1`, alates are attracted to virus-infectious
//'     plants, while `virus_attract < 1` causes them to be repelled by
//'     virus-infectious plants.
//'     Values must be `> 0`.
//' @param pseudo_repel Effect of *Pseudomonas* infection on alate alighting.
//'     Sampling weights for *Pseudomonas*-containing plants is
//'     `1 / pseudo_repel` (note difference from `virus_attract`, hence the
//'     different names!), compared to empty (no virus or *Pseudomonas*)
//'     plants whose weight is 1.
//'     Thus, when `pseudo_repel > 1`, alates are repelled by
//'     *Pseudomonas*-containing plants,
//'     while `pseudo_repel < 1` causes them to be attracted to
//'     *Pseudomonas*-containing plants.
//'     Values must be `> 0`.
//' @param p_load_alate Single numeric indicating the probability that an
//'     uninoculated alate is loaded with a virus if it interacts with an
//'     inoculated plant.
//' @param p_load_plant Single numeric indicating the probability that an
//'     uninoculated plant is loaded with a virus if it interacts with an
//'     inoculated alate.
//' @param epsilon Effect of virus infection on alate acceptance.
//'     Values `> 1` cause alates to be more likely to stay and feed
//'     (indefinitely) on virus-infected plants,
//'     while values `< 1` cause them to be less likely to stay and feed on
//'     virus-infected plants.
//'     Values must be `> 0`, and `epsilon * w` must be `< 1`.
//'     Defaults to `1`.
//' @param w Probability that an alate accepts a plant, meaning that
//'     it stays to feed on it indefinitely.
//'     Must be `> 1e-4` and `< 1`.
//'     The lower limit is because a very small value of `w` causes alates to
//'     fly to so many plants that it becomes computationally problematic.
//'     This is the same reason that I check for (and return and error) if
//'     `w*epsilon < 1e-4`.
//'     Defaults to `0.2`.
//' @param total_exp_days Single integer indicating the number of days required
//'     for a plant to transition from exposed (inoculated with virus but
//'     not able to pass it on) to infectious (able to infect other plants).
//'     A value of `0` means that an alate inoculating a plant causes that
//'     plant to be infectious the same day.
//'     Because this is not realistic and computationally troublesome,
//'     only values `>1` are allowed.
//'     Defaults to `7`, which is based on the paper "Cucumber mosaic virus
//'     isolates seedborne in *Phaseolus vulgaris*: serology, host-pathogen
//'     relationships, and seed transmission" (Davis & Hampton, 1986).
//'
//'
//' @return An `externalptr` object that points to a C++ object that can
//' be passed to [sim_plantscape()].
//'
//' @export
//'
//[[Rcpp::export]]
SEXP make_disease_ptr(const double& radius,
                       const double& virus_attract,
                       const double& pseudo_repel,
                       const double& p_load_alate,
                       const double& p_load_plant,
                       const double& epsilon = 1,
                       const double& w = 0.2,
                       const uint32& total_exp_days = 7) {


    if (radius < 1) stop("radius < 1");
    if (virus_attract < 0) stop("virus_attract < 0");
    if (pseudo_repel < 0) stop("pseudo_repel < 0");

    if (p_load_alate < 0 || p_load_alate > 1) stop("p_load_alate < 0 || p_load_alate > 1");
    if (p_load_plant < 0 || p_load_plant > 1) stop("p_load_plant < 0 || p_load_plant > 1");

    if (epsilon < 0) stop("epsilon < 0");
    if (w < 0.0001 || w > 1) stop("w < 0.0001 || w > 1");
    if ((w*epsilon) > 1) stop("w*epsilon > 1");
    if ((w*epsilon) < 0.0001) stop("w*epsilon < 0.0001");

    if (total_exp_days < 1) stop("total_exp_days < 1");
    if (total_exp_days > 1e6) stop("total_exp_days > 1e6");

    // Set a maximum on # plants an alate can fly to such that the probability
    // that it reaches this threshold < 1e-9:
    uint32 max_fly_t = 100;
    double max_leave_p = 1 - std::min(w, epsilon * w); // max Pr(leave plant)
    while (std::pow(max_leave_p, max_fly_t) > (double)1e-9) {
        max_fly_t *= 10;
        if (max_fly_t > (uint32)1000000)
            throw std::runtime_error("INTERNAL ERROR: max_fly_t too high");
    }

    XPtr<DiseaseDispersal> disease_xptr(new DiseaseDispersal(
            radius, max_fly_t, virus_attract, pseudo_repel, epsilon, w,
            p_load_alate, p_load_plant, total_exp_days), true);

    return disease_xptr;

}












/*
 ===============================================================================
 ===============================================================================
 Checks arguments, adjusts starting densities if only 1 provided,
 adjusts # threads, fills `wasp_attract` from `wasp_plant_attract` if provided and 1 / (n_x*n_y) if not.
 ===============================================================================
 ===============================================================================
 */
void check_plantscape_args(const arma::ucube& landscapes,
                           const uint32& max_t,
                           arma::cube& N0,
                           arma::cube& W0,
                           arma::cube& M0,
                           arma::vec& Y0,
                           arma::mat& wasp_attract,
                           Nullable<NumericMatrix> wasp_plant_attract,
                           const std::string& summ,
                           uint32& infect_time_n,
                           const double& aphid_gone_thresh,
                           const double& wasp_gone_thresh,
                           uint32& n_threads) {

    uint32 n_x = landscapes.n_rows;
    uint32 n_y = landscapes.n_cols;
    uint32 n_reps = landscapes.n_slices;

    if (n_x < (uint32)2) stop("nrow(landscapes) < 2");
    if (n_y < (uint32)2) stop("ncol(landscapes) < 2");
    if (n_x > (uint32)1e6) stop("nrow(landscapes) > 1e6");
    if (n_y > (uint32)1e6) stop("ncol(landscapes) > 1e6");
    if (n_reps < (uint32)1) stop("dim(landscapes)[3] < 1");
    if (n_reps > (uint32)1e6) stop("dim(landscapes)[3] > 1e6");
    if (arma::any(arma::vectorise(landscapes) > 3))
        stop("landscapes cannot contain values > 3");
    if (arma::any(arma::vectorise(N0) < 0)) stop("N0 cannot contain values < 0");
    if (arma::any(arma::vectorise(W0) < 0)) stop("W0 cannot contain values < 0");
    if (arma::any(arma::vectorise(M0) < 0)) stop("M0 cannot contain values < 0");
    if (arma::any(Y0 < 0)) stop("Y0 cannot contain values < 0");

    if (max_t == 0 || max_t > 1e6) stop("max_t == 0 || max_t > 1e6");

    if (N0.n_rows == 1) {
        if (N0.n_cols != 1 || N0.n_slices != 1)
            stop("dim(N0) must be c(1,1,1) or dim(landscapes)");
        double val = N0(0,0,0);
        N0.set_size(arma::size(landscapes));
        N0.fill(val);
    } else if (arma::size(N0) != arma::size(landscapes)) {
        stop("dim(N0) must be c(1,1,1) or dim(landscapes)");
    }

    if (W0.n_rows == 1) {
        if (W0.n_cols != 1 || W0.n_slices != 1)
            stop("dim(W0) must be c(1,1,1) or dim(landscapes)");
        double val = W0(0,0,0);
        W0.set_size(arma::size(landscapes));
        W0.fill(val);
    } else if (arma::size(W0) != arma::size(landscapes)) {
        stop("dim(W0) must be c(1,1,1) or dim(landscapes)");
    }

    if (M0.n_rows == 1) {
        if (M0.n_cols != 1 || M0.n_slices != 1)
            stop("dim(M0) must be c(1,1,1) or dim(landscapes)");
        double val = M0(0,0,0);
        M0.set_size(arma::size(landscapes));
        M0.fill(val);
    } else if (arma::size(M0) != arma::size(landscapes)) {
        stop("dim(M0) must be c(1,1,1) or dim(landscapes)");
    }

    if (Y0.n_elem == 1) {
        double val = Y0(0);
        Y0.set_size(n_reps);
        Y0.fill(val);
    } else if (Y0.n_elem != n_reps) stop("length(Y0) must be 1 or dim(landscapes)[3]");

    if (summ != "none" && summ != "time" && summ != "all") {
        stop("`summ` should be 'none', 'time', or 'all'");
    }
    // Next line is equivalent to ceiling(n_x*n_y / 2):
    if (infect_time_n == 0) infect_time_n = 1 + ((n_x * n_y - 1) / 2U);
    if (infect_time_n > n_x * n_y) stop("infect_time_n > n_x * n_y");
    if (aphid_gone_thresh <= 0) stop("aphid_gone_thresh <= 0");
    if (wasp_gone_thresh <= 0) stop("wasp_gone_thresh <= 0");

    wasp_attract.set_size(n_x, n_y);
    if (wasp_plant_attract.isNotNull()){
        NumericMatrix wpa(wasp_plant_attract);
        if (wpa.nrow() != n_x) stop("nrow(wasp_plant_attract) != nrow(landscapes)");
        if (wpa.ncol() != n_y) stop("ncol(wasp_plant_attract) != ncol(landscapes)");
        double wpa_sum = 0;
        for (uint32 x = 0; x < n_x; x++) {
            for (uint32 y = 0; y < n_y; y++) {
                if (wpa(x,y) < 0) stop("any(wasp_plant_attract < 0)");
                wasp_attract(x,y) = wpa(x,y);
                wpa_sum += wpa(x,y);
            }
        }
        if (wpa_sum != 1) wasp_attract /= wpa_sum;
    } else wasp_attract.fill(1.0 / static_cast<double>(wasp_attract.n_elem));

    thread_check(n_threads); // Check that # threads isn't too high

    // Make sure output object won't be too big for R:
    // if  (R_output) {
    if  (summ != "all") {
        uint32 n_rows = n_reps * (max_t + (uint32)1U);
        if (summ == "none") n_rows *= (n_x * n_y);
        if (n_rows > (uint32)2147483647)
            stop("This combo of parameters will produce too large of an output for R");
    }


    return;

}






//' Simulation plant landscapes with virus spread.
//'
//' @details # Summarizing
//' If `summ == "none"`, `n_x * n_y` rows are output for each rep and time point.
//' The columns are...
//' 1.  `rep`: repetition number
//' 2.  `time`: time point
//' 3.  `x`: plant x coordinate
//' 4.  `y`: plant y coordinate
//' 5.  `pseudo`: plant contains *Pseudomonas* (0 or 1); note: this column isn't
//'     included if `out_pseudo = FALSE`
//' 6.  `virus`: plant infectious with virus (0 or 1)
//' 7.  `aphids`: aphid (non-winged) density
//' 8.  `alates`: alate density
//' 9.  `parasitized`: parasitized aphid density
//' 10. `mummies`: mummy density
//' 11. `wasps`: total parasitoid density (`x` and `y` for this row will be 0)
//'
//'
//' If `summ == "time"`, one row is output per rep and time point.
//' The columns are...
//' 1.  `rep`: repetition number
//' 2.  `time`: time point
//' 3. `virus`: number of plants infectious with virus
//' 4. `aphids`: total aphid (non-winged) density summed across all plants
//' 5. `alates`: total alate density summed across all plants
//' 6. `parasitized`: total parasitized aphid density summed across all plants
//' 7. `mummies`: total mummy density summed across all plants
//' 8. `wasps` total parasitoid density
//'
//' If `summ == "all"`, one row is output per rep.
//' The columns are...
//' 1.  `rep`: repetition number
//' 2.  `p_alates`: mean total alates / total aphids
//' 3.  `log_aphids`: mean log(aphids+1)
//' 4.  `aphids`: mean aphids
//' 5.  `log_alates`: mean log(alates+1)
//' 6.  `alates`: mean alates
//' 7.  `log_alates`: mean log(parasitized aphids+1)
//' 8.  `alates`: mean parasitized aphids
//' 9.  `log_mummies`: mean log(mummies+1)
//' 10. `mummies`: mean mummies
//' 11. `log_wasps`: mean log(wasps+1)
//' 12. `wasps`: mean wasps
//' 13. `infect_time`: time it took to have `infect_time_n` plants infected
//' 14. `outbreak_size`: maximum number of plants infected with virus
//'
//'
//'
//'
//' @param landscapes Integer cube with the types of each plant in each
//'     landscape. It's assumed that rows are x the dimension,
//'     columns are the y dimension, and "slices" (i.e., `dim(landscapes)[3]`)
//'     indicate separate landscapes.
//'     Values in each cell give the state of the plant:
//'     `0` indicates nothing on plant,
//'     `1` indicates just virus on plant (infectious),
//'     `2` indicates just *Pseudomonas* on plant,
//'     `3` indicates both virus and *Pseudomonas* on plant.
//'     Values < 0 or > 3 are not allowed.
//'     Note that this array is coerced to an array of unsigned integers, so
//'     negative values will become very large integers.
//'     Hence, do not expect an error for negative numbers if you pass them here.
//' @param max_t Single integer giving the maximum time the simulations run.
//' @param N0 Numeric 3D array indicating the starting aphid (non-winged) population
//'     density for each plant and rep.
//'     To indicate separate densities for each plant and/or rep,
//'     the array should have the same dimensions as `landscapes`.
//'     The array can also be 1x1, in which case it's assumed that all plants
//'     and reps start with the same density of aphids.
//' @param W0 Numeric 3D array indicating the starting winged aphid population
//'     density for each plant and rep.
//'     To indicate separate densities for each plant and/or rep,
//'     the array should have the same dimensions as `landscapes`.
//'     The array can also be 1x1, in which case it's assumed that all plants
//'     and reps start with the same density of winged aphids.
//' @param M0 Numeric 3D array indicating the starting mummy population
//'     density for each plant and rep.
//'     To indicate separate densities for each plant and/or rep,
//'     the array should have the same dimensions as `landscapes`.
//'     The array can also be 1x1, in which case it's assumed that all plants
//'     and reps start with the same density of mummies.
//' @param Y0 Numeric vector indicating the starting parasitoid population
//'     density for each rep.
//'     The vector can also be of length 1, in which case it's assumed that
//'     all reps start with the same density of parasitoids.
//' @param insect_ptr External pointer to a C++ object with insect population
//'     information, output from function [make_insects_ptr()].
//' @param disease_ptr External pointer to a C++ object with disease (and
//'     dispersal) information, output from function [make_disease_ptr()].
//' @param wasp_plant_attract Relative attractiveness of plants to wasps.
//'     This affects the proportion of wasps that immigrate from the dispersal
//'     pool to each plant.
//'     It doesn't change the number of wasps that leave plants.
//'     This must be `NULL` or a matrix the same dimensions as each slice
//'     of `landscapes`.
//'     If `NULL`, all plants are equally attractive to wasps.
//'     If a matrix is provided, then the values are divided by their sum
//'     (to make it sum to 1), then those values are used as the proportion of
//'     wasps from the dispersal pool that immigrate to each plant.
//'     Defaults to `NULL`.
//' @param summ Single string to indicate how to summarize output.
//'     If `summ == "none"`, then no summarizing is done, so output is separate
//'     by plant.
//'     If `summ == "time"`, then output is summarized across all plants, so
//'     will only be separated by time and rep.
//'     If `summ == "all"`, then output is summarized across all plants and
//'     time points, so will only be separated by rep.
//'     See "Summarizing" section below for details on output columns.
//'     Defaults to `"none"`.
//' @param infect_time_n Single integer specifying the number of plants
//'     to use when calculating the `infect_time` column when `summ == "all"`.
//'     Ignored when `summ != "all"` except for error checking.
//'     Note that since this is coerced to an unsigned integer, using
//'     a negative number here could cause an error to occur.
//'     Defaults to `0`, which results in at least half the total number of
//'     plants (i.e., `ceiling(prod(dim(landscapes)[1:2]) / 2)`) being used.
//' @param aphid_gone_thresh Single numeric specifying the threshold for aphid
//'     abundance (all stages + parasitized, summed across all plants)
//'     below which it's considered gone when calculating the
//'     `aphid_gone_n` column when `summ == "all"`.
//'     Ignored when `summ != "all"` except for error checking.
//'     Defaults to `1`.
//' @param wasp_gone_thresh Single numeric specifying the threshold for wasp
//'     abundance below which it's considered gone when calculating the
//'     `wasp_gone_n` column when `summ == "all"`.
//'     Ignored when `summ != "all"` except for error checking.
//'     Defaults to `1`.
//' @param infect_stop Single logical for whether to stop simulations
//'     when all plants are infected with virus.
//'     Defaults to `TRUE`.
//' @param out_pseudo Single logical for whether to include *Pseudomonas*
//'     presence in output. Ignored if `summ != "none"`.
//'     Defaults to `FALSE`.
//' @param out_attack_surv Single logical for whether to include aphid survival
//'     from parasitoid attack (weighted mean by abundance for each aphid stage)
//'     in output. Ignored if `summ != "none"`.
//'     Defaults to `FALSE`.
//' @param out_stages Single logical for whether to separate output for aphids
//'     by juvenile vs adults.
//'     Defaults to `FALSE`.
//' @param out_dispersals Single logical for whether to output a column list
//'     containing matrices with the number of alate dispersals connecting
//'     plants. The column indicates the plant the alate came from,
//'     and the row indicates the plant the alate dispersed to.
//'     This argument only does something when `summ %in% c("time", "all")`.
//'     Defaults to `FALSE`.
//' @param show_progress Single logical for whether to show progress bar.
//'     Defaults to `FALSE`.
//' @param n_threads Single integer for the number of threads to use.
//'     Ignored if `dim(landscapes)[3] == 1`.
//'     Defaults to `0L`, in which case it takes the value from
//'     `options("mc.cores")` if assigned, and `1L` if not assigned.
//'
//'
//' @export
//'
//' @return A tibble with columns following the description in the
//' "Summarizing" section.
//'
//'
//' @importFrom tibble as_tibble
//'
//'
//[[Rcpp::export]]
DataFrame sim_plantscape(const arma::ucube& landscapes,
                         arma::cube N0,
                         arma::cube W0,
                         arma::cube M0,
                         arma::vec Y0,
                         SEXP insect_ptr,
                         SEXP disease_ptr,
                         Nullable<NumericMatrix> wasp_plant_attract = R_NilValue,
                         const uint32& max_t = 100,
                         const std::string& summ = "none",
                         uint32 infect_time_n = 0,
                         const double& aphid_gone_thresh = 1,
                         const double& wasp_gone_thresh = 1,
                         const bool& infect_stop = false,
                         const bool& out_pseudo = false,
                         const bool& out_attack_surv = false,
                         const bool& out_stages = false,
                         const bool& out_dispersals = false,
                         const bool& show_progress = false,
                         uint32 n_threads = 0) {

    arma::mat wasp_attract;
    check_plantscape_args(landscapes, max_t, N0, W0, M0, Y0, wasp_attract,
                          wasp_plant_attract, summ, infect_time_n,
                          aphid_gone_thresh, wasp_gone_thresh, n_threads);

    uint32 n_reps = landscapes.n_slices;

    // Base for all insect populations to start with:
    XPtr<InsectPops> insect_xptr(insect_ptr);
    const InsectPops& insects(*insect_xptr);
    // Disease (and dispersal) info:
    XPtr<DiseaseDispersal> disease_xptr(disease_ptr);
    const DiseaseDispersal& disease(*disease_xptr);

    std::vector<std::vector<uint64>> seeds = mt_seeds(n_reps);
    std::vector<ScapeSimmer> simmers;
    simmers.reserve(n_reps);

    for (uint32 i = 0; i < n_reps; i++) {
        simmers.emplace_back(landscapes.slice(i), disease, insects,
                             N0.slice(i), W0.slice(i), M0.slice(i),
                             Y0(i), wasp_attract, seeds[i],
                             summ, max_t, out_dispersals);
    }

    RcppThread::ProgressBar prog_bar(n_reps * max_t, 1);

    if (n_threads > 1U && n_reps > 1U) {
        auto job = [&] (ScapeSimmer& simmer) {
            simmer.run(infect_stop, prog_bar, show_progress);
        };
        RcppThread::parallelForEach(simmers, job, n_threads);
    } else {
        for (uint32 i = 0; i < n_reps; i++) {
            simmers[i].run(infect_stop, prog_bar, show_progress);
        }
    }

    // Produce output dataframe:
    DataFrame out_df;

    if (summ == "none") {
        ps_out_none(out_df, simmers, landscapes, max_t,
                    out_pseudo, out_attack_surv, out_stages);
    } else if (summ == "time") {
        ps_out_time(out_df, simmers, landscapes, max_t, out_stages);
    } else if (summ == "all") {
        ps_out_all(out_df, simmers, landscapes, max_t, out_stages,
                   infect_time_n, aphid_gone_thresh, wasp_gone_thresh);
    } else {
        stop("INTERNAL ERROR: `! summ %in% c('none', 'time', 'all')`");
    }

    return out_df;

}








