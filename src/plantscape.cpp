#include <RcppArmadillo.h>
#include <vector>
#include <math.h>
#include <algorithm>
#include <random>
#include <pcg/pcg_random.hpp>   // pcg prng


#include "aeonia_types.hpp"     // integer types
#include "plantscape.hpp"       // PlantScape class
#include "plantscape-outputs.hpp" // Functions to create outputs
#include "pcg.hpp"              // mt_seeds fxn
#include "util.hpp"             // thread_check fxn



using namespace Rcpp;





/*
 Does most of the work of sim_plantscape functions, including checking
 validity of arguments.
 The only thing it doesn't do is produce output.
 */
std::vector<PlantScape> sim_plantscape_cpp(const bool& R_output,
                                           const arma::ucube& landscapes,
                                           const uint32& max_t,
                                           SEXP insect_ptr,
                                           arma::cube& N0,
                                           arma::cube& W0,
                                           arma::vec& Y0,
                                           const double& alpha,
                                           const double& beta,
                                           const double& epsilon,
                                           const double& delta_a,
                                           const double& delta_p,
                                           const uint32& total_exp_days,
                                           const double& w,
                                           const double& radius,
                                           Nullable<NumericMatrix> wasp_plant_attract,
                                           const std::string& summ,
                                           uint32& infect_time_n,
                                           const bool& infect_stop,
                                           const bool& out_pseudo,
                                           const bool& show_progress,
                                           uint32& n_threads) {

    uint32 n_x = landscapes.n_rows;
    uint32 n_y = landscapes.n_cols;
    uint32 n_reps = landscapes.n_slices;

    if (n_x < (uint32)2) stop("nrow(landscapes) < 2");
    if (n_y < (uint32)2) stop("ncol(landscapes) < 2");
    if (n_reps < (uint32)1) stop("dim(landscapes)[3] < 1");
    if (arma::any(arma::vectorise(landscapes) > 3))
        stop("landscapes cannot contain values > 3");

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

    if (Y0.n_elem == 1) {
        double val = Y0(0);
        Y0.set_size(n_reps);
        Y0.fill(val);
    } else if (Y0.n_elem != n_reps) stop("length(Y0) must be 1 or dim(landscapes)[3]");

    if (epsilon < 0) stop("epsilon < 0");
    if (delta_a < 0 || delta_a > 1) stop("delta_a < 0 || delta_a > 1");
    if (delta_p < 0 || delta_p > 1) stop("delta_p < 0 || delta_p > 1");
    if (total_exp_days < 1) stop("total_exp_days < 1");
    if (total_exp_days > 1e6) stop("total_exp_days > 1e6");
    if (w < 0.0001 || w > 1) stop("w < 0.0001 || w > 1");
    if ((w*epsilon) > 1) stop("w*epsilon > 1");
    if ((w*epsilon) < 0.0001) stop("w*epsilon < 0.0001");
    if (radius < 1) stop("radius < 1");
    if (summ != "none" && summ != "pseudo" && summ != "time" && summ != "all") {
        stop("`summ` should be 'none', 'pseudo', 'time', or 'all'");
    }
    if (infect_time_n == 0) infect_time_n = n_x * n_y;
    if (infect_time_n > n_x * n_y) stop("infect_time_n > n_x * n_y");

    arma::mat wasp_attract(n_x, n_y, arma::fill::none);
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
    } else wasp_attract.fill(1.0 / static_cast<double>(n_x * n_y));

    thread_check(n_threads); // Check that # threads isn't too high

    // Set a maximum on # plants an alate can fly to such that the probability
    // that it reaches this threshold < 1e-9:
    uint32 max_fly_t = 100;
    double max_leave_p = 1 - std::min(w, epsilon * w); // max Pr(leave plant)
    while (std::pow(max_leave_p, max_fly_t) > (double)1e-9) {
        max_fly_t *= 10;
        if (max_fly_t > (uint32)1000000) stop("INTERNAL ERROR: max_fly_t too high");
    }

    // Make sure output object won't be too big for R:
    if  (R_output) {
        uint32 n_rows = n_reps * (max_t + (uint32)1U);
        if (summ == "none") n_rows *= (n_x * n_y);
        if (summ == "pseudo") n_rows *= (uint32)2U;
        if (n_rows > (uint32)2147483647)
            stop("This combo of parameters will produce too large of an output for R");
    }

    // Base for all insect populations to start with:
    XPtr<InsectPops> insect_xptr(insect_ptr);
    const InsectPops& insects0(*insect_xptr);

    std::vector<std::vector<uint64>> seeds = mt_seeds(n_reps);
    std::vector<PlantScape> plantscapes;
    plantscapes.reserve(n_reps);

    for (uint32 i = 0; i < n_reps; i++) {
        plantscapes.push_back(PlantScape(max_t, summ, max_fly_t,
                                         landscapes.slice(i), radius, alpha,
                                         beta, epsilon, w, wasp_attract,
                                         delta_a, delta_p,
                                         total_exp_days, insects0,
                                         N0.slice(i), W0.slice(i), Y0(i),
                                         seeds[i]));
    }

    RcppThread::ProgressBar prog_bar(n_reps * max_t, 1);

    if (n_threads > 1U && n_reps > 1U) {
        auto job = [&] (PlantScape& plantscape) {
            plantscape.run(infect_stop, prog_bar, show_progress);
        };
        RcppThread::parallelForEach(plantscapes, job, n_threads);
    } else {
        for (uint32 i = 0; i < n_reps; i++) {
            plantscapes[i].run(infect_stop, prog_bar, show_progress);
        }
    }

    return plantscapes;

}






//' Simulation plant landscapes with virus spread.
//'
//' @details # Radius
//' From "The Role of Aphid Behaviour in the Epidemiology of Potato Virus Y:
//' a Simulation Study" by Thomas Nemecek (1993; p. 72), dispersal distances
//' follow a Weibull distribution with shape = 0.6569 and scale = 9.613.
//'
//' The default for the `radius` argument uses the median of this
//' distribution.
//' I'm dividing by 0.75 to convert from meters to plant locations that are
//' 0.75 meters apart (typical spacing for pea):
//' `radius = qweibull(0.5, 0.6569, 9.613) / 0.75`.
//'
//' # Summarizing
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
//' If `summ == "pseudo"`, two rows are output for each rep and time point.
//' The columns are...
//' 1.  `rep`: repetition number
//' 2.  `time`: time point
//' 3. `pseudo`: plant contains *Pseudomonas* (0 or 1)
//' 4. `n`: total plants of this type
//' 5. `virus`: number of plants infectious with virus
//' 6. `aphids`: total aphid (non-winged) density summed across all plants of
//'    this type
//' 7. `alates`: total alate density summed across all plants of this type
//' 8. `parasitized`: total parasitized aphid density summed across all plants
//'    of this type
//' 9. `mummies`: total mummy density summed across all plants of this type
//' 10.`wasps`:  total parasitoid density (`type` for this row will be `2`)
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
//' @param max_t Single integer giving the maximum time the simulations run.
//' @param insect_ptr External pointer to a C++ object with insect population
//'     information, output from function [make_insect_ptr()].
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
//' @param Y0 Numeric vector indicating the starting parasitoid population
//'     density for each rep.
//'     The vector can also be of length 1, in which case it's assumed that
//'     all reps start with the same density of parasitoids.
//' @param alpha Effect of virus infection on alate alighting.
//'     Values `> 0` cause alates to be attracted to virus-infected plants,
//'     while values `< 0` cause them to be repelled by virus-infected plants.
//' @param beta Effect of *Pseudomonas* infection on alate alighting.
//'     Values `> 0` cause alates to be attracted to *Pseudomonas*-infected plants,
//'     while values `< 0` cause them to be repelled by *Pseudomonas*-infected plants.
//' @param epsilon Effect of virus infection on alate acceptance.
//'     Values `> 1` cause alates to be more likely to stay and feed
//'     (indefinitely) on virus-infected plants,
//'     while values `< 1` cause them to be less likely to stay and feed on
//'     virus-infected plants.
//'     Values must be `> 0`, and `epsilon * w` must be `< 1`.
//' @param delta_a Single numeric indicating the probability that an
//'     uninoculated alate is loaded with a virus if it interacts with an
//'     inoculated plant.
//' @param delta_p Single numeric indicating the probability that an
//'     uninoculated plant is loaded with a virus if it interacts with an
//'     inoculated alate.
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
//' @param w Probability that an alate accepts a plant, meaning that
//'     it stays to feed on it indefinitely.
//'     Must be `> 1e-4` and `< 1`.
//'     The lower limit is because a very small value of `w` causes alates to
//'     fly to so many plants that it becomes computationally problematic.
//'     This is the same reason that I check for (and return and error) if
//'     `w*epsilon < 1e-4`.
//'     Defaults to `0.2`.
//' @param radius Max distance that alates will travel between plants.
//'     Defaults to `7.336451`, which is based on previous work.
//'     See "Radius" section below for details.
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
//'     If `summ == "pseudo"`, then output is summarized by whether plants
//'     contain *Pseudomonas*.
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
//'     Defaults to `0`, which results in the total number of plants
//'     (i.e., `prod(dim(landscapes)[1:2])`) being used.
//' @param infect_stop Single logical for whether to stop simulations
//'     when all plants are infected with virus.
//'     Defaults to `TRUE`.
//' @param out_pseudo Single logical for whether to include *Pseudomonas*
//'     presence in output. Ignored if `summ != "none"`.
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
//[[Rcpp::export]]
DataFrame sim_plantscape(const arma::ucube& landscapes,
                         const uint32& max_t,
                         SEXP insect_ptr,
                         arma::cube N0,
                         arma::cube W0,
                         arma::vec Y0,
                         const double& alpha,
                         const double& beta,
                         const double& epsilon,
                         const double& delta_a,
                         const double& delta_p,
                         const uint32& total_exp_days = 7,
                         const double& w = 0.2,
                         const double& radius = 7.336451,
                         Nullable<NumericMatrix> wasp_plant_attract = R_NilValue,
                         const std::string& summ = "none",
                         uint32 infect_time_n = 0,
                         const bool& infect_stop = true,
                         const bool& out_pseudo = false,
                         const bool& show_progress = false,
                         uint32 n_threads = 0) {

    std::vector<PlantScape> plantscapes =
        sim_plantscape_cpp(true, landscapes, max_t, insect_ptr, N0, W0, Y0,
                           alpha, beta, epsilon, delta_a, delta_p,
                           total_exp_days, w, radius, wasp_plant_attract, summ,
                           infect_time_n,
                           infect_stop, out_pseudo, show_progress, n_threads);

    // Produce output dataframe:
    DataFrame out_df;

    if (summ == "none") {
        ps_out_none(out_df, plantscapes, landscapes, max_t, out_pseudo);
    } else if (summ == "pseudo") {
        ps_out_pseudo(out_df, plantscapes, landscapes, max_t);
    } else if (summ == "time") {
        ps_out_time(out_df, plantscapes, landscapes, max_t);
    } else if (summ == "all") {
        ps_out_all(out_df, plantscapes, landscapes, max_t, infect_time_n);
    } else {
        stop("INTERNAL ERROR: `! summ %in% c('none', 'pseudo', 'time', 'all')`");
    }

    return out_df;

}




