
#include "sim-plant-types.hpp"




using namespace Rcpp;



//' Simulate field(s) of plant types.
//'
//' Simulate locations of plants of different types
//' (virus-infected or -uninfected, *Pseudomonas*-inhabited or uninhabited)
//' along an evenly spaced grid of integers, where the placement of one
//' target can affect subsequent placement of other targets.
//' Target locations are drawn from all combinations of `1` to
//' `n_x` and `1` to `n_y`.
//'
//' @details # Weighting
//' All `wt_*` parameters refer to how placement of one type of plant
//' (virus-infected or -uninfected, *Pseudomonas*-inhabited or uninhabited)
//' affects the sampling weight of subsequent plants of different types.
//' Weight values above 1 cause neighboring locations to be more likely to be
//' sampled later, while values below 1 cause them to be less likely
//' sampled.
//' These weights apply to all plants directly next to the plant of a given
//' type, and they accrue multiplicatively.
//' For example, if `wt_vp = 2` and *Pseudomonas* was place on a plant at
//' location `1,1` (x,y), then plants at locations `1,1`; `1,2`; `2,2`; and `2,1`
//' would all be twice as likely to be sampled as a location for a virus to be
//' placed later on.
//' In the same situation, if *Pseudomonas* was next placed at location `2,1`,
//' then those same locations would now be four times as likely to be chosen
//' for a virus; locations `3,1` and `3,2` would now be twice as likely.
//'
//'
//' @param n_x Single integer indicating x dimension of search area.
//'     Locations will be drawn from `1` to `n_x`.
//'     Must be at least `2` and less than `1e6`.
//' @param n_y Single integer indicating y dimension of search area.
//'     Locations will be drawn from `1` to `n_y`.
//'     Must be at least `2` and less than `1e6`.
//' @param wt_vv A single number indicating how virus placement
//'     affects subsequent placement of virus.
//'     See section "Weighting" above for how these weights work.
//'     Must be >= 0. Defaults to `1`.
//' @param wt_pp A single number indicating how *Pseudomonas* placement
//'     affects subsequent placement of *Pseudomonas*.
//'     See section "Weighting" above for how these weights work.
//'     Must be >= 0. Defaults to `1`.
//' @param wt_vp A single number indicating how virus placement
//'     affects subsequent placement of *Pseudomonas*.
//'     See section "Weighting" above for how these weights work.
//'     Must be >= 0. Defaults to `1`.
//' @param wt_pv A single number indicating how *Pseudomonas* placement
//'     affects subsequent placement of virus.
//'     See section "Weighting" above for how these weights work.
//'     Must be >= 0. Defaults to `1`.
//' @param n_virus Single integer indicating the number of plants that should
//'     be virus infected. This value must be `>=1` and `<= n_x * n_y`.
//' @param n_pseudo Single integer indicating the number of plants that should
//'     be *Pseudomonas* infected.
//'     This value must be `>=1` and `<= n_x * n_y`.
//' @param n_lands Single integer indicating the number of independent
//'     landscapes to simulate. Each landscape will have a number of samples
//'     per type according to `n_virus` and `n_pseudo`.
//'     Must be `> 0` and `< 1e6`.
//'     Defaults to `1`.
//' @param show_progress Single logical for whether to show progress bar.
//'     Defaults to `FALSE`.
//' @param n_threads Single integer for the number of threads to use.
//'     Multithreading happens across different landscapes, so this is
//'     ignored if `n_lands == 1`.
//'     Defaults to `0L`, in which case it takes the value from
//'     `options("mc.cores")` if assigned, and `1L` if not assigned.
//'
//' @return A cube where rows indicate the x coordinate, columns indicate the
//'     y coordinate, and slices indicate the landscape number.
//'     For each item,
//'     `0` indicates nothing,
//'     `1` indicates a virus infected plant (no *Pseudomonas*),
//'     `2` indicates a *Pseudomonas*-infected plant (no virus), and
//'     `3` indicates a plant with both virus and *Pseudomonas*.
//'     See [land_cube2list()] for converting this to a list of length `n_lands`,
//'     where each item is a [`tibble`][tibble::tbl_df] with columns
//'     `x`, `y`, and `type`.
//'
//'
//' @export
//'
//[[Rcpp::export]]
arma::icube sim_plant_types(const uint32& n_x,
                            const uint32& n_y,
                            const uint32& n_virus,
                            const uint32& n_pseudo,
                            const double& wt_vp = 1,
                            const double& wt_pv = 1,
                            const double& wt_vv = 1,
                            const double& wt_pp = 1,
                            const uint32& n_lands = 1,
                            Nullable<IntegerMatrix> virus_starts = R_NilValue,
                            Nullable<IntegerMatrix> pseudo_starts = R_NilValue,
                            const bool& show_progress = false,
                            uint32 n_threads = 0) {

    if (n_x < 1 || n_x > 1e6) stop("n_x must be >= 1 and <= 1e6");
    if (n_y < 1 || n_y > 1e6) stop("n_y must be >= 1 and <= 1e6");
    if (n_x * n_y < 2) stop("n_x * n_y must be >= 2");
    if ((n_virus + n_pseudo) < 1) {
        std::string err_msg = "n_virus + n_pseudo cannot be < 1; ";
        err_msg += "just use an array of zeros instead of using this function.";
        stop(err_msg.c_str());
    }
    if (n_virus > n_x * n_y) stop("n_virus cannot be > n_x * n_y");
    if (n_pseudo > n_x * n_y) stop("n_pseudo cannot be > n_x * n_y");
    if (n_lands < 1) stop("n_lands must be >= 1");
    if (n_lands >= 1e6) stop("n_lands must be < 1e6");
    thread_check(n_threads); // Check that # threads isn't too high

    arma::mat wt_mat(2, 2, arma::fill::none);
    wt_mat(0,0) = wt_vv;
    wt_mat(1,1) = wt_pp;
    wt_mat(0,1) = wt_vp;
    wt_mat(1,0) = wt_pv;

    RcppThread::ProgressBar prog_bar(n_lands * (n_virus + n_pseudo), 1);

    // Extract starting positions, if provided
    // (virus_pseudo_xy0 is empty if *_starts are NULL):
    arma::umat virus_pseudo_xy0 = get_xy_starts(virus_starts, pseudo_starts,
                                                n_virus, n_pseudo, n_x, n_y);

    std::vector<OnePlantTypeSimmer> simmers;
    simmers.reserve(n_lands);
    for (uint32 i = 0; i < n_lands; i++) {
        simmers.push_back(OnePlantTypeSimmer(wt_mat, n_virus, n_pseudo, n_x, n_y,
                                             virus_pseudo_xy0));
    }


    if (n_threads > 1U && n_lands > 1U) {
        auto job = [&] (OnePlantTypeSimmer& simmer) {
            simmer.run(prog_bar, show_progress);
        };
        RcppThread::parallelForEach(simmers, job, n_threads);
    } else {
        for (uint32 i = 0; i < n_lands; i++) {
            simmers[i].run(prog_bar, show_progress);
        }
    }

    arma::icube out(n_x, n_y, n_lands);
    for (uint32 i = 0; i < n_lands; i++) {
        simmers[i].fill_output(out, i);
    }

    return out;

}



//' Convert a landscape cube to a list of dataframes.
//'
//' @param land_cube Cube containing landscape info, output from [sim_plant_types()].
//'
//' @return A list of length `dim(land_cube)[3]`, where each item is a
//'     [`tibble`][tibble::tbl_df] with columns
//'     `x`, `y`, and `type`.
//'
//[[Rcpp::export]]
List land_cube2list(const arma::icube& land_cube) {

    uint32 n_x = land_cube.n_rows;
    uint32 n_y = land_cube.n_cols;
    uint32 n_lands = land_cube.n_slices;

    if (n_x < 2) stop("nrow(land_cube) must be >= 2");
    if (n_x > 1e6) stop("nrow(land_cube) must be <= 1e6");
    if (n_y < 2) stop("ncol(land_cube) must be >= 2");
    if (n_y > 1e6) stop("ncol(land_cube) must be <= 1e6");
    if (n_lands < 1) stop("dim(land_cube)[3] must be >= 1");
    if (n_lands >= 1e6) stop("dim(land_cube)[3] must be < 1e6");

    std::vector<uint32> x;
    std::vector<uint32> y;
    std::vector<uint32> type;
    x.reserve(n_x * n_y);
    y.reserve(n_x * n_y);
    type.reserve(n_x * n_y);

    List out(n_lands);

    DataFrame out_k;

    // Do first landscape while filling x and y (and setting size for type):
    for (uint32 j = 0; j < n_y; j++) {
        for (uint32 i = 0; i < n_x; i++) {
            x.push_back(i+1);
            y.push_back(j+1);
            type.push_back(land_cube(i,j,0));
        }
    }
    out_k = DataFrame::create(_["x"] = x, _["y"] = y, _["type"] = type);
    out_k.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});
    out[0] = out_k;

    if (n_lands > 1) {
        for (uint32 k = 1; k < n_lands; k++) {
            uint32 t = 0;
            for (uint32 j = 0; j < n_y; j++) {
                for (uint32 i = 0; i < n_x; i++) {
                    type[t] = land_cube(i,j,k);
                    t++;
                }
            }
            out_k = DataFrame::create(_["x"] = x, _["y"] = y,
                                       _["type"] = type);
            out_k.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});
            out[k] = out_k;
        }
    }

    return out;

}
