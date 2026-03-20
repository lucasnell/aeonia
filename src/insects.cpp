
#include "insects.hpp"          // insect types



using namespace Rcpp;



//[[Rcpp::export]]
SEXP make_insects_ptr_cpp(SEXP aphids_xptr, SEXP wasps_xptr, SEXP mummies_xptr) {

    XPtr<AphidPop> aphids_ptr(aphids_xptr);
    const AphidPop& aphids(*aphids_ptr);

    XPtr<WaspPop> wasps_ptr(wasps_xptr);
    const WaspPop& wasps(*wasps_ptr);

    XPtr<MummyPop> mummies_ptr(mummies_xptr);
    const MummyPop& mummies(*mummies_ptr);

    XPtr<InsectPops> insects_xptr(new InsectPops(aphids, wasps, mummies), true);

    return insects_xptr;

}

// Used for testing
//[[Rcpp::export]]
arma::vec get_aphid_X0(SEXP insects_ptr) {
    XPtr<InsectPops> insects_xptr(insects_ptr);
    const AphidPop& aphids(insects_xptr->aphids);
    arma::vec X_0 = aphids.unparas_X_0();
    return X_0;
}






//' Test population dynamics for insects for a set of parameters.
//'
//' This function is just for one plant, so the probability that an
//' alate leaves the plant (argument `fly_p` in [make_insects_ptr()]) is always 1.
//'
//' @inheritParams sim_plantscape
//'
//' @export
//'
//[[Rcpp::export]]
DataFrame test_insect_pops(const uint32& max_t,
                           const double& N0,
                           const double& W0,
                           const double& M0,
                           const double& Y0,
                           SEXP insects_ptr) {

    XPtr<InsectPops> insects_xptr(insects_ptr);
    // Copy here so that changes I'll make won't propagate:
    InsectPops insects = *insects_xptr;

    insects.aphids.fly_p = 0;
    insects.aphids.refresh_abunds(N0, W0);

    insects.mummies.refresh_abunds(M0);

    insects.wasps.zeta = 0;
    insects.wasps.refresh_abunds(Y0);


    std::vector<uint32> time;
    std::vector<double> aphids;
    std::vector<double> alates;
    std::vector<double> parasitized;
    std::vector<double> mummies;
    std::vector<double> wasps;
    time.reserve(max_t+1);
    aphids.reserve(max_t+1);
    alates.reserve(max_t+1);
    parasitized.reserve(max_t+1);
    mummies.reserve(max_t+1);
    wasps.reserve(max_t+1);
    time.push_back(0);
    aphids.push_back(N0);
    alates.push_back(W0);
    parasitized.push_back(0);
    mummies.push_back(M0);
    wasps.push_back(Y0);

    pcg32 eng;
    seed_pcg(eng);

    for (uint32 t = 0; t < max_t; t++) {
        insects.iterate(eng);
        time.push_back(t+1);
        aphids.push_back(insects.aphids.apterous.total());
        alates.push_back(insects.aphids.alates.total());
        parasitized.push_back(insects.aphids.paras.total());
        mummies.push_back(insects.mummies.total());
        wasps.push_back(insects.wasps.Y);
    }

    DataFrame out_df = DataFrame::create(_["time"] = time,
                                         _["aphids"] = aphids,
                                         _["alates"] = alates,
                                         _["parasitized"] = parasitized,
                                         _["mummies"] = mummies,
                                         _["wasps"] = wasps);

    out_df.attr("class") = CharacterVector({"tbl_df", "tbl", "data.frame"});

    return out_df;

}
