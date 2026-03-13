#include <RcppArmadillo.h>      // arma namespace
#include <vector>               // vector class
#include "aeonia_types.hpp"     // integer types
#include "insects.hpp"          // insect types



using namespace Rcpp;



//[[Rcpp::export]]
SEXP make_insect_ptr_cpp(SEXP aphids_xptr, SEXP wasps_xptr, SEXP mummies_xptr) {

    XPtr<AphidPop> aphids_ptr(aphids_xptr);
    const AphidPop& aphids(*aphids_ptr);

    XPtr<WaspPop> wasps_ptr(wasps_xptr);
    const WaspPop& wasps(*wasps_ptr);

    XPtr<MummyPop> mummies_ptr(mummies_xptr);
    const MummyPop& mummies(*mummies_ptr);

    XPtr<InsectPops> insect_xptr(new InsectPops(aphids, wasps, mummies), true);

    return insect_xptr;

}
