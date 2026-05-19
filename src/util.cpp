
#include "aeonia_types.hpp"

#include <RcppArmadillo.h>
#include <pcg/pcg_random.hpp>   // pcg prng


#include "pcg.hpp"              // runif_01 fxn


using namespace Rcpp;




//' Fast nonparametric bootstrapping
//'
//' @param x Vector to bootstrap
//' @param B Number of bootstrapping replicates. Must be `>1` and `< 1e9`.
//'     Defaults to `2000L`.
//' @param alpha Confidence level is `100 * (1 - alpha)` percent.
//'     Defaults to `0.05`.
//'
//' @export
//'
//[[Rcpp::export]]
NumericVector booter(const arma::vec& x,
                     const int32& B = 2000,
                     const double& alpha = 0.05) {

    if (B < (int32)1) stop("B must be > 1");
    if (B > (int32)1000000000) stop("B must be < 10^9");
    if (alpha <= 0.0 || alpha >= 1.0) stop("alpha must be > 0 and < 1");

    uint32 n = x.n_elem;
    if (n < 2U) stop("x must contain at least two elements");

    pcg32 eng;
    seed_pcg(eng);

    uint32 uB = B;

    arma::vec means(uB, arma::fill::none);
    double n_dbl = static_cast<double>(n);

    uint32 k;
    for (uint32 i = 0; i < uB; i++) {
        means.at(i) = 0;
        for (uint32 j = 0; j < n; j++) {
            k = runif_01(eng) * n;
            means.at(i) += x.at(k);
        }
        means.at(i) /= n_dbl;
    }

    arma::vec quants = { alpha/2, 0.5, 1 - alpha/2 };

    arma::vec ci = arma::quantile(means, quants);

    NumericVector out = NumericVector::create(
        _["Lower"] = ci.at(0),
        _["Median"] = ci.at(1),
        _["Upper"] = ci.at(2));

    return out;
}


