#ifndef __AEONIA_UTIL_H
#define __AEONIA_UTIL_H

#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <RcppThread.h>         // multithreading

#include "aeonia_types.hpp"


using namespace Rcpp;


// Get # of cores from options("mc.cores"), or 1 if it's not set.
inline void assign_mc_cores(uint32& cores) {

    Rcpp::Environment base("package:base");
    Rcpp::Function get_option = base["getOption"];
    Rcpp::Nullable<Rcpp::IntegerVector> mc_cores = get_option("mc.cores");

    if (mc_cores != R_NilValue) {
        cores = as<uint32>(mc_cores);
    } else cores = 1;

    return;
}


/*
 Check that the number of threads doesn't exceed the number available.
 If `n_threads == 0`, use `options("mc.cores")`, which will assign 1 if that
 option isn't assigned.
 */
inline void thread_check(uint32& n_threads) {

    if (n_threads == 0) assign_mc_cores(n_threads);

    uint32 max_threads = std::thread::hardware_concurrency();

    if (n_threads > max_threads) {
        std::string mt_str = std::to_string(max_threads);
        std::string err_msg = "\nThe number of requested threads (" +
            std::to_string(n_threads) +
            ") exceeds the max available on the system (" + mt_str + ").";
        stop(err_msg.c_str());
    }

    return;
}


// Retrieve dataset named `dataset` of type `T` from this package
template <typename T>
T retrieve_dataset(const std::string& dataset) {

    Environment pkg("package:aeonia");
    T data = pkg[dataset];
    return data;
}





#endif
