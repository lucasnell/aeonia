#ifndef __AEONIA_CONVERT_DIMS_H
#define __AEONIA_CONVERT_DIMS_H

#include <RcppArmadillo.h>
#include <vector>

#include "aeonia_types.hpp"


using namespace Rcpp;



// Simple struct for storing x and y coordinates
struct XY {
    uint32 x;
    uint32 y;
};



// Extract the `k`th bit of `n` and cast to bool `b`:
inline bool get_bit_bool(const uint32& k, const uint32& n) {
    bool b = (n & ( 1 << k )) >> k;
    return b;
}

// Set the `k`th bit of `n` to `1`:
inline void bit_set1(const uint32& k, uint32& n) {
    n |= ((uint32)1 << k);
    return;
}






/*
 Class to convert back and forth between 2D to 1D indices.
 It assumes that x and y coordinates are sorted first by y coordinates, then
 by x coordinates.
 It also assumes 0-based indices.
 An example matrix of x and y coordinates might start like this:

 #>      x y
 #> [0,] 0 0
 #> [1,] 1 0
 #> [2,] 2 0
 #> [3,] 3 0
 #> [4,] 4 0
 #> [5,] 5 0

 */
class DimensionConverter {

    std::vector<uint32> neigh_x;
    std::vector<uint32> neigh_y;

    uint32 n_x;
    uint32 n_y;

public:

    DimensionConverter(const uint32& n_x_, const uint32& n_y_)
        : neigh_x(), neigh_y(), n_x(n_x_), n_y(n_y_) {
        neigh_x.reserve(3);
        neigh_y.reserve(3);
    }

    // Convert from 1D to 2D:
    void to_2d(uint32& x, uint32& y, const uint32& k) const {
        x = k - n_y * (k / n_y);
        y = k / n_y;
        return;
    }
    // Overloaded for signed ints (for use with Rcpp::IntegerVector)
    void to_2d(int& x, int& y, const uint32& k) const {
        x = k - n_y * (k / n_y);
        y = k / n_y;
        return;
    }
    // Convert from 2D to 1D:
    void to_1d(uint32& k, const uint32& x, const uint32& y) const {
        k = (y * n_x + x);
        return;
    }
    /*
     Fill indices (in 1D) for all "next door" neighbors based on a 1D
     input coordinate.
     NOTE:
        - It clears `indices` before adding to it
        - It also returns an index for the focal point
     */
    void nextdoor_neighbors(std::vector<uint32>& indices,
                            const uint32& k) {
        uint32 x0 = k - n_y * (k / n_y);
        uint32 y0 = k / n_y;
        indices.clear();
        neigh_x.clear();
        neigh_y.clear();
        if (x0 > 0) neigh_x.push_back(x0-1);
        neigh_x.push_back(x0);
        if (x0 < n_x-1) neigh_x.push_back(x0+1);
        if (y0 > 0) neigh_y.push_back(y0-1);
        neigh_y.push_back(y0);
        if (y0 < n_x-1) neigh_y.push_back(y0+1);
        uint32 k_out;
        for (const uint32& x : neigh_x) {
            for (const uint32& y : neigh_y) {
                to_1d(k_out, x, y);
                indices.push_back(k_out);
            }
        }
        return;
    }


};






#endif
