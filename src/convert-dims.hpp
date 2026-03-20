#ifndef __AEONIA_CONVERT_DIMS_H
#define __AEONIA_CONVERT_DIMS_H

#include "aeonia_types.hpp"

#include <RcppArmadillo.h>
#include <vector>


// To avoid many warnings from BOOST
#pragma clang diagnostic ignored "-Wlanguage-extension-token"
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
#include "boost/multi_array.hpp"
#pragma clang diagnostic warning "-Wlanguage-extension-token"
#pragma clang diagnostic warning "-Wdeprecated-declarations"



using namespace Rcpp;


template <class C>
using vMatrix = boost::multi_array<C, 2>;

typedef std::vector<uint32> vMatSize;



// Simple struct for storing x and y coordinates
struct XY {
    uint32 x;
    uint32 y;
    XY(const uint32& x_, const uint32& y_) : x(x_), y(y_) {};
};



// Extract the `k`th bit of `n` and cast to bool `b`:
inline bool get_bit_bool(const uint32& k, const uint32& n) {
    bool b = (n & ( 1 << k )) >> k;
    return b;
}
// Same but cast as int
inline int get_bit_int(const uint32& k, const uint32& n) {
    int b = (n & ( 1 << k )) >> k;
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

    inline DimensionConverter(const uint32& n_x_, const uint32& n_y_)
        : neigh_x(), neigh_y(), n_x(n_x_), n_y(n_y_) {
        neigh_x.reserve(3);
        neigh_y.reserve(3);
    }

    // Convert from 1D to 2D:
    inline void to_2d(uint32& x, uint32& y, const uint32& k) const {
        x = k - n_x * (k / n_x);
        y = k / n_x;
        return;
    }
    // Overloaded for signed ints (for use with Rcpp::IntegerVector)
    inline void to_2d(int& x, int& y, const uint32& k) const {
        x = k - n_x * (k / n_x);
        y = k / n_x;
        return;
    }
    // Convert from 2D to 1D:
    inline void to_1d(uint32& k, const uint32& x, const uint32& y) const {
        k = (y * n_x + x);
        return;
    }
    // Overloaded to return instead of assign
    inline uint32 to_1d(const uint32& x, const uint32& y) const {
        uint32 k = (y * n_x + x);
        return k;
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





// // DIY Matrix class for any class:
// template <class C>
// class vMatrix {
//
//     typedef boost::multi_array<C, 2> array_type;
//     typedef array_type::index index;
//
//     array_type data;
//
//     uint32 n_x, n_y;
//
//
// public:
//
//     // // =============================================
//     // // Used so that vMatrix(0,0) = true works.
//     // // from https://stackoverflow.com/a/51503646/5016095
//     // using value_type = C;
//     // class Proxy {
//     // public:
//     //     friend class vMatrix;
//     //     operator value_type() const {
//     //         return v_matrix.get(m_index);
//     //     }
//     //     Proxy& operator=(value_type value) {
//     //         v_matrix.set(m_index, value);
//     //         return *this;
//     //     }
//     //     // value_type operator->() { return v_matrix.get(m_index); }
//     // private:
//     //     Proxy(vMatrix& v_matrix_, const uint32& index)
//     //         : v_matrix(v_matrix_), m_index(index) {}
//     //     vMatrix& v_matrix;
//     //     uint32 m_index;
//     // };
//     //
//     // value_type operator()(const uint32& x, const uint32& y) const {
//     //     uint32 k = y * n_x + x;
//     //     return get(k);
//     // }
//     // Proxy operator()(const uint32& x, const uint32& y) {
//     //     uint32 k = y * n_x + x;
//     //     return Proxy(*this, k);
//     // }
//     // value_type get(const uint32& k) const {
//     //     if (k >= data.size())
//     //         throw std::runtime_error("INTERNAL ERROR: k beyond bounds!");
//     //     return data[k];
//     // }
//     // void set(const uint32& k, value_type value) {
//     //     if (k >= data.size())
//     //         throw std::runtime_error("INTERNAL ERROR: k beyond bounds!");
//     //     data[k] = value;
//     //     return;
//     // }
//     // // =============================================
//
//     // Constructors
//     vMatrix() : n_x(0), n_y(0), data() {}
//     vMatrix(const uint32& n_x_, const uint32& n_y_)
//         : n_x(n_x_), n_y(n_y_), data(boost::extents[n_x][n_y]) {}
//     vMatrix(const uint32& n_x_, const uint32& n_y_, const C& item)
//         : n_x(n_x_), n_y(n_y_), data(n_x_ * n_y_, item) {}
//
//     void set_size(const uint32& n_x_, const uint32& n_y_) {
//         n_x = n_x_;
//         n_y = n_y_;
//         data.resize(n_x_ * n_y_);
//         return;
//     }
//     void set_size(const uint32& n_x_,
//                   const uint32& n_y_,
//                   const C& item) {
//         n_x = n_x_;
//         n_y = n_y_;
//         data.reserve(n_x_ * n_y_);
//         for (uint32 i = 0; i < (n_x_ * n_y_); i++) data.push_back(item);
//         return;
//     }
//
//     // To fill like a vector:
//     void reserve(const uint32& n_x_,
//                  const uint32& n_y_) {
//         n_x = n_x_;
//         n_y = n_y_;
//         data.reserve(n_x_ * n_y_);
//     }
//     // Be careful using the function below because it doesn't adjust n_x and n_y
//     // Only use it if you'll set n_x and n_y to the appropriate values elsewhere
//     void push_back(const C& val) {
//         data.push_back(val);
//     }
//
//     // // Access operators (read/write)
//     // C& operator()(const uint32& x, const uint32& y) {
//     //     if (y * n_x + x >= data.size())
//     //         throw std::runtime_error("INTERNAL ERROR: x and y beyond bounds!");
//     //     return data[y * n_x + x];
//     // }
//     // // Access operator (read-only)
//     // const C& operator()(const uint32& x, const uint32& y) const {
//     //     if (y * n_x + x >= data.size())
//     //         throw std::runtime_error("INTERNAL ERROR: x and y beyond bounds!");
//     //     return data[y * n_x + x];
//     // }
//     // // Same but for a single 1D index:
//     // C& operator()(const uint32& k) {
//     //     if (k >= data.size())
//     //         throw std::runtime_error("INTERNAL ERROR: k beyond bounds!");
//     //     return data[k];
//     // }
//     // const C& operator()(const uint32& k) const {
//     //     if (k >= data.size())
//     //         throw std::runtime_error("INTERNAL ERROR: k beyond bounds!");
//     //     return data[k];
//     // }
//
//     // Member functions to get dimensions
//     uint32 n_rows() const noexcept { return n_x; }
//     uint32 n_cols() const noexcept { return n_y; }
//     uint32 size() const noexcept { return data.size(); }
//
// };






#endif
