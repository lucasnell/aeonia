# ifndef __AEONIA_TYPES_H
# define __AEONIA_TYPES_H


/*
 ********************************************************
 Basic integer types used throughout
 ********************************************************
 */

#include <RcppArmadillo.h>
#include <cstdint>

/*
 Armadillo docs:
 > the default width [of integers] is 32 bits when using Armadillo in the R
 > environment (via RcppArmadillo) on either 32-bit or 64-bit platforms
*/

typedef uint_fast8_t uint8;
typedef arma::uword uint32;
typedef arma::sword int32;
typedef uint_fast64_t uint64;
typedef int_fast64_t sint64;



#endif
