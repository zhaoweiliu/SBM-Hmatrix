#ifndef __HLIB_ALGEBRA_TYPES_HH
#define __HLIB_ALGEBRA_TYPES_HH
//
// Project     : HLib
// File        : types.hh
// Description : basic types for algebra module
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/blas/types.hh"

namespace HLIB
{

//! \enum   eval_type_t
//! \brief  defines whether point or block wise evaluation should be performed
enum eval_type_t
{
    point_wise   = 'P',
    block_wise   = 'B'
};

// import evaluation side type from BLAS
using BLAS::eval_side_t;
using BLAS::from_left;
using BLAS::from_right;

// import diagonal type from BLAS
using BLAS::diag_type_t;
using BLAS::unit_diag;
using BLAS::general_diag;

//! \enum   storage_type_t
//! \brief  options for who diagonal blocks are stored
enum storage_type_t
{
    store_normal  = 'N',   // blocks contain normal matrices
    store_inverse = 'I'    // blocks contain inverse matrices
};

//! \enum   eager_lazy_t
//! \brief  defines eager of lazy evaluation
enum eager_lazy_t : bool
{
    eval_eager = false,    // eager evaluation
    eval_lazy  = true      // lazy evaluation
};

}// namespace HLIB

#endif  // __HLIB_ALGEBRA_TYPES_HH
