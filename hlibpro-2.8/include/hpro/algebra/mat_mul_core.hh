#ifndef __HLIB_MAT_MUL_CORE_HH
#define __HLIB_MAT_MUL_CORE_HH
//
// Project     : HLib
// File        : mat_mul_cure.hh
// Description : core matrix multiplication functions for base matrix types
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/structure.hh"

namespace HLIB
{

//
// multiplication routines for leaf matrix types
//
template < typename T >
std::unique_ptr< TMatrix >
multiply  ( const T               alpha,
            const matop_t         op_A,
            const TDenseMatrix *  A,
            const matop_t         op_B,
            const TDenseMatrix *  B );

template < typename T >
std::unique_ptr< TMatrix >
multiply  ( const T               alpha,
            const matop_t         op_A,
            const TDenseMatrix *  A,
            const matop_t         op_B,
            const TRkMatrix *     B );

template < typename T >
std::unique_ptr< TMatrix >
multiply  ( const T               alpha,
            const matop_t         op_A,
            const TRkMatrix *     A,
            const matop_t         op_B,
            const TDenseMatrix *  B );

template < typename T >
std::unique_ptr< TMatrix >
multiply  ( const T               alpha,
            const matop_t         op_A,
            const TRkMatrix *     A,
            const matop_t         op_B,
            const TRkMatrix *     B );

//
// test function for this module
//

namespace MAT_MUL_CORE
{

void
test ();

}// namespace MAT_MUL_CORE

}// namespace HLIB

#endif // __HLIB_MAT_MUL_CORE_HH
