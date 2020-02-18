#ifndef __HLIB_MAT_MUL_DIAG_HH
#define __HLIB_MAT_MUL_DIAG_HH
//
// Project     : HLib
// File        : mat_mul_diag.hh
// Description : matrix multiplication with diagonal matrices (inplace)
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/base/TTruncAcc.hh"
#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TScalarVector.hh"
#include "hpro/misc/TProgressBar.hh"

namespace HLIB
{


//!
//! \ingroup Algebra_Module
//! \brief   compute C ≔ β·C + α·op(A)·op(D)·op(B) with (block) diagonal matrix D
//!
template < typename value_t >
void
multiply_diag ( const value_t  alpha,
                const matop_t  op_A, const TMatrix *  A,
                const matop_t  op_D, const TMatrix *  D,
                const matop_t  op_B, const TMatrix *  B,
                const value_t  beta, TMatrix *        C,
                const TTruncAcc & acc,
                TProgressBar * progress = NULL );

//!
//! \ingroup Algebra_Module
//! \brief   compute C ≔ C + α·op(A)·op(D)·op(B) with (block) diagonal matrix D
//!          using accumulators
//!
template < typename value_t >
void
multiply_diag_accu ( const value_t  alpha,
                     const matop_t  op_A, const TMatrix *  A,
                     const matop_t  op_D, const TMatrix *  D,
                     const matop_t  op_B, const TMatrix *  B,
                     TMatrix *      C,
                     const TTruncAcc & acc );

//!
//! \ingroup Algebra_Module
//! \brief   compute C ≔ α·op(A)·op(D)·op(B) with (block) diagonal matrix D
//!          - not all matrices must be blocked!
//!
template < typename value_t >
std::unique_ptr< TMatrix >
multiply_diag ( const value_t  alpha,
                const matop_t  op_A, const TMatrix *  A,
                const matop_t  op_D, const TMatrix *  D,
                const matop_t  op_B, const TMatrix *  B );

//!
//! \ingroup Algebra_Module
//! \brief   return number of steps for computing C ≔ C + op(A)·op(D)·op(B)
//!
size_t
multiply_diag_steps ( const matop_t  op_A, const TMatrix *  A,
                      const matop_t  op_D, const TMatrix *  D,
                      const matop_t  op_B, const TMatrix *  B,
                      const TMatrix *  C );

//!
//! \ingroup Algebra_Module
//! \brief   compute B = diag(v)·A and overwrite A
//!
void mul_diag_left  ( const TScalarVector & v,
                      TMatrix *             A );

//!
//! \ingroup Algebra_Module
//! \brief   compute B = A·diag(v) and overwrite A
//!
void mul_diag_right ( TMatrix *             A,
                      const TScalarVector & v );

//!
//! start accumulator based matrix multiplication C = C + α·A·D·B
//! with diagonal D
//!
template < typename value_t >
void
add_diag_product ( const value_t      alpha,
                   const matop_t      op_A,
                   const TMatrix *    A,
                   const matop_t      op_D,
                   const TMatrix *    D,
                   const matop_t      op_B,
                   const TMatrix *    B,
                   TMatrix *          C,
                   const TTruncAcc &  acc,
                   const bool         lazy = CFG::Arith::lazy_eval );

}// namespace HLIB

#endif  // __HLIB_MAT_MUL_DIAG_HH
