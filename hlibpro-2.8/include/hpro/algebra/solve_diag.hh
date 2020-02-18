#ifndef __HLIB_SOLVE_DIAG_HH
#define __HLIB_SOLVE_DIAG_HH
//
//! \file        solve_diag.hh
//
// Project     : HLib
// Description : solve block triangular matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/base/TTruncAcc.hh"
#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TVector.hh"
#include "hpro/algebra/solve_tri.hh"

namespace HLIB
{

/////////////////////////////////////////////////////////////////////////////////////
//
// solve upper triangular systems
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \fn      solve_diag_right
//! \brief   solve B D = A with diagonal matrix D and A known
//!          - B overwrites A
//!          - if D_inv != NULL, D_inv must contain inverted diagonal of D
//!
void
solve_diag_right  ( TMatrix *              A,
                    const matop_t          op_D,
                    const TMatrix *        D,
                    const TMatrix *        D_inv,
                    const TTruncAcc &      acc,
                    const solve_option_t & options,
                    TProgressBar *         progress  = NULL );

inline
void
solve_diag_right  ( TMatrix *              A,
                    const matop_t          op_D,
                    const TMatrix *        D,
                    const TTruncAcc &      acc,
                    const solve_option_t & options,
                    TProgressBar *         progress  = NULL )
{
    solve_diag_right( A, op_D, D, NULL, acc, options, progress );
}

//!
//! \ingroup Algebra_Module
//! \fn      solve_diag_right_steps
//! \brief   return number of steps of solve_diag_right for progress meter
//!
size_t
solve_diag_right_steps  ( const TMatrix *        A,
                          const matop_t          op_D,
                          const TMatrix *        D,
                          const solve_option_t & options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_diag_left
//! \brief   solve D B = A with diagonal matrix D and A known
//!          - B overwrites A
//!          - if D_inv != NULL, D_inv must contain inverted diagonal of D
//!
void
solve_diag_left   ( const matop_t          op_D,
                    const TMatrix *        D,
                    const TMatrix *        D_inv,
                    TMatrix *              A,
                    const TTruncAcc &      acc,
                    const solve_option_t & options,
                    TProgressBar *         progress  = NULL );

inline
void
solve_diag_left   ( const matop_t          op_D,
                    const TMatrix *        D,
                    TMatrix *              A,
                    const TTruncAcc &      acc,
                    const solve_option_t & options,
                    TProgressBar *         progress  = NULL )
{
    solve_diag_left( op_D, D, NULL, A, acc, options, progress );
}

//!
//! \ingroup Algebra_Module
//! \fn      solve_diag_left_steps
//! \brief   return number of steps of solve_diag_left for progress meter
//!
size_t
solve_diag_left_steps  ( const matop_t          op_D,
                         const TMatrix *        D,
                         const TMatrix *        A,
                         const solve_option_t & options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_upper
//! \brief   solve D x = y with diagonal D
//!
void
solve_diag        ( const matop_t          op_D,
                    const TMatrix *        D,
                    TVector *              v,
                    const solve_option_t & options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_upper
//! \brief   solve D x = y with diagonal D
//!          - if D_inv ≠ NULL, D_inv must contain the inverted diagonal blocks of D
//!
void
solve_diag        ( const matop_t          op_D,
                    const TMatrix *        D,
                    const TMatrix *        D_inv,
                    TVector *              v,
                    const solve_option_t & options );

}// namespace HLIB

#endif  // __HLIB_SOLVE_DIAG_HH
