#ifndef __HLIB_SOLVE_TRI_HH
#define __HLIB_SOLVE_TRI_HH
//
// Project     : HLib
// File        : solve_tri.hh
// Description : solve block triangular matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/base/TTruncAcc.hh"
#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TVector.hh"
#include "hpro/misc/TProgressBar.hh"
#include "hpro/algebra/solve_types.hh"

namespace HLIB
{

/////////////////////////////////////////////////////////////////////////////////////
//
// solve lower triangular systems
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \fn      solve_lower_left
//! \brief   solve op(L) B = A with lower triangular matrix L
//!          - B overwrites A
//!          - if D != NULL, D is assumed to contain inverted diagonal of L
//!
void
solve_lower_left  ( const matop_t           op_L,
                    const TMatrix *         L,
                    const TMatrix *         D,
                    TMatrix *               A,
                    const TTruncAcc &       acc,
                    const solve_option_t &  options,
                    TProgressBar *          progress = NULL );

inline
void
solve_lower_left  ( const matop_t           op_L,
                    const TMatrix *         L,
                    TMatrix *               A,
                    const TTruncAcc &       acc,
                    const solve_option_t &  options,
                    TProgressBar *          progress = NULL )
{
    solve_lower_left( op_L, L, NULL, A, acc, options, progress );
}

//!
//! \fn      solve_lower_left_steps
//! \brief   report number of steps in \see solve_lower_left for progress meter initialisation
//!
size_t
solve_lower_left_steps ( const TMatrix *         L,
                         const TMatrix *         A,
                         const solve_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_lower_right
//! \brief   solve B op(L) = A with lower triangular matrix L
//!          - B overwrites A
//!          - if D != NULL, D is assumed to contain inverted diagonal of L
//!
void
solve_lower_right ( TMatrix *               A,
                    const matop_t           op_L,
                    const TMatrix *         L,
                    const TMatrix *         D,
                    const TTruncAcc &       acc,
                    const solve_option_t &  options,
                    TProgressBar *          progress = NULL );

inline
void
solve_lower_right ( TMatrix *               A,
                    const matop_t           op_L,
                    const TMatrix *         L,
                    const TTruncAcc &       acc,
                    const solve_option_t &  options,
                    TProgressBar *          progress = NULL )
{
    solve_lower_right( A, op_L, L, NULL, acc, options, progress );
}

//!
//! \fn     solve_lower_right_steps
//! \brief  report number of steps in \see solve_lower_right for progress meter initialisation
//!
size_t
solve_lower_right_steps ( const matop_t           op_L,
                          const TMatrix *         L,
                          const TMatrix *         A,
                          const solve_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_lower
//! \brief   solve L x = y with lower triangular L
//!
void
solve_lower       ( const matop_t           op_L,
                    const TMatrix *         L,
                    TVector *               v,
                    const solve_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_lower
//! \brief   solve L x = y with lower triangular L
//!          - if D ≠ NULL, it is assumed to contain the inverted diagonal blocks of L
//!
void
solve_lower       ( const matop_t           op_L,
                    const TMatrix *         L,
                    const TMatrix *         D,
                    TVector *               v,
                    const solve_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_diag_lower_right
//! \brief   solve B D op(L) = A with diagonal D and lower triangular L
//!          - B overwrites A
//!          - if D_inv != nullptr then D_inv = D^-1
//!          - D is assume to be general diagonal matrix, e.g. NOT unit diagonal
//!          - D and L are assumed to have same eval mode, e.g. block wise or point wise
//!
void
solve_diag_lower_right ( TMatrix *               A,
                         const matop_t           op_L,
                         const TMatrix *         L,
                         const TMatrix *         D,
                         const TMatrix *         DI,
                         const TTruncAcc &       acc,
                         const solve_option_t &  options,
                         TProgressBar *          progress = NULL );

//!
//! \fn     solve_diag_lower_right_steps
//! \brief  report number of steps in \see solve_diag_lower_right for progress meter initialisation
//!
size_t
solve_diag_lower_right_steps ( const matop_t           op_L,
                               const TMatrix *         L,
                               const TMatrix *         A,
                               const solve_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_lower_diag
//! \brief   solve L D x = y with lower triangular L and diagonal D
//!          - if D_inv ≠ NULL, D_inv = D^-1
//!          - D is assume to be general diagonal matrix, e.g. NOT unit diagonal
//!          - D and L are assumed to have same eval mode, e.g. block wise or point wise
//!
void
solve_lower_diag ( const matop_t          op_L,
                   const TMatrix *        L,
                   const matop_t          op_D,
                   const TMatrix *        D,
                   const TMatrix *        D_inv,
                   TVector *              v,
                   const solve_option_t & options );

/////////////////////////////////////////////////////////////////////////////////////
//
// solve upper triangular systems
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \fn     solve_upper_left
//! \brief  solve op(U) B = A with upper triangular matrix U and A known
//!         - B overwrites A
//!         - if D != NULL, D contains inverted diagonal of U
//!
void
solve_upper_left  ( const matop_t             op_U,
                    const TMatrix *           U,
                    const TMatrix *           D,
                    TMatrix *                 A,
                    const TTruncAcc &         acc,
                    const solve_option_t      options,
                    TProgressBar *            progress = NULL );

inline
void
solve_upper_left  ( const matop_t             op_U,
                    const TMatrix *           U,
                    TMatrix *                 A,
                    const TTruncAcc &         acc,
                    const solve_option_t      options,
                    TProgressBar *            progress = NULL )
{
    solve_upper_left( op_U, U, NULL, A, acc, options, progress );
}

//!
//! \fn     solve_upper_left_steps
//! \brief  report number of steps in \see solve_upper_left for progress meter initialisation
//!
size_t
solve_upper_left_steps  ( const TMatrix *       U,
                          const TMatrix *       A,
                          const solve_option_t  options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_upper_right
//! \brief   solve B U = A with upper triangular matrix U and A known
//!          - B overwrites A
//!          - if D != NULL, D contains inverted diagonal of U
//!
void
solve_upper_right ( TMatrix *               A,
                    const TMatrix *         U,
                    const TMatrix *         D,
                    const TTruncAcc &       acc,
                    const solve_option_t &  options,
                    TProgressBar *          progress = NULL );

inline
void
solve_upper_right ( TMatrix *               A,
                    const TMatrix *         U,
                    const TTruncAcc &       acc,
                    const solve_option_t &  options,
                    TProgressBar *          progress = NULL )
{
    solve_upper_right( A, U, NULL, acc, options, progress );
}

//!
//! \fn     solve_upper_right_steps
//! \brief  report number of steps in \see solve_upper_right for progress meter initialisation
//!
size_t
solve_upper_right_steps ( const TMatrix *         L,
                          const TMatrix *         A,
                          const solve_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_upper
//! \brief   solve A x = y with upper triangular A
//!
void
solve_upper       ( const matop_t           op_U,
                    const TMatrix *         U,
                    TVector *               v,
                    const solve_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \fn      solve_upper
//! \brief   solve A x = y with upper triangular A
//!          - if D ≠ NULL, it is assumed to contain the inverted diagonal blocks of A
//!
void
solve_upper       ( const matop_t           op_U,
                    const TMatrix *         U,
                    const TMatrix *         D,
                    TVector *               v,
                    const solve_option_t &  options );

}// namespace HLIB

#endif  // __HLIB_SOLVE_TRI_HH
