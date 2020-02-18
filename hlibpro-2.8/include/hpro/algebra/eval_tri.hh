#ifndef __HLIB_EVAL_TRI_HH
#define __HLIB_EVAL_TRI_HH
//
// Project     : HLib
// File        : eval_tri.hh
// Description : evaluate block triangular matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TVector.hh"
#include "hpro/algebra/solve_types.hh"

namespace HLIB
{

/////////////////////////////////////////////////////////////////////////////////////
//
// evaluate lower triangular systems
//
/////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup Algebra_Module
//! \brief   evaluate L·U·x=y with lower triangular L and upper triangular U
//!          - L and U are both stored in A
//!          - on entry: v = x
//!          - on exit : v = y
//!
void
eval       ( const TMatrix *        A,
             TVector *              v,
             const matop_t          op,
             const eval_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \brief   evaluate A·x=y with lower triangular A
//!          - on entry: v = x
//!          - on exit : v = y
//!
void
eval_lower ( const TMatrix *        A,
             TVector *              v,
             const matop_t          op,
             const eval_option_t &  options );

//!
//! \ingroup Algebra_Module
//! \brief   evaluate A·x = y with upper triangular A
//!          - on entry: v = x
//!          - on exit : v = y
//!
void
eval_upper ( const TMatrix *        A,
             TVector *              v,
             const matop_t          op,
             const eval_option_t &  options );

}// namespace HLIB

#endif  // __HLIB_EVAL_TRI_HH
