#ifndef __HLIB_EVAL_DIAG_HH
#define __HLIB_EVAL_DIAG_HH
//
// Project     : HLib
// File        : eval_diag.hh
// Description : evaluate block diagonal matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TVector.hh"
#include "hpro/algebra/solve_types.hh"
#include "hpro/algebra/eval_tri.hh"

namespace HLIB
{

//!
//! \ingroup Algebra_Module
//! \brief   evaluate D·x=y with (block) diagonal D
//!          - on entry: v = x
//!          - on exit : v = y
//!
void
eval_diag  ( const TMatrix *        A,
             TVector *              v,
             const matop_t          op,
             const eval_option_t &  options );

}// namespace HLIB

#endif  // __HLIB_EVAL_DIAG_HH
