#ifndef __HLIB_ALG_HH
#define __HLIB_ALG_HH
//
// Project     : HLib
// File        : hlib-alg.hh
// Description : HLIBpro include file containing algebra related headers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//
//
//! \defgroup Algebra_Module Algebra
//!
//! This module provides most higher level algebra functions, e.g.
//! matrix multiplication, inversion and factorisation. See also
//! \ref TUTORAlgebraBasic and \ref TUTORMatrixFac for an introduction 
//! into \mcH-arithmetic.
//!
//! To include all algebra functions and classes add
//! \code
//! #include <hlib-alg.hh>
//! \endcode
//! to your source files.
//!

#include "hpro/algebra/approx.hh"
#include "hpro/algebra/hca_func.hh"
#include "hpro/algebra/diag_scale.hh"
#include "hpro/algebra/eval_diag.hh"
#include "hpro/algebra/eval_tri.hh"
#include "hpro/algebra/mat_add.hh"
#include "hpro/algebra/mat_conv.hh"
#include "hpro/algebra/mat_inv.hh"
#include "hpro/algebra/mat_mul.hh"
#include "hpro/algebra/mat_mul_core.hh"
#include "hpro/algebra/mat_mul_diag.hh"
#include "hpro/algebra/mat_norm.hh"
#include "hpro/algebra/mat_fac.hh"
#include "hpro/algebra/mul_vec.hh"
#include "hpro/algebra/solve_diag.hh"
#include "hpro/algebra/solve_tri.hh"
#include "hpro/algebra/TLowRankApx.hh"
#include "hpro/algebra/fft.hh"
#include "hpro/algebra/TCoarsen.hh"

#endif  // __HLIB_ALG_HH
