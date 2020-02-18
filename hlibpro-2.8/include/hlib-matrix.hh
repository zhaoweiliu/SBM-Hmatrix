#ifndef __HLIB_MATRIX_HH
#define __HLIB_MATRIX_HH
//
// Project     : HLib
// File        : hlib-matrix.hh
// Description : HLIBpro include file containing matrix related headers
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//
//! \defgroup Matrix_Module Matrix
//!
//! This modules provides all high level matrix classes, e.g.
//! dense, low rank and block matrices.
//!
//! \code
//! #include <hlib-matrix.hh>
//! \endcode
//!

#include <hpro/matrix/TMatrix.hh>
#include <hpro/matrix/TBSHMBuilder.hh>
#include <hpro/matrix/TBlockMatrix.hh>
#include <hpro/matrix/TDenseMatrix.hh>
#include <hpro/matrix/TRkMatrix.hh>
#include <hpro/matrix/THMatrix.hh>
#include <hpro/matrix/TH2Matrix.hh>
#include <hpro/matrix/TSparseMatrix.hh>
#include <hpro/matrix/TDiagMatrix.hh>
#include <hpro/matrix/TFacMatrix.hh>
#include <hpro/matrix/TFacInvMatrix.hh>
#include <hpro/matrix/TPermMatrix.hh>
#include <hpro/matrix/TMatBuilder.hh>
#include <hpro/matrix/TMatrixHierarchy.hh>
#include <hpro/matrix/TMatrixProduct.hh>
#include <hpro/matrix/TMatrixSum.hh>
#include <hpro/matrix/TJacobi.hh>
#include <hpro/matrix/TSOR.hh>
#include <hpro/matrix/TUniformMatrix.hh>
#include <hpro/matrix/BlockDiag.hh>
#include <hpro/matrix/NearField.hh>
#include <hpro/matrix/structure.hh>

#endif  // __HLIB_MATRIX_HH
