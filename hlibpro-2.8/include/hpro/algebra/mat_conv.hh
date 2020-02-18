#ifndef __HLIB_MAT_CONV_HH
#define __HLIB_MAT_CONV_HH
//
// Project     : HLib
// File        : mat_conv.hh
// Description : functions for converting matrices into various formats
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TBlockMatrix.hh"
#include "hpro/matrix/TSparseMatrix.hh"
#include "hpro/cluster/TClusterBasis.hh"

namespace HLIB
{

//!
//! Convert \a A into existing (!) matrix \a C with accuracy \a acc.
//!
void
convert  ( const TMatrix *    A,
           TMatrix *          C,
           const TTruncAcc &  acc );

//!
//! Convert \a A into dense matrix.
//!
std::unique_ptr< TDenseMatrix >
to_dense ( const TMatrix *    A );

//!
//! Convert \a A into existing (!) dense matrix \a C.
//!
void
to_dense ( const TMatrix *    A,
           TDenseMatrix *     C );

//!
//! Convert all leaf blocks in \a A into dense format
//!
void
to_dense_blocks ( TBlockMatrix * A );

//!
//! Convert \a A into low-rank matrix.
//!
std::unique_ptr< TRkMatrix >
to_rank  ( const TMatrix *    A,
           const TTruncAcc &  acc );

//!
//! Convert \a A into existing (!) low-rank matrix \a C.
//!
void
to_rank  ( const TMatrix *    A,
           TRkMatrix *        C,
           const TTruncAcc &  acc );

//!
//! Convert \a A into existing (!) block matrix \a C.
//!
void
to_block ( const TMatrix *    A,
           TBlockMatrix *     C,
           const TTruncAcc &  acc );
    
//!
//! Convert \a A into H²-matrix
//!
std::unique_ptr< TMatrix >
to_h2    ( const TMatrix *    A,
           const TTruncAcc &  acc );

//!
//! Convert (project) \a A into H²-matrix
//! (accuracy is determined by cluster basis)
//!
template <typename T>
std::unique_ptr< TMatrix >
to_h2    ( const TMatrix *             A,
           const TClusterBasis< T > *  rowcb,
           const TClusterBasis< T > *  colcb );
    
//!
//! Convert \a A into sparse matrix
//!
std::unique_ptr< TSparseMatrix >
to_sparse ( const TMatrix *   A );

//!
//! Convert \a A into sparse matrix C
//!
void
to_sparse ( const TMatrix *   A,
            TSparseMatrix *   C );

}// namespace HLIB

#endif  // __HLIB_MAT_CONV_HH
