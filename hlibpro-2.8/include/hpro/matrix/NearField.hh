#ifndef __HLIB_NEARFIELD_HH
#define __HLIB_NEARFIELD_HH
//
// Project     : HLib
// File        : Nearfield.hh
// Description : restrict H-matrices to nearfield part
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TSparseMatrix.hh"

namespace HLIB
{

//!
//! \fn     restrict_nearfield
//! \brief  Restrict given H-Matrix to nearfield part.
//! \param  A                 Matrix to be restricted
//! \param  delete_farfield   If true, farfield matrices (e.g., low-rank) will
//!                           be deleted, otherwise, they are zeroed (rank=0)
//!
void
restrict_nearfield  ( TMatrix *        A,
                      const bool       delete_farfield );

//!
//! \fn     nearfield
//! \brief  Copy nearfield part of given H-Matrix.
//! \param  A                  Matrix, of which the nearfield part is copied.
//! \param  without_farfield   If true, farfield matrices (e.g., low-rank) will not 
//!                            be copied, otherwise, they are copied with rank 0
//!
std::unique_ptr< TMatrix >
nearfield           ( const TMatrix *  A,
                      const bool       without_farfield );

//!
//! \fn     nearfield_sparse
//! \brief  Create sparse representation of nearfield part of given H-Matrix.
//! \param  A                  Matrix, of which the nearfield part is extracted.
//! \param  remove_nearfield   If true, nearfield blocks are removed from the matrix.
//!
std::unique_ptr< TSparseMatrix >
nearfield_sparse    ( TMatrix *        A,
                      const bool       remove_nearfield );

}// namespace HLIB

#endif  // __HLIB_NEARFIELD_HH
