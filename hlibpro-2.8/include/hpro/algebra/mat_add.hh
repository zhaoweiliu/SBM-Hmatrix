#ifndef __HLIB_MAT_ADD_HH
#define __HLIB_MAT_ADD_HH
//
// Project     : HLib
// File        : mat_add.hh
// Description : matrix addition functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"

namespace HLIB
{

//!
//! \{
//! \name Matrix Addition
//! Functions for matrix addition
//!

//!
//! \ingroup Algebra_Module
//! \brief   C ≔ α·A + β·C
//!
//!          The function computes the sum \f$ C := \alpha A + \beta C \f$ with
//!          a predefined accuracy \a acc. Thread parallel execution is supported.
//!
//! \param   alpha     scaling factor for \a A
//! \param   A         update matrix for \a C
//! \param   beta      scaling factor for \a C
//! \param   C         matrix to update
//! \param   acc       accuracy of summation
//!
template < typename T_value >
void
add ( const T_value      alpha,
      const TMatrix *    A,
      const T_value      beta,
      TMatrix *          C,
      const TTruncAcc &  acc );

//!
//! \ingroup Algebra_Module
//! \brief   compute A ≔ A + λ·I
//!
template < typename T_value >
void
add_identity  ( TMatrix *      A,
                const T_value  lambda );

//!
//! \ingroup Algebra_Module
//! \brief   compute and return Σ_i A_i
//!
std::unique_ptr< TMatrix >
add ( std::list< const TMatrix * > &             matrices,
      const TTruncAcc &                          acc );

std::unique_ptr< TMatrix >
add ( std::list< std::unique_ptr< TMatrix > > &  matrices,
      const TTruncAcc &                          acc );

//! \}

}// namespace

#endif  // __HLIB_MAT_ADD_HH
