#ifndef __HLIB_MUL_VEC_HH
#define __HLIB_MUL_VEC_HH
//
// Project     : HLib
// File        : mul_vec.hh
// Description : functions for matrix-vector multiplication
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/parallel/TProcSet.hh"
#include "hpro/cluster/TClusterBasis.hh"
#include "hpro/vector/TVector.hh"
#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

//!
//! \{
//! \name Matrix-Vector Multiplication
//! Functions for (parallel) matrix-vector multiplication
//!

//!
//! \ingroup Algebra_Module
//! \brief   compute y ≔ α·A·x + β·y on (possibly) distributed
//!          matrix A on all processors in ps
//! 
//! \param   ps     processor set defining processors involved in computation
//! \param   alpha  scaling factor for multiplication
//! \param   A      matrix to multiply with
//! \param   x      vector to multiply with
//! \param   beta   scaling factor of updated vector
//! \param   y      vector to update with multiplication result
//! \param   op     defines transformation of matrix, e.g. transposed, adjoint (\see matop_t)
//! 
void mul_vec  ( const TProcSet &  ps,
                const real alpha, const TMatrix * A, const TVector * x,
                const real  beta, TVector * y, const matop_t op );

//!
//! \ingroup Algebra_Module
//! \brief   same as mul_vec but with complex valued scalars (\see mul_vec)
//! 
void cmul_vec ( const TProcSet &  ps,
                const complex alpha, const TMatrix * A, const TVector * x,
                const complex  beta, TVector * y, const matop_t op );


//!
//! \ingroup Algebra_Module
//! \brief   compute y ≔ α·A·D·x + β·y with diagonal matrix D
//! 
//! \param   alpha  scaling factor for update
//! \param   op_A   transformation of matrix A, e.g. transposed, adjoint (\see matop_t)
//! \param   A      arbitrary matrix
//! \param   op_D   transformation of matrix D
//! \param   D      diagonal matrix
//! \param   x      source vector
//! \param   beta   scaling factor for destination
//! \param   y      destination vector
//! 
void mul_vec_diag  ( const real      alpha,
                     const matop_t   op_A, const TMatrix * A,
                     const matop_t   op_D, const TMatrix * D,
                     const TVector * x,
                     const real      beta, TVector * y );

//! \}

//
// template wrapper for TMatrix::mul_vec/cmul_vec
//
template < typename T_value >
void
mul_vec ( const T_value  alpha, const TMatrix * A, const TVector * x,
          const T_value  beta,  TVector * y,
          const matop_t  op );
          
template <>
inline
void
mul_vec< real > ( const real  alpha, const TMatrix * A, const TVector * x,
                  const real  beta,  TVector * y,
                  const matop_t  op )
{
    A->mul_vec( alpha, x, beta, y, op );
}

template <>
inline
void
mul_vec< complex > ( const complex  alpha, const TMatrix * A, const TVector * x,
                     const complex  beta,  TVector * y,
                     const matop_t  op )
{
    A->cmul_vec( alpha, x, beta, y, op );
}
          
}// namespace

#endif  // __HLIB_MUL_VEC_HH
