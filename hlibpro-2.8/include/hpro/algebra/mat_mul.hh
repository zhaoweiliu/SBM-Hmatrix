#ifndef __HLIB_MAT_MUL_HH
#define __HLIB_MAT_MUL_HH
//
// Project     : HLib
// File        : mat_mul.hh
// Description : matrix multiplication functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TDenseMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"

#include "hpro/parallel/dag.hh"

#include "hpro/misc/TProgressBar.hh"

#include "hpro/algebra/eval_tri.hh"

namespace HLIB
{

//!
//! \{
//! \name Matrix Multiplication
//! Functions for matrix multiplication.
//!

//!
//! \ingroup Algebra_Module
//! \brief   compute C ≔ β·C + α·op(A)·op(B)
//!
//!          The function computes to matrix product \f$C := \beta C + \alpha \tilde A \tilde B\f$
//!          where \f$\tilde A\f$ and \f$\tilde B\f$ may be the non-modified, transposed or adjoint
//!          matrices \f$A\f$ and \f$B\f$ respectively.
//!
//!          The result of the multiplication is written to \f$C\f$, whereby the block structure
//!          of \a C is not changed, e.g. the resulting block structure of the product is defined
//!          by \a C.
//!
//!          Thread-parallel execution is supported by the matrix multiplication.
//!
//!          This function is available in various versions without corresponding parameters, e.g.
//!          without \a op_A, \a op_B.
//!
//! \param   alpha     scaling factor of product
//! \param   op_A      matrix modifier for \a A
//! \param   A         first matrix factor 
//! \param   op_B      matrix modifier for \a B
//! \param   B         second matrix factor
//! \param   beta      scaling factor for \a C
//! \param   C         matrix to update
//! \param   acc       accuracy of multiplication
//! \param   progress  optional progress bar
//!
template < typename T_value >
void
multiply  ( const T_value      alpha,
            const matop_t      op_A,
            const TMatrix *    A,
            const matop_t      op_B,
            const TMatrix *    B,
            const T_value      beta ,
            TMatrix *          C,
            const TTruncAcc &  acc,
            TProgressBar *     progress = NULL );

//!
//! \ingroup Algebra_Module
//! \brief   compute C = β·C + α·A·B 
//!
//!          matrices are provided as matrix view (or in combination with matrix)
//!
template < typename T_value >
inline
void
multiply  ( const T_value      alpha,
            const TMatrix *    A,
            const TMatrix *    B,
            const T_value      beta ,
            TMatrix *          C,
            const TTruncAcc &  acc,
            TProgressBar *     progress = NULL )
{
    multiply( alpha, MATOP_NORM, A, MATOP_NORM, B, beta, C, acc, progress );
}

template < typename T_value >
inline
void
multiply  ( const T_value       alpha,
            const TMatrixView & A,
            const TMatrix *     B,
            const T_value       beta ,
            TMatrix *           C,
            const TTruncAcc &   acc,
            TProgressBar *      progress = NULL )
{
    multiply( alpha, A.op, A.M, MATOP_NORM, B, beta, C, acc, progress );
}

template < typename T_value >
inline
void
multiply  ( const T_value       alpha,
            const TMatrix *     A,
            const TMatrixView & B,
            const T_value       beta ,
            TMatrix *           C,
            const TTruncAcc &   acc,
            TProgressBar *      progress = NULL )
{
    multiply( alpha, MATOP_NORM, A, B.op, B.M, beta, C, acc, progress );
}

template < typename T_value >
inline
void
multiply  ( const T_value       alpha,
            const TMatrixView & A,
            const TMatrixView & B,
            const T_value       beta ,
            TMatrix *           C,
            const TTruncAcc &   acc,
            TProgressBar *      progress = NULL )
{
    multiply( alpha, A.op, A.M, B.op, B.M, beta, C, acc, progress );
}

//!
//! \ingroup Algebra_Module
//! \brief   compute C ≔ C + α·op(A)·op(B) using accumulators
//!
template < typename value_t >
void
multiply_accu ( const value_t      alpha,
                const matop_t      op_A,
                const TMatrix *    A,
                const matop_t      op_B,
                const TMatrix *    B,
                const value_t      beta,
                TMatrix *          C,
                const TTruncAcc &  acc );

//!
//! start accumulator based matrix multiplication C = C + α·A·B
//!
template < typename value_t >
void
add_product ( const value_t      alpha,
              const matop_t      op_A,
              const TMatrix *    A,
              const matop_t      op_B,
              const TMatrix *    B,
              TMatrix *          C,
              const TTruncAcc &  acc,
              const bool         lazy = CFG::Arith::lazy_eval );

///////////////////////////////////////////////
//
// multiplication without destination matrix
//

//!
//! \ingroup Algebra_Module
//! \brief   compute C ≔ β·C + α·op(A)·op(B)
//!
//!          The function computes to matrix product \f$\alpha \tilde A \tilde B\f$
//!          where \f$\tilde A\f$ and \f$\tilde B\f$ may be the non-modified, transposed or adjoint
//!          matrices \f$A\f$ and \f$B\f$ respectively.
//!
//!          The output format is based on the format of \f$A\f$ or \f$B\f$, which assumes that at least
//!          one of those is a non-blocked matrix. Otherwise, an exception is thrown!
//!
//! \param   alpha     scaling factor of product
//! \param   op_A      matrix modifier for \a A
//! \param   A         first matrix factor 
//! \param   op_B      matrix modifier for \a B
//! \param   B         second matrix factor
//!
template < typename T_value >
std::unique_ptr< TMatrix >
multiply   ( const T_value         alpha,
             const matop_t         op_A,
             const TMatrix *       A,
             const matop_t         op_B,
             const TMatrix *       B );


template <typename T_value>  
std::unique_ptr< TMatrix >
multiply   ( const T_value         alpha,
             const matop_t         op_A,
             const TRkMatrix *     A,
             const matop_t         op_B,
             const TMatrix *       B );
    
template <typename T_value>  
std::unique_ptr< TMatrix >
multiply   ( const T_value         alpha,
             const matop_t         op_A,
             const TMatrix *       A,
             const matop_t         op_B,
             const TRkMatrix *     B );

template <typename T_value>  
std::unique_ptr< TMatrix >
multiply   ( const T_value         alpha,
             const matop_t         op_A,
             const TDenseMatrix *  A,
             const matop_t         op_B,
             const TMatrix *       B );
    
template <typename T_value>  
std::unique_ptr< TMatrix >
multiply   ( const T_value         alpha,
             const matop_t         op_A,
             const TMatrix *       A,
             const matop_t         op_B,
             const TDenseMatrix *  B );

//!
//! \ingroup Algebra_Module
//! \brief   compute A ≔ α·L·A with lower left tridiagonal \a L
//!
//!
void multiply_ll_left  ( const real             alpha,
                         const TMatrix *        L,
                         TMatrix *              A,
                         const TTruncAcc &      acc,
                         const eval_option_t &  opts );

//!
//! \ingroup Algebra_Module
//! \brief   compute A ≔ α·A·L with lower left tridiagonal \a L
//!
//!
void multiply_ll_right ( const real             alpha,
                         TMatrix *              A,
                         const TMatrix *        L,
                         const TTruncAcc &      acc,
                         const eval_option_t &  opts );

//!
//! \ingroup Algebra_Module
//! \brief   compute C ≔ C + α·A·B with lower left tridiagonal \a B
//!
//!
void multiply_ll_right ( const real             alpha,
                         const TMatrix *        A,
                         const TMatrix *        B,
                         TMatrix *              C,
                         const TTruncAcc &      acc,
                         const eval_option_t &  opts );

//!
//! \ingroup Algebra_Module
//! \brief   compute A ≔ α·U·A with upper right tridiagonal \a U
//!
//!
void multiply_ur_left  ( const real             alpha,
                         const TMatrix *        U,
                         TMatrix *              A,
                         const TTruncAcc &      acc,
                         const eval_option_t &  opts );

//!
//! \ingroup Algebra_Module
//! \brief   compute A ≔ α·A·U with upper right tridiagonal \a U
//!
// !
void multiply_ur_right ( const real             alpha,
                         TMatrix *              A,
                         const TMatrix *        U,
                         const TTruncAcc &      acc,
                         const eval_option_t &  opts );
    
//!
//! \ingroup Algebra_Module
//! \brief   Compute A ≔ A + U·L
//!          
//!          Compute A ≔ A + U·L with lower left tridiagonal L and
//!          upper right tridiagonal U.
//!
void multiply_ur_ll    ( const TMatrix *        U,
                         const TMatrix *        L,
                         TMatrix *              A,
                         const TTruncAcc &      acc,
                         const eval_option_t &  opts_U,
                         const eval_option_t &  opts_L );
    
//!
//! \ingroup Algebra_Module
//! \brief   Compute C ≔ C + L^T·D·L.
//!          
//!          Compute C ≔ C + L^T·D·L with lower left tridiagonal L and
//!          diagonal D.
//!
void multiply_llt_d_ll ( const TMatrix *        L,
                         const TMatrix *        D,
                         TMatrix *              C,
                         const matop_t          op_L,
                         const TTruncAcc &      acc );
    
//!
//! \ingroup Algebra_Module
//! \brief   Compute B ≔ op(L)·D·B inplace, overwriting B
//!          
//!          Compute B ≔ op(L)·D·B with lower left tridiagonal L and
//!          diagonal D.
//!
void multiply_llt_d_left ( const real           alpha,
                           const TMatrix *      L,
                           const TMatrix *      D,
                           TMatrix *            B,
                           const matop_t        op_L,
                           const TTruncAcc &    acc );

//!
//! \ingroup Algebra_Module
//! \brief   Compute C ≔ C + op(L)·D·B.
//!          
//!          Compute C ≔ C + op(L)·D·B with lower left tridiagonal L,
//!          diagonal D and general B
//!
void multiply_llt_d_left ( const real           alpha,
                           const TMatrix *      L,
                           const TMatrix *      D,
                           const TMatrix *      B,
                           TMatrix *            C,
                           const matop_t        op_L,
                           const TTruncAcc &    acc );



//
// return number of steps for multiplication for progress meter initialisation
//
size_t
multiply_steps  ( const matop_t    op_A,
                  const TMatrix *  A,
                  const matop_t    op_B,
                  const TMatrix *  B,
                  const TMatrix *  C );

size_t
multiply_ll_left_steps     ( const TMatrix *  A,
                             const TMatrix *  B );

size_t
multiply_ll_right_steps    ( const TMatrix *  A,
                             const TMatrix *  B );

size_t
multiply_ur_left_steps     ( const TMatrix *  A,
                             const TMatrix *  B );

size_t
multiply_ur_right_steps    ( const TMatrix *  A,
                             const TMatrix *  B );

size_t
multiply_ur_ll_steps       ( const TMatrix *  A );

size_t
multiply_llt_d_ll_steps    ( const TMatrix *  L );

size_t
multiply_llt_d_left_steps  ( const TMatrix *  L,
                             const TMatrix *  B );

    
//! \}

}// namespace

#endif // __HLIB_MAT_MUL_HH
