#ifndef __HLIB_BLAS_ALGEBRA_HH
#define __HLIB_BLAS_ALGEBRA_HH
//
// Project     : HLib
// File        : Algebra.hh
// Description : provide linear algebra functions for BLAS matrices and vectors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#define INTERNAL_BLAS  1

#include "hpro/base/types.hh"
#include "hpro/base/traits.hh"
#include "hpro/base/System.hh"
#include "hpro/base/complex.hh"
#include "hpro/base/TTruncAcc.hh"

#include "hpro/blas/types.hh"
#include "hpro/blas/Vector.hh"
#include "hpro/blas/Matrix.hh"
#include "hpro/blas/lapack.hh"

//
// switch MKL to sequential mode before calling MKL functions
// and switch it back to standard behaviour afterwards
//
#if USE_MKL == 1 && USE_MKL_SEQ == 0

#  include <mkl_service.h>

#  define MKL_SEQ_START  auto  __old_mkl_nthreads = mkl_set_num_threads_local( 1 )
#  define MKL_SEQ_END    mkl_set_num_threads_local( __old_mkl_nthreads )

#else

#  define MKL_SEQ_START  
#  define MKL_SEQ_END    

#endif

namespace HLIB
{

namespace BLAS
{

#if HLIB_BLAS_TESTS == 1 || HLIB_DEBUG == 1
#  define HASSERT_BLAS( cond, code, fnname, msg )   { if ( ! ( cond ) ) HERROR( code, fnname, msg ); }
#else
#  define HASSERT_BLAS( cond, code, fnname, msg )
#endif

////////////////////////////////////////////////////////////////
//!
//! \{
//! \name Vector Algebra
//!
////////////////////////////////////////////////////////////////

//!
//! \ingroup  BLAS_Module
//! \brief    fill vector with constant
//!
template < typename T1,
           typename T2 >
typename enable_if< is_vector< T2 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value >::result
fill ( const T1  f,
       T2 &      x )
{
    const idx_t  n = idx_t(x.length());

    for ( idx_t  i = 0; i < n; ++i )
        x(i) = f;
}

//!
//! \ingroup  BLAS_Module
//! \brief    fill vector with repeated function evaluation, x_i = f()
//!
// template < typename T1,
//            typename T2 >
// typename enable_if< is_vector< T2 >::value >::result
// fill ( const T1  func,
//        T2 &      x )
// {
//     const idx_t  n = idx_t(x.length());

//     for ( idx_t  i = 0; i < n; ++i )
//         x(i) = func();
// }

//!
//! \ingroup  BLAS_Module
//! \brief    conjugate entries in vector
//!
template < typename T1 >
typename enable_if< is_vector< T1 >::value >::result
conj ( T1 &  x )
{
    if ( is_complex_type< typename T1::value_t >::value )
    {
        const idx_t  n = idx_t(x.length());
        
        for ( idx_t i = 0; i < n; ++i )
            x(i) = HLIB::conj( x(i) );
    }// if
}

//!
//! \ingroup  BLAS_Module
//! \brief scale vector by constant
//!
template < typename T1,
           typename T2 >
typename enable_if< is_vector< T2 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value >::result
scale ( const T1  f,
        T2 &      x )
{
    using  value_t = T1;

    MKL_SEQ_START;
    
    scal< value_t >( blas_int_t(x.length()),
                     f,
                     x.data(),
                     blas_int_t(x.stride()) );

    MKL_SEQ_END;
}

//!
//! \ingroup  BLAS_Module
//! \brief copy \a x into \a y
//!
template < typename T1,
           typename T2 >
typename enable_if< is_vector< T1 >::value  &&
                    is_vector< T2 >::value  &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
copy ( const T1 &  x,
       T2 &        y )
{
    HASSERT_BLAS( x.length() == y.length(), ERR_VEC_SIZE, "(BLAS) copy", "" );     

    using  value_t = typename T1::value_t;
    
    MKL_SEQ_START;
    
    copy< value_t >( blas_int_t(x.length()),
                     x.data(),
                     blas_int_t(x.stride()),
                     y.data(),
                     blas_int_t(y.stride()) );

    MKL_SEQ_END;
}

//!
//! \ingroup  BLAS_Module
//! \brief exchange \a x and \a y
//!
template < typename T1,
           typename T2 >
typename enable_if< is_vector< T1 >::value  &&
                    is_vector< T2 >::value  &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
swap ( T1 &  x,
       T2 &  y )
{
    HASSERT_BLAS( x.length() == y.length(), ERR_VEC_SIZE, "(BLAS) swap", "" );

    using  value_t = typename T1::value_t;

    MKL_SEQ_START;
    
    swap< value_t >( blas_int_t(x.length()),
                     x.data(),
                     blas_int_t(x.stride()),
                     y.data(),
                     blas_int_t(y.stride()) );

    MKL_SEQ_END;
}

//!
//! \ingroup  BLAS_Module
//! \brief determine index with maximal absolute value in \a x
//!
template < typename T1 >
typename enable_if_res< is_vector< T1 >::value, idx_t >::result
max_idx ( const T1 &  x )
{
    HASSERT_BLAS( x.length() > 0, ERR_VEC_SIZE, "max_idx", "" );

    using  value_t = typename T1::value_t;

    MKL_SEQ_START;
    
    auto  res = idx_t( max_idx< value_t >( blas_int_t(x.length()),
                                           x.data(),
                                           blas_int_t(x.stride()) ) - 1 );

    MKL_SEQ_END;

    return res;
}

//!
//! \ingroup  BLAS_Module
//! \brief determine index with minimax absolute value in \a x
//!
template < typename T1 >
typename enable_if_res< is_vector< T1 >::value, idx_t >::result
min_idx ( const T1 &  x )
{
    HASSERT_BLAS( x.length() > 0, ERR_VEC_SIZE, "min_idx", "" ); \

    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;
    
    const idx_t  n       = idx_t(x.length());
    real_t       min_val = Math::abs( x(0) );
    idx_t        min_idx = 0;

    for ( idx_t  i = 1; i < n; ++i )
    {
        const real_t  val = Math::abs( x(i) );
        
        if ( val < min_val )
        {
            min_val = val;
            min_idx = i;
        }// if
    }// for

    return min_idx;
}

//!
//! \ingroup  BLAS_Module
//! \brief compute y ≔ y + α·x
//!
template < typename T1,
           typename T2,
           typename T3>
typename enable_if< is_vector< T2 >::value  &&
                    is_vector< T3 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value >::result
add ( const T1    alpha,
      const T2 &  x,
      T3 &        y )
{
    HASSERT_BLAS( x.length() == y.length(), ERR_VEC_SIZE, "(BLAS) add", "" );
    
    using  value_t = T1;

    MKL_SEQ_START;
    
    axpy< value_t >( blas_int_t(x.length()),
                     alpha,
                     x.data(),
                     blas_int_t(x.stride()),
                     y.data(),
                     blas_int_t(y.stride()) );

    MKL_SEQ_END;
}

//!
//! \ingroup  BLAS_Module
//! \brief compute <x,y> = x^H · y
//!
template < typename T1,
           typename T2 >
typename enable_if_res< is_vector< T1 >::value  &&
                        is_vector< T2 >::value  &&
                        is_same_type< typename T1::value_t, typename T2::value_t >::value,
                        typename T1::value_t >::result
dot ( const T1 &  x,
      const T2 &  y )
{
    HASSERT_BLAS( x.length() == y.length(), ERR_VEC_SIZE, "(BLAS) dot", "" );

    using  value_t = typename T1::value_t;
    
    MKL_SEQ_START;
    
    auto  res = dot< value_t >( blas_int_t(x.length()),
                                x.data(),
                                blas_int_t(x.stride()),
                                y.data(),
                                blas_int_t(y.stride()) );

    MKL_SEQ_END;

    return res;
}

//!
//! \ingroup  BLAS_Module
//! \brief compute <x,y> without conjugating x, e.g. x^T · y
//!
template < typename T1,
           typename T2 >
typename enable_if_res< is_vector< T1 >::value  &&
                        is_vector< T2 >::value  &&
                        is_same_type< typename T1::value_t, typename T2::value_t >::value,
                        typename T1::value_t >::result
dotu ( const T1 &  x,
       const T2 &  y )
{
    HASSERT_BLAS( x.length() == y.length(), ERR_VEC_SIZE, "(BLAS) dot", "" );

    using  value_t = typename T1::value_t;

    MKL_SEQ_START;
    
    auto  res = dotu< value_t >( blas_int_t(x.length()),
                                 x.data(),
                                 blas_int_t(x.stride()),
                                 y.data(),
                                 blas_int_t(y.stride()) );

    MKL_SEQ_END;

    return res;
}

// helper for stable_dotu
template < typename T1 >
bool
abs_lt ( const T1  a1,
         const T1  a2 )
{
    return Math::abs( a1 ) < Math::abs( a2 );
}

//!
//! \brief  compute dot product x · y numerically stable
//! \param  x  first argument of dot product
//! \param  y  second argument of dot product
//!
template < typename T1,
           typename T2 >
typename enable_if_res< is_vector< T1 >::value  &&
                        is_vector< T2 >::value  &&
                        is_same_type< typename T1::value_t, typename T2::value_t >::value,
                        typename T1::value_t >::result
stable_dotu ( const T1 &  x,
              const T2 &  y )
{
    HASSERT_BLAS( x.length() == y.length(), ERR_VEC_SIZE, "(BLAS) dot", "" );

    using  value_t = typename T1::value_t;

    std::vector< value_t >  tmp( x.length() );
    const idx_t             n = idx_t(x.length());

    // compute conj(x_i) · y_i for each i
    for ( idx_t  i = 0; i < n; ++i )
        tmp[i] = x(i) * y(i);

    // sort tmp w.r.t. absolute value of element
    std::sort( tmp.begin(), tmp.end(), abs_lt< value_t > );

    // compute final result
    value_t  res = value_t(0);

    for ( idx_t  i = 0; i < n; ++i )
        res += tmp[i];

    return res;
}

//!
//! \brief  compute sum of elements in x numerically stable
//! \param  x  vector holding coefficients to sum up
//!
template < typename T1 >
typename enable_if_res< is_vector< T1 >::value,
                        typename T1::value_t >::result
stable_sum ( const T1 &  x )
{
    using  value_t = typename T1::value_t;
    
    std::vector< value_t >  tmp( x.length() );
    const idx_t             n = idx_t(x.length());

    // compute conj(x_i) · y_i for each i
    for ( idx_t  i = 0; i < n; ++i )
        tmp[i] = x(i);

    // sort tmp w.r.t. absolute value of element
    std::sort( tmp.begin(), tmp.end(), abs_lt< value_t > );

    // compute final result
    value_t  res = value_t(0);

    for ( idx_t  i = 0; i < n; ++i )
        res += tmp[i];

    return res;
}

//!
//! \ingroup  BLAS_Module
//! \brief compute ∥x∥₂
//!
template < typename T1 >
typename enable_if_res< is_vector< T1 >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
norm2 ( const T1 &  x )
{
    using  value_t = typename T1::value_t;
    
    MKL_SEQ_START;
    
    auto  res = norm2< value_t >( blas_int_t(x.length()),
                                  x.data(),
                                  blas_int_t(x.stride()) );

    MKL_SEQ_END;

    return res;
}

template < typename T1 >
typename enable_if_res< is_vector< T1 >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
norm_2 ( const T1 &  x )
{
    return norm2( x );
}

//! \}

////////////////////////////////////////////////////////////////
//!
//! \{
//! \name Basic Matrix Algebra
//!
////////////////////////////////////////////////////////////////

//!
//! \ingroup  BLAS_Module
//! \brief set M to f entrywise
//!
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T2 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value >::result
fill ( const T1  f,
       T2 &      M )
{
#if INTERNAL_BLAS == 1
    
    const idx_t  n = idx_t(M.nrows());
    const idx_t  m = idx_t(M.ncols());
    
    for ( idx_t j = 0; j < m; ++j )
        for ( idx_t i = 0; i < n; ++i )
            M(i,j) = f;

#else

    using  value_t = T1;
    
    MKL_SEQ_START;
    
    laset< value_t >( 'F', blas_int_t( M.nrows() ), blas_int_t( M.ncols() ),
                      f, f, M.data(), blas_int_t( M.col_stride() ) );

    MKL_SEQ_END;

#endif
}

//!
//! \ingroup  BLAS_Module
//! \brief fill M with random entries
//!
template < typename T1 >
void
fill_rand ( Matrix< T1 > &  M );

//!
//! \ingroup  BLAS_Module
//! \brief conjugate entries in vector
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
conj ( T1 &  M )
{
    using  value_t = typename T1::value_t;
    
    if ( is_complex_type< value_t >::value )
    {
        const idx_t  n = idx_t(M.nrows());
        const idx_t  m = idx_t(M.ncols());
    
        for ( idx_t j = 0; j < m; ++j )
            for ( idx_t i = 0; i < n; ++i )
                M( i, j ) = HLIB::conj( M( i, j ) );
    }// if
}

//!
//! \ingroup  BLAS_Module
//! \brief compute M ≔ f · M
//!
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T2 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value >::result
scale ( const T1  f,
        T2 &      M )
{
    const idx_t  n = idx_t(M.nrows());
    const idx_t  m = idx_t(M.ncols());
    
    for ( idx_t j = 0; j < m; ++j )
        for ( idx_t i = 0; i < n; ++i )
            M(i,j) *= f;
}

//!
//! \ingroup  BLAS_Module
//! \brief copy A to B
//!
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T1 >::value  &&
                    is_matrix< T2 >::value  &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
copy ( const T1 &  A,
       T2 &        B )
{
    HASSERT_BLAS( A.nrows() == B.nrows(), ERR_MAT_SIZE, "(BLAS) copy", "" );
    HASSERT_BLAS( A.ncols() == B.ncols(), ERR_MAT_SIZE, "(BLAS) copy", "" );

#if INTERNAL_BLAS == 1
    
    const idx_t  n = idx_t(A.nrows());
    const idx_t  m = idx_t(A.ncols());
    
    for ( idx_t j = 0; j < m; ++j )
        for ( idx_t i = 0; i < n; ++i )
            B(i,j) = A(i,j);

#else

    using  value_t = typename T1::value_t;

    MKL_SEQ_START;
    
    lacpy< value_t >( 'F', blas_int_t( A.nrows() ), blas_int_t( A.ncols() ),
                      A.data(), blas_int_t( A.col_stride() ),
                      B.data(), blas_int_t( B.col_stride() ) );

    MKL_SEQ_END;

#endif
}

//!
//! \ingroup  BLAS_Module
//! \brief transpose matrix A: A → A^T
//!        ASSUMPTION: A is square matrix
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
transpose ( T1 &  A )
{
    using  value_t = typename T1::value_t;
    
    HASSERT_BLAS( A.nrows() == A.ncols(), ERR_MAT_SIZE, "(BLAS) transpose", "" );

    const idx_t  n = idx_t(A.nrows());

    for ( idx_t i = 0; i < n-1; ++i )
    {
        Vector< value_t >  row_i( A, Range( i+1, n-1 ), i );
        Vector< value_t >  col_i( A, i, Range( i+1, n-1 ) );

        swap( row_i, col_i );
    }// for
}
            
//!
//! \ingroup  BLAS_Module
//! \brief conjugate transpose matrix A: A → A^H
//!        ASSUMPTION: A is square matrix
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
conj_transpose ( T1 &  A )
{
    using  value_t = typename T1::value_t;

    HASSERT_BLAS( A.nrows() == A.ncols(), ERR_MAT_SIZE, "(BLAS) conj_transpose", "" );

    const idx_t  n = idx_t(A.nrows());

    for ( idx_t i = 0; i < n-1; ++i )
    {
        Vector< value_t >  row_i( A, Range( i+1, n-1 ), i );
        Vector< value_t >  col_i( A, i, Range( i+1, n-1 ) );

        swap( row_i, col_i );
    }// for

    BLAS::conj( A );
}
            
//!
//! \ingroup  BLAS_Module
//! \brief determine index (i,j) with maximal absolute value in \a M
//!        and return in \a row and \a col
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
max_idx ( const T1 &  M,
          idx_t &     row,
          idx_t &     col )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;
    
    HASSERT_BLAS( M.nrows() > 0, ERR_MAT_SIZE, "(BLAS) max_idx", "" );
    HASSERT_BLAS( M.ncols() > 0, ERR_MAT_SIZE, "(BLAS) max_idx", "" );

    real_t  max_val = real_t(0);

    // init with forbidden values
    idx_t   lrow    = idx_t(M.nrows());
    idx_t   lcol    = idx_t(M.ncols());
    
    for ( idx_t  j = 0; j < idx_t(M.ncols()); ++j )
        for ( idx_t  i = 0; i < idx_t(M.nrows()); ++i )
        {
            const real_t  val = Math::abs( M(i,j) );

            if ( val > max_val )
            {
                max_val = val;
                lrow    = i;
                lcol    = j;
            }// if
        }// for

    row = lrow;
    col = lcol;
}

//!
//! \ingroup  BLAS_Module
//! \brief compute B = B + f A
//!                                                     
template < typename T1,
           typename T2,
           typename T3 >
typename enable_if< is_matrix< T2 >::value  &&
                    is_matrix< T3 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value >::result
add ( const T1    f,
      const T2 &  A,
      T3 &        B )
{
    HASSERT_BLAS( A.nrows() == B.nrows(), ERR_MAT_SIZE, "(BLAS) add", "" );
    HASSERT_BLAS( A.ncols() == B.ncols(), ERR_MAT_SIZE, "(BLAS) add", "" );

    using  value_t = T1;
    
    const idx_t  n = idx_t(A.nrows());
    const idx_t  m = idx_t(A.ncols());
    
    if (( A.col_stride() == n * A.row_stride() ) &&
        ( B.col_stride() == n * B.row_stride() ))
    {
        MKL_SEQ_START;
    
        axpy< value_t >( blas_int_t(n*m),
                         f,
                         A.data(),
                         blas_int_t(A.row_stride()),
                         B.data(),
                         blas_int_t(B.row_stride()) );

        MKL_SEQ_END;
    }// if
    else
    {
        for ( idx_t j = 0; j < m; ++j )
            for ( idx_t i = 0; i < n; ++i )
                B(i,j) += f * A(i,j);
    }// else
}

//!
//! \ingroup  BLAS_Module
//! \brief compute A ≔ A + α·x·y^H
//!
template < typename T1,
           typename T2,
           typename T3,
           typename T4 >
typename enable_if< is_vector< T2 >::value  &&
                    is_vector< T3 >::value  &&
                    is_matrix< T4 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value &&
                    is_same_type< T1, typename T4::value_t >::value >::result
add_r1 ( const T1    alpha,
         const T2 &  x,
         const T3 &  y,
         T4 &        A )
{
    HASSERT_BLAS( x.length() == A.nrows(), ERR_VEC_SIZE, "(BLAS) add_r1", "" );
    HASSERT_BLAS( y.length() == A.ncols(), ERR_VEC_SIZE, "(BLAS) add_r1", "" );

    using  value_t = T1;

    if ( A.row_stride() == 1 )
    {
        MKL_SEQ_START;
    
        ger< value_t >( blas_int_t(A.nrows()),
                        blas_int_t(A.ncols()),
                        alpha,
                        x.data(),
                        blas_int_t(x.stride()),
                        y.data(),
                        blas_int_t(y.stride()),
                        A.data(),
                        blas_int_t(A.col_stride()) );

        MKL_SEQ_END;
    }// if
    else
    {
        const idx_t  n = idx_t(A.nrows());
        const idx_t  m = idx_t(A.ncols());
    
        for ( idx_t j = 0; j < m; ++j )
        {
            const value_t  f = alpha * HLIB::conj( y(j) );
            
            for ( idx_t i = 0; i < n; ++i )
                A(i,j) += x(i) * f;
        }// for
    }// else
}

//!
//! \ingroup  BLAS_Module
//! \brief compute A ≔ A + α·x·y^T
//!
template < typename T1,
           typename T2,
           typename T3,
           typename T4 >
typename enable_if< is_vector< T2 >::value  &&
                    is_vector< T3 >::value  &&
                    is_matrix< T4 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value &&
                    is_same_type< T1, typename T4::value_t >::value >::result
add_r1u ( const T1    alpha,
          const T2 &  x,
          const T3 &  y,
          T4 &        A )
{
    HASSERT_BLAS( x.length() == A.nrows(), ERR_VEC_SIZE, "(BLAS) add_r1u", "" );
    HASSERT_BLAS( y.length() == A.ncols(), ERR_VEC_SIZE, "(BLAS) add_r1u", "" );

    using  value_t = T1;

    if ( A.row_stride() == 1 )
    {
        MKL_SEQ_START;
    
        geru< value_t >( blas_int_t(A.nrows()),
                         blas_int_t(A.ncols()),
                         alpha,
                         x.data(),
                         blas_int_t(x.stride()),
                         y.data(),
                         blas_int_t(y.stride()),
                         A.data(),
                         blas_int_t(A.col_stride()) );

        MKL_SEQ_END;
    }// if
    else
    {
        const idx_t  n = idx_t(A.nrows());
        const idx_t  m = idx_t(A.ncols());
    
        for ( idx_t j = 0; j < m; ++j )
        {
            const value_t  f = alpha * y(j);
            
            for ( idx_t i = 0; i < n; ++i )
                A(i,j) += x(i) * f;
        }// for
    }// else
}

//!
//! \ingroup  BLAS_Module
//! \brief    compute y ≔ β·y + α·A·x
//!
template < typename T1,
           typename T2,
           typename T3,
           typename T4 >
typename enable_if< is_matrix< T2 >::value  &&
                    is_vector< T3 >::value  &&
                    is_vector< T4 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value &&
                    is_same_type< T1, typename T4::value_t >::value >::result
mulvec ( const T1    alpha,
         const T2 &  A,
         const T3 &  x,
         const T1    beta,
         T4 &        y )
{
    HASSERT_BLAS( x.length() == A.ncols(), ERR_VEC_SIZE, "(BLAS) mulvec", "" );
    HASSERT_BLAS( y.length() == A.nrows(), ERR_VEC_SIZE, "(BLAS) mulvec", "" );

    using  value_t = T1;

    // if (( A.row_stride() == 1 ) && ( std::max( A.nrows(), A.ncols() ) > 8 ))
    if ( A.row_stride() == 1 )
    {
        MKL_SEQ_START;
    
        gemv< value_t >( char( A.blas_view() ),
                         blas_int_t(A.blas_nrows()),
                         blas_int_t(A.blas_ncols()),
                         alpha,
                         A.data(),
                         blas_int_t(A.col_stride()),
                         x.data(),
                         blas_int_t(x.stride()),
                         beta,
                         y.data(),
                         blas_int_t(y.stride()) );

        MKL_SEQ_END;
    }// if
    else
    {
        const idx_t  n = idx_t(A.nrows());
        const idx_t  m = idx_t(A.ncols());
        
        if ( beta == value_t(0) )
        {
            for ( idx_t i = 0; i < n; ++i )
                y(i) = value_t(0);
        }// if
        else if ( beta != value_t(1) )
        {
            for ( idx_t i = 0; i < n; ++i )
                y(i) *= beta;
        }// if
        
        for ( idx_t i = 0; i < n; ++i )
        {
            value_t  f = value_t(0);
            
            for ( idx_t j = 0; j < m; ++j )
                f += A(i,j) * x(j);
            
            y(i) += alpha * f;
        }// for
    }// else
}

//!
//! \ingroup  BLAS_Module
//! \brief    compute y ≔ α·A·x
//!
template < typename T1,
           typename T2,
           typename T3 >
typename enable_if_res< is_matrix< T2 >::value  &&
                        is_vector< T3 >::value  &&
                        is_same_type< T1, typename T2::value_t >::value &&
                        is_same_type< T1, typename T3::value_t >::value,
                        Vector< typename T2::value_t > >::result
mulvec ( const T1    alpha,
         const T2 &  A,
         const T3 &  x )
{
    HASSERT_BLAS( x.length() == A.ncols(), ERR_VEC_SIZE, "(BLAS) mulvec", "" );

    using  value_t = T1;

    Vector< value_t >  y( A.nrows() );
    
    // if (( A.row_stride() == 1 ) && ( std::max( A.nrows(), A.ncols() ) > 8 ))
    if ( A.row_stride() == 1 )
    {
        MKL_SEQ_START;
    
        gemv< value_t >( char( A.blas_view() ),
                         blas_int_t(A.blas_nrows()),
                         blas_int_t(A.blas_ncols()),
                         alpha,
                         A.data(),
                         blas_int_t(A.col_stride()),
                         x.data(),
                         blas_int_t(x.stride()),
                         value_t(1),
                         y.data(),
                         blas_int_t(y.stride()) );

        MKL_SEQ_END;
    }// if
    else
    {
        const idx_t  n = idx_t(A.nrows());
        const idx_t  m = idx_t(A.ncols());
        
        for ( idx_t i = 0; i < n; ++i )
        {
            value_t  f = value_t(0);
            
            for ( idx_t j = 0; j < m; ++j )
                f += A(i,j) * x(j);
            
            y(i) += alpha * f;
        }// for
    }// else

    return y;
}

//!
//! \ingroup  BLAS_Module
//! \brief    compute x ≔ M·x, where M is upper or lower triangular
//!           with unit or non-unit diagonal
//!
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T1 >::value  &&
                    is_vector< T2 >::value  &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
mulvec_tri ( const tri_type_t   shape,
             const diag_type_t  diag,
             const T1 &         A,
             T2 &               x )
{
    HASSERT_BLAS( x.length() == A.nrows(), ERR_VEC_SIZE, "(BLAS) mulvec_tri", "" );
    HASSERT_BLAS( x.length() == A.ncols(), ERR_VEC_SIZE, "(BLAS) mulvec_tri", "" );

    using  value_t = typename T1::value_t;

    if ( A.row_stride() == 1 )
    {
        MKL_SEQ_START;
    
        trmv< value_t >( char( shape ),
                         char( A.blas_view() ),
                         char( diag ),
                         blas_int_t(A.blas_nrows()), 
                         A.data(),
                         blas_int_t(A.col_stride()),
                         x.data(),
                         blas_int_t(x.stride()) );

        MKL_SEQ_END;
    }// if
    else
    {
        HERROR( ERR_NOT_IMPL, "", "" );
    }// else
}

//!
//! \ingroup  BLAS_Module
//! \brief solve A·x = b with known \a A and \a b; x overwrites \a b
//!        (\a A is overwritten upon exit!)
//! 
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T1 >::value  &&
                    is_vector< T2 >::value  &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
solve  ( T1 &  A,
         T2 &  b )
{
    using  value_t = typename T1::value_t;
    
    if ( A.nrows() != A.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) solve", "matrix is not square" );
    
    if ( A.ncols() != b.length() )
        HERROR( ERR_VEC_SIZE, "(BLAS) solve", "vector has wrong dimension" );
    
    std::vector< blas_int_t >  ipiv( A.nrows() );
    blas_int_t                 info = 0;

    MKL_SEQ_START;
    
    gesv< value_t >( blas_int_t( A.nrows() ), blas_int_t( 1 ),
                     A.data(), blas_int_t( A.col_stride() ), ipiv.data(),
                     b.data(), blas_int_t( 1 ), info );

    MKL_SEQ_END;
    
    if      ( info < 0 ) HERROR( ERR_ARG, "(BLAS) solve",
                                 to_string( "argument %d to LAPACK::gesv", -info ) );
    else if ( info > 0 ) HERROR( ERR_MAT_SINGULAR, "(BLAS) solve",
                                 to_string( "in LAPACK::gesv (info = %d)", info ) );
}

//!
//! \ingroup  BLAS_Module
//! \brief compute C ≔ β·C + α·A·B
//!
template < typename T1,
           typename T2,
           typename T3,
           typename T4 >
typename enable_if< is_matrix< T2 >::value  &&
                    is_matrix< T3 >::value  &&
                    is_matrix< T4 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value &&
                    is_same_type< T1, typename T4::value_t >::value >::result
prod ( const T1    alpha,
       const T2 &  A,
       const T3 &  B,
       const T1    beta,
       T4 &        C )
{
    HASSERT_BLAS( A.nrows() == C.nrows(), ERR_MAT_SIZE, "(BLAS) prod", "" );
    HASSERT_BLAS( B.ncols() == C.ncols(), ERR_MAT_SIZE, "(BLAS) prod", "" );
    HASSERT_BLAS( A.ncols() == B.nrows(), ERR_MAT_SIZE, "(BLAS) prod", "" );
    
    using  value_t = T1;

    MKL_SEQ_START;
    
    gemm< value_t >( char( A.blas_view() ),
                     char( B.blas_view() ),
                     blas_int_t(C.nrows()),
                     blas_int_t(C.ncols()),
                     blas_int_t(A.ncols()),
                     alpha, A.data(),
                     blas_int_t(A.col_stride()),
                     B.data(),
                     blas_int_t(B.col_stride()),
                     beta,
                     C.data(),
                     blas_int_t(C.col_stride()) );

    MKL_SEQ_END;
}

template < typename T1,
           typename T2,
           typename T3,
           typename T4 >
typename enable_if< is_matrix< T2 >::value  &&
                    is_matrix< T3 >::value  &&
                    is_matrix< T4 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value &&
                    is_same_type< T1, typename T4::value_t >::value >::result
prod_par ( const T1    alpha,
           const T2 &  A,
           const T3 &  B,
           const T1    beta,
           T4 &        C )
{
    HASSERT_BLAS( A.nrows() == C.nrows(), ERR_MAT_SIZE, "(BLAS) prod", "" );
    HASSERT_BLAS( B.ncols() == C.ncols(), ERR_MAT_SIZE, "(BLAS) prod", "" );
    HASSERT_BLAS( A.ncols() == B.nrows(), ERR_MAT_SIZE, "(BLAS) prod", "" );
    
    using  value_t = T1;

    gemm< value_t >( char( A.blas_view() ),
                     char( B.blas_view() ),
                     blas_int_t(C.nrows()),
                     blas_int_t(C.ncols()),
                     blas_int_t(A.ncols()),
                     alpha, A.data(),
                     blas_int_t(A.col_stride()),
                     B.data(),
                     blas_int_t(B.col_stride()),
                     beta,
                     C.data(),
                     blas_int_t(C.col_stride()) );
}

//!
//! \ingroup  BLAS_Module
//! \brief compute C ≔ α·A·B
//!
template < typename T1,
           typename T2,
           typename T3 >
typename enable_if_res< is_matrix< T2 >::value  &&
                        is_matrix< T3 >::value  &&
                        is_same_type< T1, typename T2::value_t >::value &&
                        is_same_type< T1, typename T3::value_t >::value,
                        Matrix< typename T2::value_t > >::result
prod ( const T1    alpha,
       const T2 &  A,
       const T3 &  B )
{
    HASSERT_BLAS( A.ncols() == B.nrows(), ERR_MAT_SIZE, "(BLAS) prod", "" );
    
    using  value_t = T1;

    Matrix< value_t >  C( A.nrows(), B.ncols() );

    MKL_SEQ_START;
    
    gemm< value_t >( char( A.blas_view() ),
                     char( B.blas_view() ),
                     blas_int_t(C.nrows()),
                     blas_int_t(C.ncols()),
                     blas_int_t(A.ncols()),
                     alpha, A.data(),
                     blas_int_t(A.col_stride()),
                     B.data(),
                     blas_int_t(B.col_stride()),
                     value_t(1),  // C is initialised with zero
                     C.data(),
                     blas_int_t(C.col_stride()) );

    MKL_SEQ_END;

    return C;
}

template < typename T1,
           typename T2,
           typename T3 >
typename enable_if_res< is_matrix< T2 >::value  &&
                        is_matrix< T3 >::value  &&
                        is_same_type< T1, typename T2::value_t >::value &&
                        is_same_type< T1, typename T3::value_t >::value,
                        Matrix< typename T2::value_t > >::result
prod_par ( const T1    alpha,
           const T2 &  A,
           const T3 &  B )
{
    HASSERT_BLAS( A.ncols() == B.nrows(), ERR_MAT_SIZE, "(BLAS) prod", "" );
    
    using  value_t = T1;

    Matrix< value_t >  C( A.nrows(), B.ncols() );

    gemm< value_t >( char( A.blas_view() ),
                     char( B.blas_view() ),
                     blas_int_t(C.nrows()),
                     blas_int_t(C.ncols()),
                     blas_int_t(A.ncols()),
                     alpha, A.data(),
                     blas_int_t(A.col_stride()),
                     B.data(),
                     blas_int_t(B.col_stride()),
                     value_t(1),  // C is initialised with zero
                     C.data(),
                     blas_int_t(C.col_stride()) );

    return C;
}

//!
//! \ingroup  BLAS_Module
//! \brief compute B ≔ A⊙B, e.g. Hadamard product
//!
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T1 >::value  &&
                    is_matrix< T2 >::value  &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
hadamard_prod ( const T1 &  A,
                T2 &        B )
{
    HASSERT_BLAS( A.nrows() == B.nrows(), ERR_MAT_SIZE, "(BLAS) hadamard_prod", "" );
    HASSERT_BLAS( A.ncols() == B.ncols(), ERR_MAT_SIZE, "(BLAS) hadamard_prod", "" );

    const idx_t  nrows = idx_t(A.nrows());
    const idx_t  ncols = idx_t(A.ncols());
    
    for ( idx_t  j = 0; j < ncols; ++j )
        for ( idx_t  i = 0; i < nrows; ++i )
            B(i,j) *= A(i,j);
}

//!
//! \ingroup  BLAS_Module
//! \brief    compute B ≔ α·A·B or B ≔ α·B·A with triangular matrix A
//! 
template < typename T1,
           typename T2,
           typename T3 >
typename enable_if< is_matrix< T2 >::value  &&
                    is_matrix< T3 >::value  &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value >::result
prod_tri  ( const eval_side_t  side,
            const tri_type_t   uplo,
            const diag_type_t  diag,
            const T1           alpha,
            const T2 &         A,
            T3 &               B )
{
    HASSERT_BLAS( A.nrows() == A.ncols(), ERR_MAT_SIZE, "(BLAS) prod_tri", "" );
    HASSERT_BLAS( ( side == from_left  ) && ( A.ncols() == B.nrows() ),
             ERR_MAT_SIZE, "(BLAS) prod_tri", "" );
    HASSERT_BLAS( ( side == from_right ) && ( A.nrows() == B.ncols() ),
             ERR_MAT_SIZE, "(BLAS) prod_tri", "" );

    using  value_t = T1;
    
    MKL_SEQ_START;
    
    trmm< value_t >( char(side),
                     char(uplo),
                     char( A.blas_view() ),
                     char(diag),
                     blas_int_t( A.nrows() ),
                     blas_int_t( A.ncols() ),
                     alpha, A.data(),
                     blas_int_t( A.col_stride() ),
                     B.data(),
                     blas_int_t( B.col_stride() ) );

    MKL_SEQ_END;
}

//!
//! \ingroup  BLAS_Module
//! \brief    multiply \a k columns of \a M with diagonal matrix \a D,
//!           e.g. compute M ≔ M·D
//!
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T1 >::value  &&
                    is_vector< T2 >::value  &&
                    is_same_type< typename real_type< typename T1::value_t >::type_t,
                                  typename real_type< typename T2::value_t >::type_t >::value >::result
prod_diag ( T1 &         M,
            const T2 &   D,
            const idx_t  k )
{
    using  value_t = typename T1::value_t;
    
    const Range  row_is( 0, idx_t( M.nrows() )-1 );
    
    for ( idx_t  i = 0; i < k; ++i )
    {
        Vector< value_t >  Mi( M, row_is, i );

        scale( value_t(D(i)), Mi );
    }// for
}

//!
//! \ingroup  BLAS_Module
//! \brief    multiply \a k rows of \a M with diagonal matrix \a D,
//!           e.g. compute M ≔ D·M
//!
template < typename T1,
           typename T2 >
typename enable_if< is_vector< T1 >::value  &&
                    is_matrix< T2 >::value  &&
                    is_same_type< typename real_type< typename T1::value_t >::type_t,
                                  typename real_type< typename T2::value_t >::type_t >::value >::result
prod_diag ( const T1 &   D,
            T2 &         M,
            const idx_t  k )
{
    using  value_t = typename T1::value_t;
    
    const Range  col_is( 0, idx_t( M.ncols() )-1 );
    
    for ( idx_t  i = 0; i < k; ++i )
    {
        Vector< value_t >  Mi( M, i, col_is );

        scale( value_t(D(i)), Mi );
    }// for
}

//!
//! \ingroup  BLAS_Module
//! \brief return spectral norm of M
//!
template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
norm2 ( const T1 & M );

template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
norm_2 ( const T1 & M )
{
    return norm2( M );
}

//!
//! \ingroup  BLAS_Module
//! \brief return Frobenius norm of M
//!
template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
normF ( const T1 & M )
{
    using  value_t = typename T1::value_t;
    
    const idx_t  n = idx_t(M.nrows());
    const idx_t  m = idx_t(M.ncols());

    if ( M.col_stride() == n * M.row_stride() )
    {
        return norm2( blas_int_t(n*m),
                      M.data(),
                      blas_int_t(M.row_stride()) );
    }// if
    else
    {
        value_t  f = value_t(0);
    
        for ( idx_t j = 0; j < m; ++j )
            for ( idx_t i = 0; i < n; ++i )
                f += HLIB::conj( M(i,j) ) * M(i,j);
        
        return HLIB::re( Math::sqrt( f ) );
    }// else
}

template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
norm_F ( const T1 & M )
{
    return normF( M );
}

//!
//! \brief compute Frobenius norm of A-B
//! \ingroup  BLAS_Module
//!
template < typename T1,
           typename T2 >
typename enable_if_res< is_matrix< T1 >::value  &&
                        is_matrix< T2 >::value  &&
                        is_same_type< typename T1::value_t, typename T2::value_t >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
diff_normF ( const T1 &  A,
             const T2 &  B )
{
    HASSERT_BLAS( A.nrows() == B.nrows(), ERR_MAT_SIZE, "(BLAS) diff_normF", "matrix size of A and B differs" );
    HASSERT_BLAS( A.ncols() == B.ncols(), ERR_MAT_SIZE, "(BLAS) diff_normF", "matrix size of A and B differs" );
    
    using  value_t = typename T1::value_t;

    const idx_t  n = idx_t(A.nrows());
    const idx_t  m = idx_t(A.ncols());
    value_t      f = value_t(0);
    
    for ( idx_t j = 0; j < m; ++j )
    {
        for ( idx_t i = 0; i < n; ++i )
        {
            const value_t  a_ij = A(i,j) - B(i,j);
            
            f += HLIB::conj( a_ij ) * a_ij;
        }// for
    }// for

    return HLIB::re( Math::sqrt( f ) );
}

template < typename T1,
           typename T2 >
typename enable_if_res< is_matrix< T1 >::value  &&
                        is_matrix< T2 >::value  &&
                        is_same_type< typename T1::value_t, typename T2::value_t >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
diff_norm_F ( const T1 &  A,
              const T2 &  B )
{
    return diff_normF( A, B );
}

//!
//! \brief compute Frobenius norm of M=A·B^H
//! \ingroup  BLAS_Module
//!
template < typename T1,
           typename T2 >
typename enable_if_res< is_matrix< T1 >::value  &&
                        is_matrix< T2 >::value  &&
                        is_same_type< typename T1::value_t, typename T2::value_t >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
lr_normF ( const T1 &  A,
           const T2 &  B )
{
    HASSERT_BLAS( A.ncols() != B.ncols(), ERR_MAT_SIZE, "(BLAS) lr_normF", "rank of A and B differs" );
    
    using  value_t = typename T1::value_t;

    //
    // ∑_ij (M_ij)² = ∑_ij (∑_k A_ik B_jk')²
    //              = ∑_ij (∑_k A_ik B_jk') (∑_l A_il B_jl')'
    //              = ∑_ij ∑_k ∑_l A_ik B_jk' A_il' B_jl
    //              = ∑_k ∑_l ∑_i A_ik A_il' ∑_j B_jk' B_jl
    //              = ∑_k ∑_l (A_l)^H · A_k  B_k^H · B_l
    //              = ∑_kl ( A^H·A ⦻ B^H·B )    { ⦻ : pointwise product }
    //

    const idx_t        k = idx_t(A.ncols());
    Matrix< value_t >  AtA( k, k );
    Matrix< value_t >  BtB( k, k );

    prod( value_t(1), adjoint( A ), A, value_t(0), AtA );
    prod( value_t(1), adjoint( B ), B, value_t(0), BtB );

    value_t  t = 0;
    
    for ( idx_t  j = 0; j < k; ++j )
        for ( idx_t  i = 0; i < k; ++i )
            t += AtA( i, j ) * BtB( i, j );

    return re( Math::sqrt( t ) );
}

template < typename T1,
           typename T2 >
typename enable_if_res< is_matrix< T1 >::value  &&
                        is_matrix< T2 >::value  &&
                        is_same_type< typename T1::value_t, typename T2::value_t >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
lr_norm_F ( const T1 &  A,
            const T2 &  B )
{
    return lr_normF( A, B );
}

//!
//! \ingroup  BLAS_Module
//! \brief return condition of M
//!        - if M ≡ 0 or size(M) = 0, 0 is returned
//!
template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value,
                        typename real_type< typename T1::value_t >::type_t >::result
cond ( const T1 &  M );

//!
//! \ingroup  BLAS_Module
//! \brief    make given matrix symmetric, e.g. copy lower left part to
//!           upper right part
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
make_symmetric ( T1 &  A )
{
    const idx_t  m = idx_t( A.ncols() );
    
    for ( idx_t  j = 1; j < m; j++ )
        for ( idx_t  i = 0; i < j; i++ )
            A( i, j ) = A( j, i );
}

//!
//! \ingroup  BLAS_Module
//! \brief    make given matrix hermitian, e.g. copy conjugated lower left part
//!           to upper right part and make diagonal real
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
make_hermitian ( T1 &  A )
{
    using  value_t = typename T1::value_t;
    using  real_t  = typename real_type< value_t >::type_t;

    const idx_t  n = idx_t( A.nrows() );
    const idx_t  m = idx_t( A.ncols() );
    
    for ( idx_t  j = 1; j < m; j++ )
        for ( idx_t  i = 0; i < j; i++ )
            A( i, j ) = HLIB::conj( A( j, i ) );

    // reset imaginary part of diagonal (assuming that it is small !!!)
    for ( idx_t  i = 0; i < std::min( n, m ); i++ )
    {
        const value_t  z = A( i, i );

        if ( im( z ) != real_t(0) )
            A( i, i ) = real_t(0.5) * ( z + HLIB::conj(z) );
    }// for
}

//! \}

////////////////////////////////////////////////////////////////
//!
//! \{
//! \name Advanced Matrix Algebra
//!
////////////////////////////////////////////////////////////////

//!
//! \ingroup  BLAS_Module
//! \brief solve A·X = B with known \a A and \a B; X overwrites \a B
//! 
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T1 >::value &&
                    is_matrix< T2 >::value &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
solve  ( T1 &  A,
         T2 &  B )
{
    using  value_t = typename T1::value_t;

    if ( A.nrows() != A.ncols() )
        HERROR( ERR_MAT_SIZE, "(BLAS) solve", "matrix is not square" );
    
    if ( A.ncols() != B.nrows() )
        HERROR( ERR_MAT_SIZE, "(BLAS) solve", "#column(A) != #rows(B)" );
    
    std::vector< blas_int_t >  ipiv( A.nrows() );
    blas_int_t                 info = 0;

    MKL_SEQ_START;
    
    gesv< value_t >( blas_int_t( A.nrows() ), blas_int_t( B.ncols() ),
                     A.data(), blas_int_t( A.col_stride() ), ipiv.data(),
                     B.data(), blas_int_t( B.col_stride() ), info );

    MKL_SEQ_END;
    
    if      ( info < 0 ) HERROR( ERR_ARG, "(BLAS) solve",
                                 to_string( "argument %d to LAPACK::gesv", -info ) );
    else if ( info > 0 ) HERROR( ERR_MAT_SINGULAR, "(BLAS) solve",
                                 to_string( "in LAPACK::gesv (info = %d)", info ) );
}

//!
//! \ingroup  BLAS_Module
//! \brief    solve A·x = b with known \a A and \a b; x overwrites \a b
//! 
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T1 >::value &&
                    is_vector< T2 >::value &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
solve_tri  ( const tri_type_t   uplo,
             const diag_type_t  diag,
             const T1 &         A,
             T2 &               b )
{
    HASSERT_BLAS( A.nrows() == A.ncols(),  ERR_MAT_SIZE, "(BLAS) solve_tri", "matrix is not square" );
    HASSERT_BLAS( A.nrows() == b.length(), ERR_MAT_SIZE, "(BLAS) solve_tri", "matrix and vector size differ" );

    using  value_t = typename T1::value_t;

    MKL_SEQ_START;
    
    trsv< value_t >( char(uplo),
                     char( A.blas_view() ),
                     char(diag),
                     blas_int_t( A.nrows() ),
                     A.data(),
                     blas_int_t( A.col_stride() ),
                     b.data(),
                     blas_int_t( b.stride() ) );

    MKL_SEQ_END;
}

//!
//! \ingroup  BLAS_Module
//! \brief    solve A·X = α·B or X·A· = α·B with known \a A and
//!           \a B; X overwrites \a B
//! 
template < typename T1,
           typename T2,
           typename T3 >
typename enable_if< is_matrix< T2 >::value &&
                    is_matrix< T3 >::value &&
                    is_same_type< T1, typename T2::value_t >::value &&
                    is_same_type< T1, typename T3::value_t >::value >::result
solve_tri  ( const eval_side_t  side,
             const tri_type_t   uplo,
             const diag_type_t  diag,
             const T1           alpha,
             const T2 &         A,
             T3 &               B )
{
    HASSERT_BLAS( A.nrows() == A.ncols(), ERR_MAT_SIZE, "(BLAS) solve_tri", "matrix is not square" );

    using  value_t = T1;

    MKL_SEQ_START;
    
    trsm< value_t >( char(side),
                     char(uplo),
                     char( A.blas_view() ),
                     char(diag),
                     blas_int_t( B.nrows() ),
                     blas_int_t( B.ncols() ),
                     alpha, A.data(),
                     blas_int_t( A.col_stride() ),
                     B.data(),
                     blas_int_t( B.col_stride() ) );

    MKL_SEQ_END;
}

//!
//! \ingroup  BLAS_Module
//! \brief invert matrix \a A; \a A will be overwritten with A^-1 upon exit
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
invert ( T1 &  A );

//!
//! \ingroup  BLAS_Module
//! \brief invert lower or upper triangular matrix \a A with unit or non-unit diagonal;
//!        \a A will be overwritten with A^-1 upon exit
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
invert ( T1 &               A,
         const tri_type_t   tri_type,
         const diag_type_t  diag_type );

//!
//! \ingroup  BLAS_Module
//! \brief    compute pseudo inverse of matrix \a A with precision \a acc
//!
//! \detail   Compute pseudo inverse B of matrix \a A up to precision \a acc,
//!           e.g. \f$\|A-B\|\le \epsilon\f$ with \f$\epsilon\f$ defined by \a acc.
//!
template <typename T>
void
pseudo_invert ( Matrix< T > &      A,
                const TTruncAcc &  acc );

//!
//! \ingroup  BLAS_Module
//! \brief compute LU factorisation of the n×m matrix \a A with
//!        n×min(n,m) unit diagonal lower triangular matrix L and
//!        min(n,m)xm upper triangular matrix U; \a A will be
//!        overwritten with L <b>and</b> U upon exit
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
lu     ( T1 &  A );

//!
//! \ingroup  BLAS_Module
//! \brief compute L·L^T factorisation of given symmetric n×n matrix \a A
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
llt    ( T1 &  A );

//!
//! \ingroup  BLAS_Module
//! \brief compute L·L^H factorisation of given hermitian n×n matrix \a A
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
llh    ( T1 &  A );

//!
//! \ingroup  BLAS_Module
//! \brief compute L·D·L^T factorisation of given symmetric n×n matrix \a A
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
ldlt   ( T1 &  A );

//!
//! \ingroup  BLAS_Module
//! \brief compute L·D·L^H factorisation of given hermitian n×n matrix \a A
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
ldlh   ( T1 &  A );

//!
//! \ingroup  BLAS_Module
//! \brief Compute QR factorisation of the matrix \a A.
//!
//!        Compute QR factorisation of the n×m matrix \a A with
//!        n×m matrix Q and mxm matrix R (n >= m); \a A will be
//!        overwritten with Q upon exit
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
qr     ( T1 &                              A,
         Matrix< typename T1::value_t > &  R );

//!
//! \ingroup  BLAS_Module
//! \brief Compute QR factorisation of the tall-and-skinny matrix \a A.
//!
//!        Compute QR factorisation of the n×m matrix \a A, m ≪ n, with
//!        n×m matrix Q and mxm matrix R (n >= m); \a A will be
//!        overwritten with Q upon exit
//!        - ntile defines tile size (0: use internal default tile size)
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
tsqr   ( T1 &                              A,
         Matrix< typename T1::value_t > &  R,
         const size_t                      ntile = 0 );

//!
//! \ingroup  BLAS_Module
//! \brief Compute QR factorisation with column pivoting of the matrix \a A.
//!
//!        Compute QR factorisation with column pivoting of the n×m matrix \a A,
//!        i.e., \f$ AP = QR \f$, with n×m matrix Q and mxm matrix R (n >= m).
//!        \a P_i = j means, column i of AP was column j of A; \a A will be
//!        overwritten with Q upon exit
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
qrp    ( T1 &                              A,
         Matrix< typename T1::value_t > &  R,
         std::vector< blas_int_t > &       P );

//!
//! \ingroup  BLAS_Module
//! \brief compute eigenvalues and eigenvectors of matrix \a M
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
eigen ( T1 &                              M,
        Vector< typename T1::value_t > &  eig_val,
        Matrix< typename T1::value_t > &  eig_vec );

//!
//! \ingroup  BLAS_Module
//! \brief compute selected (by \a eig_range) eigenvalues and eigenvectors 
//!        of the \b symmetric matrix \a M
//!        - the lower half of \a M is accessed
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
eigen ( T1 &                              M,
        const Range &                     eig_range,
        Vector< typename T1::value_t > &  eig_val,
        Matrix< typename T1::value_t > &  eig_vec );

//!
//! \ingroup  BLAS_Module
//! \brief compute eigenvalues and eigenvectors 
//!        of the \b symmetric, \b tridiagonal matrix defines by diagonal
//!        coefficients in \a diag and off-diagonal coefficients \a subdiag
//!
template < typename T1,
           typename T2 >
typename enable_if< is_vector< T1 >::value &&
                    is_vector< T2 >::value &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
eigen ( T1 &  diag,
        T2 &  subdiag,
        Vector< typename T1::value_t > &  eig_val,
        Matrix< typename T1::value_t > &  eig_vec );

//!
//! \ingroup  BLAS_Module
//! \brief compute SVD decomposition \f$ A = U·S·V^H \f$ of the nxm matrix \a A with
//!        n×min(n,m) matrix U, min(n,m)×min(n,m) matrix S (diagonal)
//!        and m×min(n,m) matrix V; \a A will be overwritten with U upon exit
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
svd    ( T1 &                                                            A,
         Vector< typename real_type< typename T1::value_t >::type_t > &  S,
         Matrix< typename T1::value_t > &                                V );

//!
//! \ingroup  BLAS_Module
//! \brief compute SVD decomposition \f$ A = U·S·V^H \f$ of the nxm matrix \a A
//!        but return only the left/right singular vectors and the
//!        singular values S ∈ ℝ^min(n,m);
//!        upon exit, \a A will be contain the corresponding sing. vectors
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
svd    ( T1 &                                                            A,
         Vector< typename real_type< typename T1::value_t >::type_t > &  S,
         const bool                                                      left = true );

//!
//! \ingroup  BLAS_Module
//! \brief compute SVD decomposition \f$ A = U·S·V^H \f$ of the nxm 
//!        matrix \a A but return only the singular values S ∈ ℝ^min(n,m);
//!        \a A will be overwritten with U upon exit
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
sv     ( T1 &                                                            A,
         Vector< typename real_type< typename T1::value_t >::type_t > &  S );

//!
//! \ingroup  BLAS_Module
//! \brief compute SVD decomposition \f$ M = A·B^H = U·S·V^H \f$ of the nxm 
//!        low-rank matrix \a M but return only the singular values S ∈ ℝ^min(n,m);
//!        \a A and \a B will be overwritten upon exit
//!
template < typename T1,
           typename T2 >
typename enable_if< is_matrix< T1 >::value &&
                    is_matrix< T2 >::value &&
                    is_same_type< typename T1::value_t, typename T2::value_t >::value >::result
sv     ( T1 &                                                            A,
         T2 &                                                            B,
         Vector< typename real_type< typename T1::value_t >::type_t > &  S );

//!
//! \ingroup  BLAS_Module
//! \brief approximate given dense matrix \a M by low rank matrix
//!        according to accuracy \a acc. The low rank matrix will
//!        be stored in \a A and \a B
//!
template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value, size_t >::result
approx ( T1 &                              M,
         const TTruncAcc &                 acc,
         Matrix< typename T1::value_t > &  A,
         Matrix< typename T1::value_t > &  B );

//!
//! \ingroup  BLAS_Module
//! \brief approximate given dense matrix \a M by low rank matrix
//!        according to accuracy \a acc using SVD. The low rank matrix
//!        will be stored in \a A and \a B
//!
template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value, size_t >::result
approx_svd ( T1 &                              M,
             const TTruncAcc &                 acc,
             Matrix< typename T1::value_t > &  A,
             Matrix< typename T1::value_t > &  B );

template < typename T >
std::pair< Matrix< T >, Matrix< T > >
approx_svd ( Matrix< T > &                     M,
             const TTruncAcc &                 acc );

//!
//! \ingroup  BLAS_Module
//! \brief approximate given dense matrix \a M by low rank matrix
//!        according to accuracy \a acc using RRQR. The 
//!        low rank matrix will be stored in \a A and \a B. 
//!
template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value, size_t >::result
approx_rrqr ( T1 &                              M,
              const TTruncAcc &                 acc,
              Matrix< typename T1::value_t > &  A,
              Matrix< typename T1::value_t > &  B );

//!
//! \ingroup  BLAS_Module
//! \brief approximate given dense matrix \a M by low rank matrix
//!        according to accuracy \a acc using randomized SVD. The 
//!        low rank matrix will be stored in \a A and \a B. 
//!
template < typename T1 >
typename enable_if_res< is_matrix< T1 >::value, size_t >::result
approx_randsvd ( T1 &                              M,
                 const TTruncAcc &                 acc,
                 Matrix< typename T1::value_t > &  A,
                 Matrix< typename T1::value_t > &  B,
                 const uint                        power_steps  = CFG::BLAS::power_steps,
                 const uint                        oversampling = CFG::BLAS::oversampling );

//!
//! \ingroup  BLAS_Module
//! \brief truncate \a A · \a B^H based on accuracy \a acc.
//!
//!        Truncate given \a A · \a B^H low rank matrix matrix 
//!        (\a A being n×k and \a B being m×k) with respect to
//!        given accuracy \a acc; store truncated matrix in
//!        A(:,0:k-1) and B(:,0:k-1) where k is the returned
//!        new rank after truncation
//!
template <typename T>
size_t
truncate ( Matrix< T > &      A,
           Matrix< T > &      B,
           const TTruncAcc &  acc );

//!
//! \ingroup  BLAS_Module
//! \brief truncate \a A · \a B^H based on accuracy \a acc using SVD.
//!
template <typename T>
size_t
truncate_svd  ( Matrix< T > &      A,
                Matrix< T > &      B,
                const TTruncAcc &  acc );

template <typename T>
std::pair< Matrix< T >, Matrix< T > >
truncate2_svd  ( const Matrix< T > &  A,
                 const Matrix< T > &  B,
                 const TTruncAcc &    acc );

//!
//! \ingroup  BLAS_Module
//! \brief truncate \a A · \a B^H based on accuracy \a acc using RRQR.
//!
template <typename T>
size_t
truncate_rrqr ( Matrix< T > &      A,
                Matrix< T > &      B,
                const TTruncAcc &  acc );

//!
//! \ingroup  BLAS_Module
//! \brief truncate \a A · \a B^H based on accuracy \a acc using randomized SVD.
//!
template <typename T>
size_t
truncate_rand ( Matrix< T > &      A,
                Matrix< T > &      B,
                const TTruncAcc &  acc );

//!
//! \ingroup  BLAS_Module
//! \brief construct factorisation A = Q·R of \a A, with orthonormal Q
//! - A is overwritten with Q upon exit
//!
template < typename T1 >
typename enable_if< is_matrix< T1 >::value >::result
factorise_ortho ( T1 &                              A,
                  Matrix< typename T1::value_t > &  R );

//!
//! \ingroup  BLAS_Module
//! \brief construct approximate factorisation A = Q·R of \a A, with orthonormal Q
//!        - approximation quality is definded by \a acc
//!        - A is overwritten with Q upon exit
//!
template < typename T >
void
factorise_ortho ( Matrix< T > &      A,
                  Matrix< T > &      R,
                  const TTruncAcc &  acc );

//!
//! \ingroup  BLAS_Module
//! \brief print statistics for Algebra functions
//!
void
print_statistics ();

//!
//! \ingroup  BLAS_Module
//! \brief reset statistics for Algebra functions
//!
void
reset_statistics ();

//! \}

}// namespace BLAS

}// namespace HLIB

#endif  // __HLIB_BLAS_ALGEBRA_HH
