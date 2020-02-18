#ifndef __HLIB_MAT_APPROX_HH
#define __HLIB_MAT_APPROX_HH
//
// Project     : HLib
// File        : mat_approx.hh
// Description : matrix approxition functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <list>
#include <utility>

#include "hpro/matrix/TLinearOperator.hh"
#include "hpro/matrix/TRkMatrix.hh"

namespace HLIB
{

//
// compute low-rank approximation of a given matrix sum
// Σ_i M_i using SVD
//
std::unique_ptr< TMatrix >
approx_svd      ( std::list< const TMatrix * > &          M,
                  const TTruncAcc &                       acc,
                  const bool                              seq = false );

//
// compute low-rank approximation of a sum Σ_i U_i V_i^H using SVD
// ASSUMPTIONS: - #cols(U_i) = #cols(V_i)
//              - #rows(U_i) = #rows(U_j)
//              - #rows(V_i) = #rows(V_j)
//
template< typename T >
std::pair< BLAS::Matrix< T >, BLAS::Matrix< T > >
approx_svd      ( const std::list< BLAS::Matrix< T > > &  U,
                  const std::list< BLAS::Matrix< T > > &  V,
                  const TTruncAcc &                       acc );

//
// same as bove but with faster pair-wise SVD not yielding best approximation
//
template< typename T >
std::pair< BLAS::Matrix< T >, BLAS::Matrix< T > >
approx_svd_pair ( const std::list< BLAS::Matrix< T > > &  U,
                  const std::list< BLAS::Matrix< T > > &  V,
                  const TTruncAcc &                       acc );

//
// compute low-rank approximation of a given sum Σ_i M_i using
// randomised SVD
// - only need matrix-vector evaluation of given operators
//
std::unique_ptr< TRkMatrix >
approx_randsvd  ( const TBlockIndexSet &                  bis,
                  const std::list< TLinearOperator * > &  M,
                  const TTruncAcc &                       acc );

template < typename container_t >
std::unique_ptr< TRkMatrix >
approx_randsvd  ( const TBlockIndexSet &                  bis,
                  const container_t &                     M,
                  const TTruncAcc &                       acc )
{
    std::list< TLinearOperator * >  linop_M;

    for ( auto  M_i : M )
        linop_M.push_back( M_i );

    return approx_randsvd( bis, linop_M, acc );
}

template < typename value_t >
std::pair< BLAS::Matrix< value_t >, BLAS::Matrix< value_t > >
approx_randsvd  ( const std::list< TLinearOperator * > &  M,
                  const TTruncAcc &                       acc );

template < typename value_t, typename container_t >
std::pair< BLAS::Matrix< value_t >, BLAS::Matrix< value_t > >
approx_randsvd  ( const container_t &                     M,
                  const TTruncAcc &                       acc )
{
    std::list< TLinearOperator * >  linop_M;

    for ( auto  M_i : M )
        linop_M.push_back( M_i );

    return approx_randsvd< value_t >( linop_M, acc );
}

template < typename value_t >
std::pair< BLAS::Matrix< value_t >, BLAS::Matrix< value_t > >
approx_randsvd  ( const BLAS::Matrix< value_t > &  M,
                  const TTruncAcc &                acc );

//
// compute randomized low-rank approximation of a given sum Σ_i M_i
// - direct randomised method without orthogonal column basis
// - only need matrix-vector evaluation of given operators
//
std::unique_ptr< TRkMatrix >
approx_randlr  ( const TBlockIndexSet &                  bis,
                 const std::list< TLinearOperator * > &  M,
                 const TTruncAcc &                       acc );

template < typename container_t >
std::unique_ptr< TRkMatrix >
approx_randlr  ( const TBlockIndexSet &                  bis,
                 const container_t &                     M,
                 const TTruncAcc &                       acc )
{
    std::list< TLinearOperator * >  linop_M;

    for ( auto  M_i : M )
        linop_M.push_back( M_i );

    return approx_randlr( bis, linop_M, acc );
}

//
// compute low-rank approximation of a sum Σ_i M_i using randomized LR
// without index sets (only block dimension of interest)
//
template < typename value_t >
std::pair< BLAS::Matrix< value_t >, BLAS::Matrix< value_t > >
approx_randlr  ( const std::list< TLinearOperator * > &  M,
                 const TTruncAcc &                       acc );

template < typename value_t, typename container_t >
std::pair< BLAS::Matrix< value_t >, BLAS::Matrix< value_t > >
approx_randlr  ( const container_t &                     M,
                 const TTruncAcc &                       acc )
{
    std::list< TLinearOperator * >  linop_M;

    for ( auto  M_i : M )
        linop_M.push_back( M_i );

    return approx_randlr< value_t >( linop_M, acc );
}

//
// compute low-rank approximation of a sum Σ_i U_i V_i^H using randomized LR
// ASSUMPTIONS: - #cols(U_i) = #cols(V_i)
//              - #rows(U_i) = #rows(U_j)
//              - #rows(V_i) = #rows(V_j)
//
template< typename T >
std::pair< BLAS::Matrix< T >, BLAS::Matrix< T > >
approx_randlr  ( const std::list< BLAS::Matrix< T > > &  U,
                 const std::list< BLAS::Matrix< T > > &  V,
                 const TTruncAcc &                       acc );

//
// compute low-rank approximation of a dense matrix M using randomized LR
//
template< typename T >
std::pair< BLAS::Matrix< T >, BLAS::Matrix< T > >
approx_randlr  ( const BLAS::Matrix< T > &  M,
                 const TTruncAcc &          acc );

//
// compute low-rank approximation U·V^H of a given sum Σ_i M_i using ACA
// - only matrix-vector evaluation of given operators is used
// - if \a pivots != nullptr, the pivot elements are stored in it
//
template < typename value_t >
std::pair< BLAS::Matrix< value_t >, BLAS::Matrix< value_t > >
approx_aca  ( const std::list< TLinearOperator * > &    M,
              const TTruncAcc &                         acc,
              std::list< std::pair< idx_t, idx_t > > *  pivots = nullptr );

template < typename value_t, typename container_t >
std::pair< BLAS::Matrix< value_t >, BLAS::Matrix< value_t > >
approx_aca  ( const container_t &                       M,
              const TTruncAcc &                         acc,
              std::list< std::pair< idx_t, idx_t > > *  pivots = nullptr )
{
    std::list< TLinearOperator * >  linop_M;

    for ( auto  M_i : M )
        linop_M.push_back( M_i );

    return approx_aca< value_t >( linop_M, acc, pivots );
}

}// namespace HLIB

#endif  // __HLIB_MAT_APPROX_HH
