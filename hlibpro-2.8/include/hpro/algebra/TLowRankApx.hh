#ifndef __HLIB_TLOWRANKAPX_HH
#define __HLIB_TLOWRANKAPX_HH
//
// Project     : HLib
// File        : TLowRankApx.hh
// Description : classes for computing a low rank approximation of a dense matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <deque>
#include <atomic>

#include "hpro/base/error.hh"
#include "hpro/base/TPoint.hh"
#include "hpro/blas/Matrix.hh"
#include "hpro/cluster/TCluster.hh"
#include "hpro/cluster/TGeomCluster.hh"
#include "hpro/matrix/TCoeffFn.hh"
#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TRkMatrix.hh"
#include "hpro/parallel/TMutex.hh"

namespace HLIB
{
    
//!
//! \{
//! \name Low Rank Approximation
//!       Classes and functions related to low rank approximation of dense matrices.
//!

//!
//! \ingroup Algebra_Module
//! \class   TLowRankApx
//! \brief   base class for all low rank approximation techniques
//!
class TLowRankApx
{
public:
    //
    // statistical data of computations
    //
    struct stat_t
    {
        // number of requested matrix coefficients
        std::atomic< size_t >  ncoeff;

        // ctor nullifying data
        stat_t ()
                : ncoeff(0)
        {}
    };
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TLowRankApx () {}

    virtual ~TLowRankApx () {}

    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block cluster \a bct with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockCluster *   bct,
                              const TTruncAcc &       acc ) const
    {
        if ( bct == nullptr )
            HERROR( ERR_ARG, "(TLowRankApx) build", "block cluster is null" );
        
        const auto  rowis = * ( bct->rowcl() );
        const auto  colis = * ( bct->colcl() );

        return build( bis( rowis, colis ), acc );
    }

    //! build low rank matrix for block index set \a bis with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockIndexSet &  bis,
                              const TTruncAcc &       acc ) const = 0;

    //! indicate if algorithm provides statistics
    virtual bool  has_statistics () const { return false; }
};

//!
//! \ingroup Algebra_Module
//! \class   TZeroLRApx
//! \brief   Approximate all low-rank blocks by zero, e.g. for nearfield only.
//!
//!          If only the near field blocks of a matrix should be computed,
//!          TZeroLRApx will do that by approximating all far field blocks
//!          by zero, e.g. a 0-rank low rank matrix.
//!
class TZeroLRApx : public TLowRankApx
{
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TZeroLRApx () {}

    virtual ~TZeroLRApx () {}
    
    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! return low rank matrix of rank 0
    virtual TMatrix * build ( const TBlockIndexSet &  bis,
                              const TTruncAcc &       acc ) const;
    using TLowRankApx::build;
};
    
//!
//! \ingroup Algebra_Module
//! \class   TDenseLRApx
//! \brief   Computes dense matrix block without approximation.
//!
//!          Instead of performing approximation for a matrix block, the whole block
//!          is computed and returned as a dense matrix.
//!
//!          This is usually used for debugging or accuracy tests.
//!
template < typename T >
class TDenseLRApx : public TLowRankApx
{
public:
    //
    // template type as public member type
    //

    using  value_t    = T;
    using  real_t     = typename real_type< value_t >::type_t;
    using  coeff_fn_t = TCoeffFn< value_t >;
    
protected:
    // function return coefficient for index-pair
    const coeff_fn_t *  _coeff_fn;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TDenseLRApx ( const coeff_fn_t *  coeff_fn );

    virtual ~TDenseLRApx () {}
    
    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block index set \a bis with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockIndexSet &  bis,
                              const TTruncAcc &       acc ) const;
    using TLowRankApx::build;
};

//!
//! \ingroup Algebra_Module
//! \class   TSVDLRApx
//! \brief   Uses exact SVD to compute low rank approximation (WARNING: O(n³) complexity)
//!
//!          TSVDLRApx uses singular value decomposition to approximate a given
//!          matrix block. The resulting low rank matrix is the best approximation
//!          with respect to the given accuracy and rank. However, the computational
//!          costs are cubic in the dimension of the matrix block.
//!
template < typename T >
class TSVDLRApx : public TLowRankApx
{
public:
    //
    // template type as public member type
    //

    using  value_t    = T;
    using  real_t     = typename real_type< value_t >::type_t;
    using  coeff_fn_t = TCoeffFn< value_t >;
    
protected:
    // function return coefficient for index-pair
    const coeff_fn_t *  _coeff_fn;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TSVDLRApx ( const coeff_fn_t *  coeff_fn );

    virtual ~TSVDLRApx () {}
    
    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block cluster \a bcl with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockCluster *   bcl,
                              const TTruncAcc &       acc ) const;

    //! build low rank matrix for block index set \a block_is with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockIndexSet &  block_is,
                              const TTruncAcc &       acc ) const;

};

//!
//! \ingroup Algebra_Module
//! \class   TRandSVDLRApx
//! \brief   Uses randomized SVD to compute low rank approximation (WARNING: O(n²) complexity)
//!
//!          TRandSVDLRApx uses randomized singular value decomposition to approximate 
//!          a given matrix block. For the approximation, the complete matrix block has
//!          to be evaluated, hence complexity is O(n²).
//!
template < typename T >
class TRandSVDLRApx : public TLowRankApx
{
public:
    //! @cond
    
    //
    // template type as public member type
    //

    using  value_t    = T;
    using  real_t     = typename real_type< value_t >::type_t;
    using  coeff_fn_t = TCoeffFn< value_t >;
    
private:
    // function return coefficient for index-pair
    const coeff_fn_t *  _coeff_fn;

    // number of power iteration steps
    const uint          _power_steps;
    
    // oversampling for randomized approximation
    const uint          _oversampling;
    
    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TRandSVDLRApx ( const coeff_fn_t *  coeff_fn,
                    const uint          power_steps  = 0,
                    const uint          oversampling = 0 );

    virtual ~TRandSVDLRApx () {}
    
    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block cluster \a bcl with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockCluster *   bcl,
                              const TTruncAcc &       acc ) const;

    //! build low rank matrix for block index set \a block_is with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockIndexSet &  block_is,
                              const TTruncAcc &       acc ) const;

};

//!
//! \ingroup Algebra_Module
//! \class   TRRQRLRApx
//! \brief   Uses rank-revealing QR to compute low rank approximation
//!
//!          TRRQRLRApx uses rank-revealing QR to approximate a given matrix block.
//!          For the approximation, the complete matrix block has to be evaluated,
//!          hence complexity is at least O(n²).
//!
template < typename T >
class TRRQRLRApx : public TLowRankApx
{
public:
    //! @cond
    
    //
    // template type as public member type
    //

    using  value_t    = T;
    using  real_t     = typename real_type< value_t >::type_t;
    using  coeff_fn_t = TCoeffFn< value_t >;
    
private:
    // function return coefficient for index-pair
    const coeff_fn_t *  _coeff_fn;

    //! @endcond
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TRRQRLRApx ( const coeff_fn_t *  coeff_fn );

    virtual ~TRRQRLRApx () {}
    
    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block cluster \a bcl with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockCluster *   bcl,
                              const TTruncAcc &       acc ) const;

    //! build low rank matrix for block index set \a block_is with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockIndexSet &  block_is,
                              const TTruncAcc &       acc ) const;

};

//!
//! \ingroup Algebra_Module
//! \class   TACA
//! \brief   Defines interface for all ACA algorithms and implements classical ACA.
//!
//!          Adaptive cross approximation (ACA) is a heuristic for computing a
//!          low rank approximation of a given dense matrix by successively
//!          removing specific pairs of rows and columns (crosses) from the
//!          matrix until the rest is below some threshold (defined by block-wise
//!          accuracy).
//!
//!          Due to the algorithm, only the matrix coefficients in form of a
//!          TCoeffFn are needed, permitting the straightforward adaption of existing 
//!          implementations for the construction of \mcH-matrices.
//!
//!          The costs are linear in the dimension of the block and quadratic in the rank. 
//!
template < typename T >
class TACA : public TLowRankApx
{
public:
    //
    // template type as public member type
    //

    using  value_t    = T;
    using  real_t     = typename real_type< value_t >::type_t;

    //
    // internally used types
    //
    
    using  coeff_fn_t   = TCoeffFn< value_t >;
    using  vec_list_t   = std::deque< BLAS::Vector< value_t > >;
    using  idx_pair_t   = std::pair< idx_t, idx_t >;
    using  pivot_list_t = std::deque< idx_pair_t >;
    using  stat_t       = TLowRankApx::stat_t;
    
protected:
    // function return coefficient for index-pair
    const coeff_fn_t *  _coeff_fn;

    // statistics
    stat_t *            _stat;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TACA ( const coeff_fn_t *  coeff_fn );

    virtual ~TACA () {}
    
    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block cluster \a bct with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockCluster *   cl,
                              const TTruncAcc &       acc ) const;

    //! build low rank matrix for block index set \a block_is with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockIndexSet &  block_is,
                              const TTruncAcc &       acc ) const;

protected:
    // wrap ACA algorithm and convert result to matrix
    virtual TMatrix * aca_build ( const TBlockIndexSet &  block_is,
                                  const TTruncAcc &       acc ) const;
    
    // actual ACA algorithm
    virtual bool aca  ( const TIndexSet &  rowcl,
                        const TIndexSet &  colcl,
                        vec_list_t &       A,
                        vec_list_t &       B,
                        pivot_list_t &     pivots,
                        const size_t       rank,
                        const real_t       eps,
                        vec_list_t *       rows_orig,
                        vec_list_t *       cols_orig ) const;

public:
    // compute single corrected row of matrix
    void compute_row  ( const TIndexSet &          rowcl,
                        const TIndexSet &          colcl,
                        const vec_list_t &         A,
                        const vec_list_t &         B,
                        const idx_t                i,
                        BLAS::Vector< value_t > &  row,
                        BLAS::Vector< value_t > *  row_orig = nullptr ) const;
    
    // compute single corrected column of matrix
    void compute_col  ( const TIndexSet &          rowcl,
                        const TIndexSet &          colcl,
                        const vec_list_t &         A,
                        const vec_list_t &         B,
                        const idx_t                j,
                        BLAS::Vector< value_t > &  col,
                        BLAS::Vector< value_t > *  col_orig = nullptr ) const;

    // set statistics object
    void  set_stat  ( stat_t *  stat ) { _stat = stat; }
    
    //! indicate if algorithm provides statistics
    virtual bool  has_statistics () const { return true; }
};

//!
//! \ingroup Algebra_Module
//! \class   TACAPlus
//! \brief   Implements ACA+, which corrects some of the deficits of the original ACA algorithm.
//!
//!          The crucial point in the ACA algorithm is the search for suitable crosses (pairs of
//!          rows and columns) to be subtracted from the matrix. To keep costs linear in the
//!          dimension of the block, only some matrix entries may be checked, which may lead
//!          to a breakdown of the algorithm before the user defined accuracy was reached.
//!
//!          ACA+ implements a more sophisticated pivot search, which eliminates some of the
//!          cases where the original ACA fails to compute adequate low-rank approximations.
//!          Nevertheless, the computational costs are only slightly higher than ACA.
//!
template < typename T >
class TACAPlus : public TACA< T >
{
public:
    //
    // template type as public member type
    //

    using  value_t    = T;
    using  real_t     = typename real_type< value_t >::type_t;

    //
    // derived types
    //
    
    using  coeff_fn_t   = typename TACA< value_t >::coeff_fn_t;
    using  vec_list_t   = typename TACA< value_t >::vec_list_t;
    using  idx_pair_t   = typename TACA< value_t >::idx_pair_t;
    using  pivot_list_t = typename TACA< value_t >::pivot_list_t;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TACAPlus ( const coeff_fn_t *  coeff_fn )
            : TACA< value_t >( coeff_fn )
    {}

    virtual ~TACAPlus () {}
    
protected:
    // actual ACA algorithm
    virtual bool  aca  ( const TIndexSet &  rowcl,
                         const TIndexSet &  colcl,
                         vec_list_t &       A,
                         vec_list_t &       B,
                         pivot_list_t &     pivots,
                         const size_t       rank,
                         const real_t       eps,
                         vec_list_t *       rows_orig,
                         vec_list_t *       cols_orig ) const;
};

//!
//! \ingroup Algebra_Module
//! \class   TACAFull
//! \brief   ACA with full pivot search (complexity: O(n²))
//!
//!          Full ACA tests <em>all</em> matrix coefficients in the search for the
//!          best cross. This results in a guaranteed approximation within the given
//!          accuracy. However, the costs are now quadratic in the size of the block.
//!
template < typename T >
class TACAFull : public TACA< T >
{
public:
    //
    // template type as public member type
    //

    using  value_t    = T;
    using  real_t     = typename real_type< value_t >::type_t;

    //
    // derived types
    //
    
    using  coeff_fn_t   = typename TACA< value_t >::coeff_fn_t;
    using  vec_list_t   = typename TACA< value_t >::vec_list_t;
    using  idx_pair_t   = typename TACA< value_t >::idx_pair_t;
    using  pivot_list_t = typename TACA< value_t >::pivot_list_t;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TACAFull ( const coeff_fn_t *  coeff_fn )
            : TACA< value_t >( coeff_fn )
    {}

    virtual ~TACAFull () {}
    
protected:
    // actual ACA algorithm
    virtual bool  aca  ( const TIndexSet &  rowcl,
                         const TIndexSet &  colcl,
                         vec_list_t &       A,
                         vec_list_t &       B,
                         pivot_list_t &     pivots,
                         const size_t       rank,
                         const real_t       eps,
                         vec_list_t *       rows_orig,
                         vec_list_t *       cols_orig ) const;
};

//!
//! \ingroup Algebra_Module
//! \class   THCA
//! \brief   uses hybrid cross approximation (HCA) for computing low rank approximation
//!
//!          THCA provides a low rank approximation algorithm with a guaranteed
//!          approximation quality. It is based on the <em>generator function</em> \f$\gamma(x,y)\f$
//!          of a BEM kernel function \f$k(x,y)\f$ and it's derivatives \f$ D_x \gamma(x, y_{l}) \f$
//!          and \f$ D_y \gamma(x, y_{l}) \f$.
//!
//!          The class THCA needs a user implemented generator function of type TGeneratorFn, in which
//!          the function itself and the integrals of the corresponding derivates are defined.
//!
//!          Furthermore, HCA is based on interpolation, of which the order defines the accuracy of
//!          the final result. This interpolation order together with an accuracy for the approximation
//!          of the generator function is specific to the given problem, and hence, user defined.
//!
template < typename T >
class THCA : public TLowRankApx
{
public:
    //
    // template arguments as internal types
    //

    using  value_t    = T;
    using  real_t     = typename real_type< value_t >::type_t;
    using  coeff_fn_t = TCoeffFn< value_t >;

    //!
    //! \class  TGeneratorFn
    //! \brief  class defining kernel generator function used by HCA
    //!
    class TGeneratorFn
    {
    public:
        //
        // constructor and destructor
        //

        TGeneratorFn () {}
        
        virtual ~TGeneratorFn () {}

        //! indicate complex nature of function
        virtual bool     is_complex    () const { return is_complex_type< value_t >::value; }
        
        //!
        //! Evaluate generator function γ at (\a x, \a y).
        //!
        virtual value_t  eval          ( const T3Point &  x,
                                         const T3Point &  y ) const = 0;

        //!
        //! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
        //! and points \f$ y_{l} \f$ defined by \a pts.
        //! Store results in \a matrix at index (i,l).
        //!
        virtual void     integrate_dx  ( const TIndexSet &               is,
                                         const std::vector< T3Point > &  pts,
                                         BLAS::Matrix< value_t > &       matrix ) const = 0;
        //!
        //! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
        //! and points \f$ x_{l} \f$ defined by \a pts. 
        //! Store results in \a matrix at index (j,l).
        //!
        virtual void     integrate_dy  ( const TIndexSet &               is,
                                         const std::vector< T3Point > &  pts,
                                         BLAS::Matrix< value_t > &       matrix ) const = 0;
    };

    //!
    //! statistics
    //! - timing will only work if enabled internally
    //!
    struct stat_t : public TLowRankApx::stat_t
    {
        //! time for pivot computation
        double  time_pivot;

        //! time for fall back in case of "full rank"
        double  time_fallback;

        //! time for computing U and V matrices
        double  time_UV;

        //! indicates, whether timing data is collected
        bool    has_timing;

        //! ctor nullifying everything
        stat_t ()
                : TLowRankApx::stat_t()
                , time_pivot( 0 )
                , time_fallback( 0 )
                , time_UV( 0 )
                , has_timing( false )
        {}
    };
    
protected:
    //
    // multiindex of dimension three
    //
    struct  idx3_t
    {
        idx_t  idxs[3];

        idx_t &  operator []  ( idx_t  i ) { return idxs[i]; }
    };
    
    //
    // store tensor grid, e.g. grid defined by (x_i,y_j,z_k)
    //
    struct  tensor_grid_t
    {
        std::vector< double >  x, y, z;

        //! return point at index (i,j,k)
        T3Point  grid_point ( const idx_t  i,
                              const idx_t  j,
                              const idx_t  k ) const
        {
            return T3Point( x[i], y[j], z[k] );
        }

        //! return point at (multi-) index \a midx = (i,j,k)
        T3Point  grid_point ( const idx3_t  midx ) const
        {
            return T3Point( x[ midx.idxs[0] ], y[ midx.idxs[1] ], z[ midx.idxs[2] ] );
        }
    };
    
protected:
    // interpolation order to use
    uint                  _ipol_order;

    // unit interval interpolation points
    std::vector< std::vector< double > > _ipol_points;

    // function for matrix evaluation
    const coeff_fn_t *    _coeff;
    
    // function for kernel evaluation
    const TGeneratorFn *  _generator_fn;

    // accuracy of ACA inside HCA
    const real_t          _aca_eps;

    // statistics
    stat_t *              _stat;
    
public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct HCA based on matrix coefficient function \a coeff
    //! and kernel generator function \a generator.
    //!
    //! The accuracy of the HCA algorithm is defined by \a aca_eps,
    //! the relative accuracy of the internal ACA algorithm and the
    //! interpolation order \a ipol_order.
    THCA ( const coeff_fn_t *    coeff,
           const TGeneratorFn *  generator,
           const real_t          aca_eps,
           const uint            ipol_order );

    virtual ~THCA () {}

    //////////////////////////////////////
    //
    // build low-rank matrix
    //

    //! build low rank matrix for block cluster \a bct with
    //! rank defined by accuracy \a acc
    virtual TMatrix * build ( const TBlockCluster *   bc,
                              const TTruncAcc &       acc ) const;

    virtual TMatrix * build ( const TBlockIndexSet & ,
                              const TTruncAcc & ) const
    {
        HERROR( ERR_NOT_IMPL, "(THCA) build", "" );
    }
    
    //////////////////////////////////////
    //
    // misc.
    //

    // set statistics object
    void  set_stat  ( stat_t *  stat ) { _stat = stat; }
    
    //! indicate if algorithm provides statistics
    virtual bool  has_statistics () const { return _stat != nullptr; }
    
protected:
    //
    // actual HCA algorithm
    //
    TMatrix *  hca         ( const TBlockCluster *   bc,
                             const TTruncAcc &       acc ) const;

    //
    // compute pivot points for generator function approximation
    //
    size_t     comp_pivot  ( const tensor_grid_t &    rowcl_grid,
                             const tensor_grid_t &    colcl_grid,
                             const real_t             eps,
                             std::vector< idx_t > &   row_pivots,
                             std::vector< idx_t > &   col_pivots,
                             const uint               ipol_order ) const;

    //
    // interpolation management
    //
    
    //! build interpolation points for unit cube in \f$ R^d \f$
    void  setup_interpolation    ();

    //! transform interpolation points to local cluster
    void  transform_ipol_points  ( const TGeomCluster *  cl,
                                   tensor_grid_t &       grid,
                                   const uint            ipol_order ) const;

    //
    // evaluate kernel generator matrix
    //
    
    //! compute single corrected row of kernel generator matrix
    void  compute_row            ( const idx_t                      row_idx,
                                   const tensor_grid_t &            rowcl_grid,
                                   const tensor_grid_t &            colcl_grid,
                                   const BLAS::Matrix< value_t > &  A,
                                   const BLAS::Matrix< value_t > &  B,
                                   const size_t                     rank,
                                   BLAS::Vector< value_t > &        row,
                                   const uint                       ipol_order ) const;

    //! compute single corrected row of kernel generator matrix
    void  compute_col            ( const idx_t                      col_idx,
                                   const tensor_grid_t &            rowcl_grid,
                                   const tensor_grid_t &            colcl_grid,
                                   const BLAS::Matrix< value_t > &  A,
                                   const BLAS::Matrix< value_t > &  B,
                                   const size_t                     rank,
                                   BLAS::Vector< value_t > &        col,
                                   const uint                       ipol_order ) const;

    //
    // compute collocation matrices
    //
    
    //! compute matrix U
    void  compute_U              ( const TIndexSet &             rowis,
                                   const size_t                  rank,
                                   BLAS::Matrix< value_t > &     U,
                                   const std::vector< idx_t > &  col_pivot,
                                   const tensor_grid_t &         colcl_grid,
                                   const uint                    ipol_order ) const;
    
    //! compute matrix V
    void  compute_V              ( const TIndexSet &             colis,
                                   const size_t                  rank,
                                   BLAS::Matrix< value_t > &     V,
                                   const std::vector< idx_t > &  row_pivot,
                                   const tensor_grid_t &         rowcl_grid,
                                   const uint                    ipol_order ) const;

    DISABLE_COPY_OP( THCA );
};

//!
//! \ingroup Algebra_Module
//! \class   TPermHCAGeneratorFn
//! \brief   base class for HCA generator functions using row/column permutations
//!
//!          Provides basic permutation management for evaluating the integrals
//!          over the derivatives of the kernel generator function.
//!
template < typename T_val >
class TPermHCAGeneratorFn : public THCA< T_val >::TGeneratorFn
{
public:
    //
    // template types as internal types
    //
    using  value_t = T_val;

protected:
    
    // mapping from int. to ext. numbering
    const TPermutation *  _row_perm_i2e;
    const TPermutation *  _col_perm_i2e;
    
public:
    //!
    //! constructor
    //! - \a order defines (maximal) quadrature order for
    //!   evaluating the integrals 
    //!
    TPermHCAGeneratorFn ( const TPermutation *  row_perm_i2e,
                          const TPermutation *  col_perm_i2e )
            : _row_perm_i2e( row_perm_i2e ),
              _col_perm_i2e( col_perm_i2e )
    {}
    
    //!
    //! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
    //! and points \f$ y_{l} \f$ defined by \a pts.
    //! Store results in \a matrix at index (i,l).
    //!
    virtual void  integrate_dx  ( const TIndexSet &               is,
                                  const std::vector< T3Point > &  pts,
                                  BLAS::Matrix< value_t > &       matrix ) const
    {
        //
        // permute indices
        //

        const bool            has_perm  = (this->_row_perm_i2e != NULL);
        std::vector< idx_t >  idxs( is.size() );
        size_t                pos = 0;

        for ( auto idx : is )
        {
            const idx_t  ex_idx = (has_perm ? this->_row_perm_i2e->permute( idx ) : idx );

            idxs[ pos++ ] = ex_idx;
        }// for

        integrate_dx_perm( idxs, pts, matrix );
    }
    
    //!
    //! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
    //! and points \f$ x_{l} \f$ defined by \a pts. 
    //! Store results in \a matrix at index (j,l).
    //!
    virtual void  integrate_dy  ( const TIndexSet &               is,
                                  const std::vector< T3Point > &  pts,
                                  BLAS::Matrix< value_t > &       matrix ) const
    {
        //
        // permute indices
        //

        const bool            has_perm  = (this->_col_perm_i2e != NULL);
        std::vector< idx_t >  idxs( is.size() );
        size_t                pos = 0;

        for ( auto idx : is )
        {
            const idx_t  ex_idx = (has_perm ? this->_col_perm_i2e->permute( idx ) : idx );

            idxs[ pos++ ] = ex_idx;
        }// for

        integrate_dy_perm( idxs, pts, matrix );
    }

    //!
    //! Evaluate \f$ \int \phi_i(x) D_x \gamma(x, y_{l}) dx\f$ for i ∈ \a is
    //! and points \f$ y_{l} \f$ defined by \a pts.
    //! Indices in \a idxs are given in external order.
    //!
    virtual void  integrate_dx_perm  ( const std::vector< idx_t > &    idxs,
                                       const std::vector< T3Point > &  pts,
                                       BLAS::Matrix< value_t > &       matrix ) const = 0;
    //!
    //! Evaluate \f$ \int \phi_j(y) D_y \gamma(x_{l}, y) dy\f$ for j ∈ \a is 
    //! and points \f$ x_{l} \f$ defined by \a pts. 
    //! Indices in \a idxs are given in external order.
    //!
    virtual void  integrate_dy_perm  ( const std::vector< idx_t > &    idxs,
                                       const std::vector< T3Point > &  pts,
                                       BLAS::Matrix< value_t > &       matrix ) const = 0;
};

//! \}

}// namespace HLIB

#endif  // __HLIB_TLOWRANKAPX_HH
