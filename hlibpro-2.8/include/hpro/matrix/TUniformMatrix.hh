#ifndef __HLIB_TUNIFORMMATRIX_HH
#define __HLIB_TUNIFORMMATRIX_HH
//
// Project     : HLib
// File        : TUniformMatrix.hh
// Description : class for uniform lowrank matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/blas/Matrix.hh"
#include "hpro/cluster/TClusterBasis.hh"
#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TUniformMatrix );

//!
//! \ingroup Matrix_Module
//! \class   TUniformMatrix
//! \brief   Represents low rank matrices as uniform matrix: \f$M = V S W^H\f$,
//!          where \f$V\f$ and \f$W\f$ are cluster bases and \f$S\f$ holds the 
//!          corresponding coefficients.
//!
class TUniformMatrix : public TMatrix
{
protected:
    // coefficient matrices for real and complex storage
    BLAS::Matrix< real >              _rcoeff;
    BLAS::Matrix< complex >           _ccoeff;

    // row and column cluster bases, real
    const TClusterBasis< real > *     _rrow_cb;
    const TClusterBasis< real > *     _rcol_cb;
    // and complex valued
    const TClusterBasis< complex > *  _crow_cb;
    const TClusterBasis< complex > *  _ccol_cb;
    
    // number of rows and columns
    size_t                            _rows, _cols;
    
    // row and column rank
    size_t                            _row_rank, _col_rank;
    
public:
    /////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! default constructor
    TUniformMatrix ();

    //! construct uniform matrix over \a block_is with
    //! clusterbases \a row_cb and \a col_cb and coefficients
    //! stored in \a S
    TUniformMatrix ( const TBlockIndexSet &            block_is,
                     const TClusterBasis< real > *     row_cb,
                     const TClusterBasis< real > *     col_cb,
                     const BLAS::Matrix< real > &      S );

    //! construct uniform matrix over \a block_is with
    //! clusterbases \a row_cb and \a col_cb and coefficients
    //! stored in \a S
    TUniformMatrix ( const TBlockIndexSet &            block_is,
                     const TClusterBasis< complex > *  row_cb,
                     const TClusterBasis< complex > *  col_cb,
                     const BLAS::Matrix< complex > &   S );

    //! copy constructor
    TUniformMatrix ( const TUniformMatrix &            A );

    //! dtor
    ~TUniformMatrix ()
    {}

    /////////////////////////////////////////////////
    //
    // access local variables
    //

    //
    // access size information
    //
    
    //! return number of rows 
    virtual size_t  rows         () const { return _rows; }

    //! return number of columns 
    virtual size_t  cols         () const { return _cols; }

    //
    // change size
    //
    
    //! set size of matrix
    void            set_size     ( const size_t  n,
                                   const size_t  m )
    {
        _rows = n;
        _cols = m;
    }

    //
    // rank and cluster bases
    //
    
    //! return row rank of matrix
    size_t          row_rank     () const { return _row_rank; }

    //! return column rank of matrix
    size_t          col_rank     () const { return _col_rank; }

    //! return row cluster basis
    const TClusterBasis< real > *     rrow_cb () const { return _rrow_cb; }

    //! return column cluster basis
    const TClusterBasis< real > *     rcol_cb () const { return _rcol_cb; }
        
    //! return row cluster basis
    const TClusterBasis< complex > *  crow_cb () const { return _crow_cb; }

    //! return column cluster basis
    const TClusterBasis< complex > *  ccol_cb () const { return _ccol_cb; }
        
    //! set ranks of matrix
    //! - has to be identical with current or future cluster bases rank!
    void            set_rank     ( const size_t  row_rank,
                                   const size_t  col_rank );

    //! assign cluster bases
    //! - rank and value type of bases must be identical to corresponding
    //!   dimension and value type of coefficient matrix
    void            assign_cb    ( const TClusterBasis< real > *     row_cb,
                                   const TClusterBasis< real > *     col_cb );
    void            assign_cb    ( const TClusterBasis< complex > *  row_cb,
                                   const TClusterBasis< complex > *  col_cb );
    
    //
    // access matrix coefficients
    //
    
    //! return real valued coefficent (\a i, \a j) of matrix
    real            entry        ( const idx_t  i, const idx_t j ) const;

    //! return complex valued coefficent (\a i, \a j) of matrix
    const complex   centry       ( const idx_t  i, const idx_t j ) const;

    //! return basis coefficients
    const BLAS::Matrix< real > &     rcoeff () const { return _rcoeff; }
    BLAS::Matrix< real > &           rcoeff ()       { return _rcoeff; }
    
    //! return basis coefficients
    const BLAS::Matrix< complex > &  ccoeff () const { return _ccoeff; }
    BLAS::Matrix< complex > &        ccoeff ()       { return _ccoeff; }
    
    //
    // change field type 
    //
    
    //! switch to real valued storage if possible, e.g. all imaginary components zero
    virtual void    to_real      ();

    //! switch to complex valued storage
    virtual void    to_complex   ();
    
    /////////////////////////////////////////////////
    //
    // manage Rk-matrices
    //

    //! transpose matrix
    virtual void    transpose    ();
    
    //! conjugate matrix coefficients
    virtual void    conjugate    ();
    
    //! truncate matrix  w.r.t. accuracy \a acc
    virtual void    truncate     ( const TTruncAcc & acc );

    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! scale matrix by constant factor
    virtual void         scale      ( const real       f );
    
    //! compute y ≔ α·M·x + β·y, where M is either this, this^T or this^H depending on \a op
    virtual void         mul_vec    ( const real       alpha,
                                      const TVector *  x,
                                      const real       beta,
                                      TVector *        y,
                                      const matop_t    op = MATOP_NORM ) const;

    //! compute this = this + α·A without truncation
    virtual void         add        ( const real       alpha,
                                      const TMatrix *  A );
        
    //! compute matrix product α·op_A(this)·op_B(B) 
    virtual TMatrix *    mul_right  ( const real       alpha,
                                      const TMatrix *  B,
                                      const matop_t    op_A,
                                      const matop_t    op_B ) const;

    //! compute matrix product α·op_A(A)·op_B(this) 
    virtual TMatrix *    mul_left   ( const real       alpha,
                                      const TMatrix *  A,
                                      const matop_t    op_A,
                                      const matop_t    op_B ) const;

    /////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! scale matrix by constant factor
    virtual void        cscale      ( const complex    f );
    
    //! compute y ≔ α·M·x + β·y, where M is either this, this^T or this^H depending on \a op
    virtual void        cmul_vec    ( const complex    alpha,
                                      const TVector *  x,
                                      const complex    beta,
                                      TVector *        y,
                                      const matop_t    op_A = MATOP_NORM ) const;

    //! compute this = this + α·A without truncation
    virtual void        cadd        ( const complex    a,
                                      const TMatrix *  matrix );
        
    //! compute matrix product α·op_A(this)·op_B(B) 
    virtual TMatrix *   cmul_right  ( const complex    alpha,
                                      const TMatrix *  B,
                                      const matop_t    op_A,
                                      const matop_t    op_B ) const;
    
    //! compute matrix product α·op_A(A)·op_B(this) 
    virtual TMatrix *   cmul_left   ( const complex    alpha,
                                      const TMatrix *  A,
                                      const matop_t    op_A,
                                      const matop_t    op_B ) const;

    /////////////////////////////////////////////////
    //
    // misc methods
    //

    //
    // virtual constructors
    //

    //! return object of same type
    virtual auto  create       () const -> std::unique_ptr< TMatrix >
    {
        return std::make_unique< TUniformMatrix >();
    }

    //! return copy of matrix
    virtual auto  copy         () const -> std::unique_ptr< TMatrix >;

    //! return copy matrix wrt. given accuracy; if \a do_coarsen is set, perform coarsening
    virtual auto  copy         ( const TTruncAcc &  acc,
                                 const bool         do_coarsen = false ) const -> std::unique_ptr< TMatrix >;

    //! return structural copy of matrix
    virtual auto  copy_struct  () const -> std::unique_ptr< TMatrix >;

    // copy matrix data to \a A
    virtual void  copy_to      ( TMatrix *          A ) const;

    // copy matrix data to \a A and truncate w.r.t. \acc with optional coarsening
    virtual void  copy_to      ( TMatrix *          A,
                                 const TTruncAcc &  acc,
                                 const bool         do_coarsen = false ) const;
    
    //! return size in bytes used by this object
    virtual size_t byte_size   () const;

    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TUniformMatrix, TMatrix )

    //
    // serialisation
    //

    //! read matrix from byte stream
    virtual void read  ( TByteStream & s );

    //! construct matrix from byte stream
    virtual void build ( TByteStream & s );
    
    //! write matrix into byte stream
    virtual void write ( TByteStream & s ) const;

    //! return size of object in a bytestream
    virtual size_t bs_size () const;
};

}// namespace HLIB

#endif  // __HLIB_TUNIFORMMATRIX_HH
