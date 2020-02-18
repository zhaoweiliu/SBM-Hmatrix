#ifndef __HLIB_TJACOBI_HH
#define __HLIB_TJACOBI_HH
//
// Project     : HLib
// File        : TJacobi.hh
// Description : class representing the Jacobi iteration matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/algebra/solve_types.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TJacobi );

//!
//! \ingroup Matrix_Module
//! \class   TJacobi
//! \brief   implements Jacobi preconditioner
//!
//!          TJacobi provides application of a Jacobi type preconditioner for a given
//!          matrix A, e.g. \f$ D^{-1} \f$ with \f$D\f$ being the diagonal of A. Here, 
//!          A may be of any type, e.g. dense, low-rank, block or sparse matrix. The
//!          diagonal elements are assumed to be non-zero.
//!
//!          Also support block-wise mode, i.e., diagonal blocks are inverted. For this,
//!          the size of diagonal block is defined during construction. If this equals 1
//!          point-wise Jacobi is used. Please note, that for block-wise mode, the leaf
//!          blocks define the smallest block size.
//!
class TJacobi : public TLinearOperator
{
private:
    //! @cond
    
    // the considered matrix
    const TMatrix *             _mat_A;

    // damping factor
    real                        _damping;

    // point- or blockwise
    eval_type_t                 _eval_type;

    // block diagonal inverse (internal use)
    std::unique_ptr< TMatrix >  _blockdiag;
    
    //! @endcond
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct Jacobi preconditioner based on matrix \a A
    //! with damping factor \a damping
    TJacobi ( const TMatrix *    A,
              const real         damping    = 1.0,
              const size_t       block_size = 1,
              const TTruncAcc &  acc_block  = acc_machine );

    virtual ~TJacobi ();
    
    /////////////////////////////////////////////////
    //
    // access internal data
    //

    //! return sparse matrix
    const TMatrix *  matrix  () const  { return _mat_A; }

    //! return damping factor
    real  damping_factor     () const  { return _damping; }
                                       
    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    bool  is_complex      () const
    {
        return _mat_A->is_complex();
    }
    
    //! return true, of operator is self adjoint
    bool  is_self_adjoint () const
    {
        return _mat_A->is_hermitian();
    }
    
    /////////////////////////////////////////////////
    //
    // linear operator mapping
    //

    //!
    //! mapping function of linear operator \f$A\f$, e.g. \f$ y := A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply       ( const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;

    //!
    //! mapping function with update: \f$ y := y + \alpha A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply_add   ( const real       alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;
    virtual void  capply_add  ( const complex    alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;

    virtual void  apply_add   ( const real       alpha,
                                const TMatrix *  X,
                                TMatrix *        Y,
                                const matop_t    op = apply_normal ) const;
    
    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const real                       alpha,
                                const BLAS::Vector< real > &     x,
                                BLAS::Vector< real > &           y,
                                const matop_t                    op = apply_normal ) const;
    virtual void  apply_add   ( const complex                    alpha,
                                const BLAS::Vector< complex > &  x,
                                BLAS::Vector< complex > &        y,
                                const matop_t                    op = apply_normal ) const;

    ///////////////////////////////////////////////////////////
    //
    // access to vector space elements
    //

    //! return dimension of domain
    virtual size_t  domain_dim () const { return _mat_A->range_dim(); }
    
    //! return dimension of range
    virtual size_t  range_dim  () const { return _mat_A->domain_dim(); }

    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector > { return _mat_A->row_vector(); }

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector > { return _mat_A->col_vector(); }
    
    HLIB_RTTI_DERIVED( TJacobi, TLinearOperator );
};

}// namespace HLIB

#endif  // __HLIB_TJACOBI_HH
