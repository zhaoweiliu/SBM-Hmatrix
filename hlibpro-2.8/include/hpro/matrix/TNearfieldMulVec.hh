#ifndef __HLIB_TNEARFIELDMULVEC_HH
#define __HLIB_TNEARFIELDMULVEC_HH
//
// Project     : HLib
// File        : TNearfieldMulVec.hh
// Description : class performing matrix-vector multiplication with nearfield part
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/THMatrix.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TNearfieldMulVec );

//!
//! \ingroup  Matrix_Module
//! \class    TNearfieldMulVec
//! \brief    Implements matrix-vector multiplication with nearfield part of H-matrix.
//!
class TNearfieldMulVec : public TLinearOperator
{
private:
    // the referenced H matrix
    const THMatrix *  _matrix;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    TNearfieldMulVec ( const THMatrix *  A );

    virtual ~TNearfieldMulVec ();
    
    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    virtual bool  is_complex      () const { return _matrix->is_complex(); }
    
    //! return true, of operator is self adjoint
    virtual bool  is_self_adjoint () const { return _matrix->is_self_adjoint(); }
    
    /////////////////////////////////////////////////
    //
    // linear operator mapping
    //

    //!
    //! mapping function of linear operator A, e.g. y ≔ A(x).
    //! Depending on \a op, either A, A^T or A^H is applied.
    //!
    virtual void  apply       ( const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;

    //!
    //! mapping function with update: \a y ≔ \a y + \a α \a A( \a x ).
    //! Depending on \a op, either A, A^T or A^H is applied.
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
    virtual size_t  domain_dim () const { return _matrix->domain_dim(); }
    
    //! return dimension of range
    virtual size_t  range_dim  () const { return _matrix->range_dim(); }

    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector > { return _matrix->col_vector(); }

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector > { return _matrix->row_vector(); }


    HLIB_RTTI_DERIVED( TNearfieldMulVec, TLinearOperator );
};

}// namespace

#endif  // __HLIB_TNEARFIELDMULVEC_HH
