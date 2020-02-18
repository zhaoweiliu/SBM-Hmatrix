#ifndef __HLIB_TPERMMATRIX_HH
#define __HLIB_TPERMMATRIX_HH
//
// Project     : HLib
// File        : TPermMatrix.hh
// Description : class for a permuted matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <utility>

#include "hpro/base/types.hh"

#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TPermMatrix );

//!
//! implements a permuted matrix, e.g., P·A·R
//!
//! For a given, general matrix A and two permutations R, P, the product
//! \f$P \cdot A \cdot R\f$ is represented. If P/R is nullptr, it is treated
//! as identity.
//!
class TPermMatrix : public TLinearOperator
{
private:
    //! @cond

    const TPermutation *     _P;
    const TLinearOperator *  _A;
    const TPermutation *     _R;

    // unique pointers in case of ownership transfer
    std::unique_ptr< TPermutation >     _P_ptr;
    std::unique_ptr< TLinearOperator >  _A_ptr;
    std::unique_ptr< TPermutation >     _R_ptr;
    
    //! @endcond

private:

    // disable default constructor
    TPermMatrix ();
    
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    template < typename T_perm_P,
               typename T_linop,
               typename T_perm_R >
    TPermMatrix ( T_perm_P &&  P,
                  T_linop  &&  A,
                  T_perm_R &&  R )
    {
        set_perm_P( std::forward< T_perm_P >( P ) );
        set_linop(  std::forward< T_linop  >( A ) );
        set_perm_R( std::forward< T_perm_R >( R ) );

        if ( _A == nullptr )
            HERROR( ERR_ARG, "(TPermMatrix) ctor", "A is nullptr" );
    }

    virtual ~TPermMatrix () {}

    ///////////////////////////////////////////////////////////
    //
    // access to internal data
    //

    const TPermutation *     perm_P () const { return _P; }
    const TLinearOperator *  matrix () const { return _A; }
    const TPermutation *     perm_R () const { return _R; }

    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    virtual bool  is_complex      () const { return _A->is_complex(); }
    
    //! return true, of operator is self adjoint
    virtual bool  is_self_adjoint () const { return _A->is_self_adjoint(); }
    
    ///////////////////////////////////////////////////////////
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
    virtual size_t  domain_dim () const { return _A->domain_dim(); }
    
    //! return dimension of range
    virtual size_t  range_dim  () const { return _A->range_dim(); }

    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector > { return _A->domain_vector(); }

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector > { return _A->range_vector(); }

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TPermMatrix, TLinearOperator );

private:
    ///////////////////////////////////////////////////////////
    //
    // for perfect forwarding
    //

    // WITHOUT ownership transfer
    void  set_perm_P ( const TPermutation *     P ) { _P = P; }
    void  set_linop  ( const TLinearOperator *  A ) { _A = A; }
    void  set_perm_R ( const TPermutation *     R ) { _R = R; }


    // WITH ownership transfer
    void  set_perm_P ( std::unique_ptr< TPermutation > &&  P )
    {
        _P_ptr = std::move( P );
        _P     = _P_ptr.get();
    }
    void  set_linop  ( std::unique_ptr< TLinearOperator > &&  A )
    {
        _A_ptr = std::move( A );
        _A     = _A_ptr.get();
    }
    void  set_perm_R ( std::unique_ptr< TPermutation > &&  R )
    {
        _R_ptr = std::move( R );
        _R     = _R_ptr.get();
    }
};

}// namespace HLIB

#endif  // __HLIB_TPERMMATRIX_HH
