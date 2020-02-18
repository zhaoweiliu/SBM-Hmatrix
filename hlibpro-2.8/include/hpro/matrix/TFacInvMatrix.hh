#ifndef __HLIB_TFACINVMATRIX_HH
#define __HLIB_TFACINVMATRIX_HH
//
// Project     : HLib
// File        : TFacInvMatrix.hh
// Description : classes for representing the inverse of factorised matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <utility>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/algebra/solve_types.hh"

namespace HLIB
{

// local RTTI types
DECLARE_TYPE( TFacInvMatrix );
DECLARE_TYPE( TLUInvMatrix );
DECLARE_TYPE( TLDUInvMatrix );
DECLARE_TYPE( TLDLInvMatrix );
DECLARE_TYPE( TLLInvMatrix );

//////////////////////////////////////////////////////////////////////////////////
//
// TFacInvMatrix
//
//////////////////////////////////////////////////////////////////////////////////

// local matrix type
// DECLARE_TYPE( TFacInvMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TFacInvMatrix
//! \brief    Baseclass for representing the inverse of factorised matrices.
//!
//!           TFacInvMatrix and the derived classes represent operators based on
//!           factorised matrices providing the ability to evaluate the inverse
//!           of the corresponding matrix (in factorised form), e.g. for \f$ A = LU \f$
//!           the inverse \f$ (LU)^{-1} \f$ is applied.
//!
//!           Optionally, diagonal scalings from left and right may have been
//!           applied before the factorisation and now taken into account in the
//!           evaluation.
//!
class TFacInvMatrix : public TLinearOperator
{
private:
    //! @cond
    
    //! factorised matrix
    const TMatrix *             _fac_matrix;

    //! auto pointer in case of internal memory management
    std::unique_ptr< TMatrix >  _mat_ptr;
    
    //! indicates usage of diagonal scaling
    bool                        _has_scaling;
    
    //! diagonal scaling factors
    TScalarVector               _D_left, _D_right;

    //! determines pointwise of blockwise evaluation
    eval_type_t                 _eval_type;

    //! determines type of storage for diagonal blocks
    storage_type_t              _storage_type;

    //! signal if operator is self adjoint
    bool                        _self_adjoint;

    // accuracy for matrix arithmetic
    TTruncAcc                   _acc;
    
    //! @endcond
    
public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct inverse operator with factorised matrix \a afac_matrix
    //! (no diagonal scaling applied)
    template < typename T_mat >
    TFacInvMatrix ( T_mat &&               afac_matrix,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : _has_scaling( false )
            , _eval_type(    aeval_type )
            , _storage_type( astorage_type )
            , _acc( acc_exact )
    {
        set_matrix( std::forward< T_mat >( afac_matrix ) );
        
        if ( _fac_matrix == nullptr )
            HERROR( ERR_ARG, "(TFacInvMatrix)", "LU factor is nullptr" );
        
        init();
    }

    //! construct inverse operator with factorised matrix \a afac_matrix
    //! with accuracy for matrix arithmetic (without diagonal scaling)
    template < typename T_mat >
    TFacInvMatrix ( T_mat &&               afac_matrix,
                    const TTruncAcc &      aacc,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : _has_scaling( false )
            , _eval_type(    aeval_type )
            , _storage_type( astorage_type )
            , _acc( aacc )
    {
        set_matrix( std::forward< T_mat >( afac_matrix ) );
        
        if ( _fac_matrix == nullptr )
            HERROR( ERR_ARG, "(TFacInvMatrix)", "LU factor is nullptr" );
        
        init();
    }

    //! construct inverse operator with factorised matrix \a afac_matrix
    //! and additional row and column scaling factors \a D1 and \a D2,
    //! e.g. A = D1·Ã·D2, with Ã = \a afac_matrix
    template < typename T_mat >
    TFacInvMatrix ( const TScalarVector *  D1,
                    T_mat &&               afac_matrix,
                    const TScalarVector *  D2,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : _has_scaling( true )
            , _eval_type(   aeval_type )
            , _storage_type( astorage_type )
            , _acc( acc_exact )
    {
        set_matrix( std::forward< T_mat >( afac_matrix ) );
        
        if ( _fac_matrix == nullptr )
            HERROR( ERR_ARG, "(TFacInvMatrix)", "matrix is nullptr" );
    
        if (( D1 == nullptr ) || ( D2 == nullptr ))
            HERROR( ERR_ARG, "(TFacInvMatrix)", "scaling factors are nullptr" );

        _D_left.copy_from( D1 );
        _D_right.copy_from( D2 );
        
        init();
    }

    //! dtor
    virtual ~TFacInvMatrix () {}

    ////////////////////////////////////////////////////////
    //
    // structure of matrix
    //

    //! access factorised matrix
    const TMatrix *        matrix        () const { return _fac_matrix; }
    
    //! return true if diagonal scaling is used
    bool                   has_scaling   () const { return _has_scaling; }

    //! return left scaling matrix
    const TScalarVector &  scaling_left  () const { return _D_left; }

    //! return right scaling matrix
    const TScalarVector &  scaling_right () const { return _D_right; }
    
    //! return evaluation type (pointwise/blockwise)
    eval_type_t            eval_type     () const { return _eval_type; }

    //! return evaluation type (pointwise/blockwise)
    storage_type_t         storage_type  () const { return _storage_type; }

    //! return internal accuracy for matrix arithmetic
    const TTruncAcc        accuracy      () const { return _acc; }
    
    //! set internal accuracy for matrix arithmetic
    void                   set_accuracy  ( const TTruncAcc &  aacc ) { _acc = aacc; }
    
    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    bool  is_complex      () const
    {
        return ( _fac_matrix->is_complex() || _D_left.is_complex() || _D_right.is_complex() );
    }
    
    //! return true, of operator is self adjoint
    bool  is_self_adjoint () const
    {
        return _self_adjoint;
    }
    
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
    virtual size_t  domain_dim () const { return _fac_matrix->range_dim(); }
    
    //! return dimension of range
    virtual size_t  range_dim  () const { return _fac_matrix->domain_dim(); }
    
    //! return vector in domain space
    virtual auto domain_vector  () const -> std::unique_ptr< TVector > { return _fac_matrix->row_vector(); }

    //! return vector in range space
    virtual auto range_vector   () const -> std::unique_ptr< TVector > { return _fac_matrix->col_vector(); }

    ///////////////////////////////////////////////////////////
    //
    // internal solve method for factors
    //

    virtual void solve ( TVector *      x,
                         const matop_t  op ) const = 0;

    virtual void solve ( TMatrix *      X,
                         const matop_t  op ) const = 0;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TFacInvMatrix, TLinearOperator );

private:
    ///////////////////////////////////////////////////////////
    //
    // internal initialisation
    //

    void  init ();
    
    ///////////////////////////////////////////////////////////
    //
    // for perfect forwarding
    //

    // WITHOUT ownership transfer
    void  set_matrix ( const TMatrix *  A )
    {
        _fac_matrix = A;
    }

    // WITH ownership transfer
    void  set_matrix ( std::unique_ptr< TMatrix > &&  A )
    {
        _mat_ptr    = std::move( A );
        _fac_matrix = _mat_ptr.get();
    }

};

//////////////////////////////////////////////////////////////////////////////////
//
// TLUInvMatrix
//
//////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup  Matrix_Module
//! \class    TLUInvMatrix
//! \brief    Represents the inverse of a LU factored matrix
//!
//!           Evaluates \f$ (LU)^{-1} \f$.
//!
class TLUInvMatrix : public TFacInvMatrix
{
public:
    ////////////////////////////////////////////////////////
    //
    // constructors and destructor
    //

    //! construct inverse operator with LU factorised matrix \a afac_matrix
    //! (no diagonal scaling applied)
    template < typename T_mat >
    TLUInvMatrix ( T_mat &&               afac_matrix,
                   const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                   const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( std::forward< T_mat >( afac_matrix ), aeval_type, astorage_type )
    {}

    //! construct inverse operator with LU factorised matrix \a afac_matrix
    //! with accuracy for matrix arithmetic (without diagonal scaling)
    template < typename T_mat >
    TLUInvMatrix ( T_mat &&               afac_matrix,
                   const TTruncAcc &      aacc,
                   const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                   const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( std::forward< T_mat >( afac_matrix ), aacc, aeval_type, astorage_type )
    {}

    //! construct inverse operator with LU factorised matrix \a afac_matrix
    //! and additional row and column scaling factors \a D1 and \a D2
    template < typename T_mat >
    TLUInvMatrix ( const TScalarVector *  D1,
                   T_mat &&               afac_matrix,
                   const TScalarVector *  D2,
                   const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                   const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( D1, std::forward< T_mat >( afac_matrix ), D2, aeval_type, astorage_type )
    {}

    //! dtor
    virtual ~TLUInvMatrix ()
    {}
    
    //!
    //! solve op(A) x = y with given \a y (supplied in \a x)
    //!
    virtual void solve ( TVector *      x,
                         const matop_t  op ) const;

    virtual void solve ( TMatrix *      X,
                         const matop_t  op ) const;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TLUInvMatrix, TFacInvMatrix );
};

//////////////////////////////////////////////////////////////////////////////////
//
// TLDUInvMatrix
//
//////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup  Matrix_Module
//! \class    TLDUInvMatrix
//! \brief    Represents the inverse of a LDU factored matrix
//!
//!           Evaluates \f$ (LDU)^{-1} \f$.
//!
class TLDUInvMatrix : public TFacInvMatrix
{
public:
    ////////////////////////////////////////////////////////
    //
    // constructors and destructor
    //

    //! construct inverse operator with LDU factorised matrix \a afac_matrix
    //! (no diagonal scaling applied)
    template < typename T_mat >
    TLDUInvMatrix ( T_mat &&               afac_matrix,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( std::forward< T_mat >( afac_matrix ), aeval_type, astorage_type )
    {}

    //! construct inverse operator with LDU factorised matrix \a afac_matrix
    //! with accuracy for matrix arithmetic (without diagonal scaling)
    template < typename T_mat >
    TLDUInvMatrix ( T_mat &&               afac_matrix,
                    const TTruncAcc &      aacc,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( std::forward< T_mat >( afac_matrix ), aacc, aeval_type, astorage_type )
    {}

    //! construct inverse operator with LDU factorised matrix \a afac_matrix
    //! and additional row and column scaling factors \a D1 and \a D2
    template < typename T_mat >
    TLDUInvMatrix ( const TScalarVector *  D1,
                    T_mat &&               afac_matrix,
                    const TScalarVector *  D2,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( D1, std::forward< T_mat >( afac_matrix ), D2, aeval_type, astorage_type )
    {}

    //! dtor
    virtual ~TLDUInvMatrix ()
    {}
    
    //!
    //! solve op(A) x = y with given \a y (supplied in \a x)
    //!
    virtual void solve ( TVector *      x,
                         const matop_t  op ) const;

    virtual void solve ( TMatrix *      X,
                         const matop_t  op ) const;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TLDUInvMatrix, TFacInvMatrix );
};

//////////////////////////////////////////////////////////////////////////////////
//
// TLDLInvMatrix
//
//////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup  Matrix_Module
//! \class    TLDLInvMatrix
//! \brief    Represents the inverse of a LDL factored matrix
//!
//!           Evaluates \f$ (LDL^H)^{-1} \f$ or \f$ (LDL^T)^{-1} \f$.
//!
class TLDLInvMatrix : public TFacInvMatrix
{
private:
    // matrix format of original matrix before factorisation
    const matform_t  _format_A;
    
public:
    ////////////////////////////////////////////////////////
    //
    // constructors and destructor
    //

    //! construct inverse operator with LDL factorised matrix \a afac_matrix
    //! which has \a format (symmetric, etc.)
    //! (no diagonal scaling applied)
    template < typename T_mat >
    TLDLInvMatrix ( T_mat &&               afac_matrix,
                    const matform_t        aformat,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( std::forward< T_mat >( afac_matrix ), aeval_type, astorage_type )
            , _format_A( aformat )
    {}

    //! construct inverse operator with LDL factorised matrix \a afac_matrix
    //! with accuracy for matrix arithmetic (without diagonal scaling)
    template < typename T_mat >
    TLDLInvMatrix ( T_mat &&               afac_matrix,
                    const matform_t        aformat,
                    const TTruncAcc &      aacc,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( std::forward< T_mat >( afac_matrix ), aacc, aeval_type, astorage_type )
            , _format_A( aformat )
    {}

    //! construct inverse operator with LDL factorised matrix \a afac_matrix
    //! which has \a format (symmetric, etc.) and additional row and column
    //! scaling factors \a D1 and \a D2
    template < typename T_mat >
    TLDLInvMatrix ( const TScalarVector *  D1,
                    T_mat &&               afac_matrix,
                    const TScalarVector *  D2,
                    const matform_t        aformat,
                    const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                    const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( D1, std::forward< T_mat >( afac_matrix ), D2, aeval_type, astorage_type )
            , _format_A( aformat )
    {}

    //! dtor
    virtual ~TLDLInvMatrix ()
    {}
    
    //!
    //! solve op(A) x = y with given \a y (supplied in \a x)
    //!
    virtual void solve ( TVector *      x,
                         const matop_t  op ) const;

    virtual void solve ( TMatrix *      X,
                         const matop_t  op ) const;

    //
    // return matrix format before factorisation
    //
    matform_t  matrix_format () const { return _format_A; }

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TLDLInvMatrix, TFacInvMatrix );
};

//////////////////////////////////////////////////////////////////////////////////
//
// TLLInvMatrix
//
//////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup  Matrix_Module
//! \class    TLLInvMatrix
//! \brief    Represents the inverse of a Cholesky factored matrix
//!
//!           Evaluates \f$ (LL^H)^{-1} \f$ or \f$ (LL^T)^{-1} \f$.
//!
class TLLInvMatrix : public TFacInvMatrix
{
private:
    //! @cond
    
    // matrix format of original matrix before factorisation
    const matform_t  _format_A;
    
    //! @endcond
    
public:
    ////////////////////////////////////////////////////////
    //
    // constructors and destructor
    //

    //! construct inverse operator with Cholesky factorised matrix \a afac_matrix
    //! which has \a format (symmetric, etc.)
    //! (no diagonal scaling applied)
    template < typename T_mat >
    TLLInvMatrix ( T_mat &&               afac_matrix,
                   const matform_t        aformat,
                   const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( std::forward< T_mat >( afac_matrix ), point_wise, astorage_type )
            , _format_A( aformat )
    {}

    //! construct inverse operator with LL factorised matrix \a afac_matrix
    //! with accuracy for matrix arithmetic (without diagonal scaling)
    template < typename T_mat >
    TLLInvMatrix ( T_mat &&               afac_matrix,
                   const matform_t        aformat,
                   const TTruncAcc &      aacc,
                   const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                   const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( std::forward< T_mat >( afac_matrix ), aacc, aeval_type, astorage_type )
            , _format_A( aformat )
    {}

    //! construct inverse operator with Cholesky factorised matrix \a afac_matrix
    //! which has \a format (symmetric, etc.) and additional row and column
    //! scaling factors \a D1 and \a D2
    template < typename T_mat >
    TLLInvMatrix ( const TScalarVector *  D1,
                   T_mat &&               afac_matrix,
                   const TScalarVector *  D2,
                   const matform_t        format,
                   const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacInvMatrix( D1, std::forward< T_mat >( afac_matrix ), D2, point_wise, astorage_type )
            , _format_A( format )
    {}

    //! dtor
    virtual ~TLLInvMatrix ()
    {}
    
    //!
    //! solve op(A) x = y with given \a y (supplied in \a x)
    //!
    virtual void solve ( TVector *      x,
                         const matop_t  op ) const;

    virtual void solve ( TMatrix *      X,
                         const matop_t  op ) const;

    //
    // return matrix format before factorisation
    //
    matform_t  matrix_format () const { return _format_A; }

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TLLInvMatrix, TFacInvMatrix );
};

}// namespace

#endif  // __HLIB_TFACINVMATRIX_HH
