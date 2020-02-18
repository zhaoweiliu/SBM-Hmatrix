#ifndef __HLIB_TFACMATRIX_HH
#define __HLIB_TFACMATRIX_HH
//
// Project     : HLib
// File        : TFacMatrix.hh
// Description : classes for representing factorised matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <utility>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/algebra/solve_types.hh"

namespace HLIB
{

// local RTTI types
DECLARE_TYPE( TFacMatrix );
DECLARE_TYPE( TLUMatrix );
DECLARE_TYPE( TLDUMatrix );
DECLARE_TYPE( TLDLMatrix );
DECLARE_TYPE( TLLMatrix );

//////////////////////////////////////////////////////////////////////////////////
//
// TFacMatrix
//
//////////////////////////////////////////////////////////////////////////////////

// local matrix type
// DECLARE_TYPE( TFacMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TFacMatrix
//! \brief    Baseclass for representing factorised matrices.
//!
//!           TFacMatrix and the derived classes represent operators based on
//!           factorised matrices providing the ability to evaluate the corresponding
//!           matrix (in factorised form), e.g. for \f$ A = LU \f$ provide f$ (LU) x = y \f$.
//!
//!           Optionally, diagonal scalings from left and right may have been
//!           applied before the factorisation and now taken into account in the
//!           evaluation.
//!
class TFacMatrix : public TLinearOperator
{
private:

    //! @cond
    
    //! factorised matrix
    const TMatrix *           _fac_matrix;

    //! if true, matrix has diagonal scaling applied
    bool                      _has_scaling;

    //! diagonal scaling factor from left
    TScalarVector             _D_left;

    //! diagonal scaling factor from right
    TScalarVector             _D_right;

    //! determines pointwise of blockwise evaluation
    eval_type_t               _eval_type;

    //! determines type of storage for diagonal blocks
    storage_type_t            _storage_type;

    //! signal if operator is self adjoint
    bool                      _self_adjoint;

    //! auto pointer in case of ownership transfer
    std::unique_ptr< const TMatrix >  _mat_ptr;
    
    //! @endcond
    
public:
    //!
    //! constructor without diagonal scaling
    //!
    template < typename T_mat >
    TFacMatrix ( T_mat &&               afac_matrix,
                 const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                 const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : _has_scaling( false )
            , _eval_type(    aeval_type )
            , _storage_type( astorage_type )
    {
        set_matrix( std::forward< T_mat >( afac_matrix ) );
        
        if ( _fac_matrix == nullptr )
            HERROR( ERR_ARG, "(TFacMatrix)", "LU factor is nullptr" );
        
        init();
    }

    //!
    //! constructor with diagonal scaling applied
    //!
    template < typename T_mat >
    TFacMatrix ( const TScalarVector *  D1,
                 T_mat &&               afac_matrix,
                 const TScalarVector *  D2,
                 const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                 const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : _has_scaling( true )
            , _eval_type(   aeval_type )
            , _storage_type( astorage_type )
    {
        set_matrix( std::forward< T_mat >( afac_matrix ) );
        
        if ( _fac_matrix == nullptr )
            HERROR( ERR_ARG, "(TFacMatrix)", "matrix is nullptr" );
    
        if (( D1 == nullptr ) || ( D2 == nullptr ))
            HERROR( ERR_ARG, "(TFacMatrix)", "scaling factors are nullptr" );

        _D_left.copy_from( D1 );
        _D_right.copy_from( D2 );
        
        init();
    }

    //! dtor
    virtual ~TFacMatrix () {}

    ////////////////////////////////////////////////////////
    //
    // structure of matrix
    //

    //! access internal factorised matrix
    const TMatrix *        matrix        () const { return _fac_matrix; }

    //! return info about presence of scaling matrices
    bool                   has_scaling   () const { return _has_scaling; }

    //! return left scaling matrix
    const TScalarVector &  scaling_left  () const { return _D_left; }

    //! return right scaling matrix
    const TScalarVector &  scaling_right () const { return _D_right; }
    
    //! return evaluation type (pointwise/blockwise)
    eval_type_t            eval_type     () const { return _eval_type; }

    //! return evaluation type (pointwise/blockwise)
    storage_type_t         storage_type  () const { return _storage_type; }

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
    virtual size_t  domain_dim () const { return _fac_matrix->domain_dim(); }
    
    //! return dimension of range
    virtual size_t  range_dim  () const { return _fac_matrix->range_dim(); }
    
    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector > { return _fac_matrix->col_vector(); }

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector > { return _fac_matrix->row_vector(); }

    ///////////////////////////////////////////////////////////
    //!
    //! internal evaluation method for factors
    //!

    virtual void eval ( TVector *      x,
                        const matop_t  op ) const = 0;

    virtual void eval ( TMatrix *      X,
                        const matop_t  op ) const = 0;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TFacMatrix, TLinearOperator );

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
// TLUMatrix
//
//////////////////////////////////////////////////////////////////////////////////

// local matrix type
// DECLARE_TYPE( TLUMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TLUMatrix
//! \brief    Represents operator for LU factored matrices
//!
//!           Evaluates \f$ L \cdot U \f$.
//!
class TLUMatrix : public TFacMatrix
{
public:
    //!
    //! constructor without diagonal scaling
    //!
    template < typename T_mat >
    TLUMatrix ( T_mat &&               afac_matrix,
                const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacMatrix( std::forward< T_mat >( afac_matrix ), aeval_type, astorage_type )
    {}

    //!
    //! constructor with diagonal scaling
    //!
    template < typename T_mat >
    TLUMatrix ( const TScalarVector *  D1,
                T_mat &&               afac_matrix,
                const TScalarVector *  D2,
                const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacMatrix( D1, std::forward< T_mat >( afac_matrix ), D2, aeval_type, astorage_type )
    {}

    //! dtor
    virtual ~TLUMatrix () {}
    
    //!
    //! internal evaluation method for factors
    //!

    virtual void eval ( TVector *      x,
                        const matop_t  op ) const;

    virtual void eval ( TMatrix *      X,
                        const matop_t  op ) const;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TLUMatrix, TFacMatrix );
};

//////////////////////////////////////////////////////////////////////////////////
//
// TLDUMatrix
//
//////////////////////////////////////////////////////////////////////////////////

// local matrix type
// DECLARE_TYPE( TLDUMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TLDUMatrix
//! \brief    Represents operator for LDU factored matrices
//!
//!           Evaluates \f$ L \cdot D \cdot U \f$.
//!
class TLDUMatrix : public TFacMatrix
{
public:
    //!
    //! constructor without diagonal scaling
    //!
    template < typename T_mat >
    TLDUMatrix ( T_mat &&               afac_matrix,
                 const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                 const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacMatrix( std::forward< T_mat >( afac_matrix ), aeval_type, astorage_type )
    {}

    //!
    //! constructor with diagonal scaling
    //!
    template < typename T_mat >
    TLDUMatrix ( const TScalarVector *  D1,
                 T_mat &&               afac_matrix,
                 const TScalarVector *  D2,
                 const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                 const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacMatrix( D1, std::forward< T_mat >( afac_matrix ), D2, aeval_type, astorage_type )
    {}

    //! dtor
    virtual ~TLDUMatrix () {}
    
    //!
    //! internal evaluation method for factors
    //!

    virtual void eval ( TVector *      x,
                        const matop_t  op ) const;

    virtual void eval ( TMatrix *      X,
                        const matop_t  op ) const;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TLDUMatrix, TFacMatrix );
};

//////////////////////////////////////////////////////////////////////////////////
//
// TLDLMatrix
//
//////////////////////////////////////////////////////////////////////////////////

// local matrix type
// DECLARE_TYPE( TLDLMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TLDLMatrix
//! \brief    Represents operator for LDL factored matrices
//!
//!           Evaluates \f$ L \cdot D \cdot L^H \f$ or \f$ L \cdot D \cdot L^T \f$.
//!
class TLDLMatrix : public TFacMatrix
{
private:
    // matrix format of original matrix before factorisation
    const matform_t  _format_A;
    
public:
    //!
    //! constructor without diagonal scaling
    //!
    template < typename T_mat >
    TLDLMatrix ( T_mat &&               afac_matrix,
                 const matform_t        format,
                 const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                 const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacMatrix( std::forward< T_mat >( afac_matrix ), aeval_type, astorage_type )
            , _format_A( format )
    {}

    //!
    //! constructor with diagonal scaling
    //!
    template < typename T_mat >
    TLDLMatrix ( const TScalarVector *  D1,
                 T_mat &&               afac_matrix,
                 const TScalarVector *  D2,
                 const matform_t        format,
                 const eval_type_t      aeval_type    = CFG::Arith::eval_type,
                 const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacMatrix( D1, std::forward< T_mat >( afac_matrix ), D2, aeval_type, astorage_type )
            , _format_A( format )
    {}

    //! dtor
    virtual ~TLDLMatrix () {}
    
    //!
    //! internal evaluation method for factors
    //!

    virtual void eval ( TVector *      x,
                        const matop_t  op ) const;

    virtual void eval ( TMatrix *      X,
                        const matop_t  op ) const;

    //
    // return matrix format before factorisation
    //
    matform_t  matrix_format () const { return _format_A; }

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //
    
    HLIB_RTTI_DERIVED( TLDLMatrix, TFacMatrix );
};

//////////////////////////////////////////////////////////////////////////////////
//
// TLLMatrix
//
//////////////////////////////////////////////////////////////////////////////////

// local matrix type
// DECLARE_TYPE( TLLMatrix );

//!
//! \ingroup  Matrix_Module
//! \class    TLLMatrix
//! \brief    Represents operator for Cholesky factored matrices
//!
//!           Evaluates \f$ L \cdot L^H \f$ or \f$ L \cdot L^T \f$.
//!
class TLLMatrix : public TFacMatrix
{
private:
    //! @cond

    // matrix format of original matrix before factorisation
    const matform_t  _format_A;

    //! @endcond
    
public:
    //!
    //! constructor without diagonal scaling
    //!
    template < typename T_mat >
    TLLMatrix ( T_mat &&               afac_matrix,
                const matform_t        format,
                const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacMatrix( std::forward< T_mat >( afac_matrix ), point_wise, astorage_type )
            , _format_A( format )
    {}

    //!
    //! constructor with diagonal scaling
    //!
    template < typename T_mat >
    TLLMatrix ( const TScalarVector *  D1,
                T_mat &&               afac_matrix,
                const TScalarVector *  D2,
                const matform_t        format,
                const storage_type_t   astorage_type = CFG::Arith::storage_type )
            : TFacMatrix( D1, std::forward< T_mat >( afac_matrix ), D2, point_wise, astorage_type )
            , _format_A( format )
    {}

    //! dtor
    virtual ~TLLMatrix () {}
    
    //!
    //! internal evaluation method for factors
    //!

    virtual void eval ( TVector *      x,
                        const matop_t  op ) const;

    virtual void eval ( TMatrix *      X,
                        const matop_t  op ) const;

    //
    // return matrix format before factorisation
    //
    matform_t  matrix_format () const { return _format_A; }

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TLLMatrix, TFacMatrix );
};

}// namespace HLIB

#endif  // __HLIB_TFACMATRIX_HH
