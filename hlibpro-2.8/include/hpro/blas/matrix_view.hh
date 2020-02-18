#ifndef __HLIB_BLAS_MATRIX_VIEW_HH
#define __HLIB_BLAS_MATRIX_VIEW_HH
//
// Project     : HLib
// File        : Matrix.hh
// Description : provides transposed and adjoint matrix views
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/base/traits.hh"

#include "hpro/blas/types.hh"
#include "hpro/blas/MatrixBase.hh"

namespace HLIB
{

namespace BLAS
{

//!
//! \ingroup  BLAS_Module
//! \brief    return conjugate of given matrix operation
//! \param    op  matrix op. to be conjugated
//!
inline
matop_t
conjugate ( const matop_t   op )
{
    switch ( op )
    {
        case MATOP_NORM  : HERROR( ERR_CONSISTENCY, "conjugate", "conjugate not supported" );
        case MATOP_TRANS : return MATOP_ADJ ;
        case MATOP_ADJ   : return MATOP_TRANS;
        default          : HERROR( ERR_CONSISTENCY, "conjugate", "unknown matrix operation" );
    }// switch
}

//!
//! \ingroup  BLAS_Module
//! \brief    return transposed of given matrix operation
//! \param    op  matrix op. to be transposed
//!
inline
matop_t
transposed ( const matop_t   op )
{
    switch ( op )
    {
        case MATOP_NORM  : return MATOP_TRANS ;
        case MATOP_TRANS : return MATOP_NORM ;
        case MATOP_ADJ   : HERROR( ERR_CONSISTENCY, "transposed", "conjugate not supported" );
        default          : HERROR( ERR_CONSISTENCY, "transposed", "unknown matrix operation" );
    }// switch
}

//!
//! \ingroup  BLAS_Module
//! \brief    return adjoint of given matrix operation
//! \param    op  matrix op. to be adjoint
//!
inline
matop_t
adjoint ( const matop_t   op )
{
    switch ( op )
    {
        case MATOP_NORM  : return MATOP_ADJ ;
        case MATOP_TRANS : HERROR( ERR_CONSISTENCY, "adjoint", "conjugate not supported" );
        case MATOP_ADJ   : return MATOP_NORM;
        default          : HERROR( ERR_CONSISTENCY, "adjoint", "unknown matrix operation" );
    }// switch
}

//!
//! \{
//! \name Matrix Modifiers
//! Classes for matrix modifiers, e.g. transposed and adjoint view.
//!

//!
//! \ingroup  BLAS_Module
//! \class    TransposeView
//! \brief    Provide transposed view of a matrix.
//!
template < typename T_matrix >
class TransposeView : public MatrixBase< TransposeView< T_matrix > >
{
public:
    using  value_t = typename T_matrix::value_t;

private:
    const T_matrix &  _mat;

public:
    TransposeView ( const T_matrix &  M ) noexcept
            : _mat( M )
    {}

    value_t *        data         () const noexcept { return _mat.data(); }

    size_t           nrows        () const noexcept { return _mat.ncols(); }
    size_t           ncols        () const noexcept { return _mat.nrows(); }

    size_t           blas_nrows   () const noexcept { return _mat.nrows(); }
    size_t           blas_ncols   () const noexcept { return _mat.ncols(); }

    size_t           row_stride   () const noexcept { return _mat.row_stride(); }
    size_t           col_stride   () const noexcept { return _mat.col_stride(); }

    blasview_t       blas_view    () const
    {
        switch ( _mat.blas_view() )
        {
        case BLAS_NORMAL     :
            return BLAS_TRANSPOSED;
            
        case BLAS_TRANSPOSED :
            return BLAS_NORMAL;
            
        case BLAS_ADJOINT    :
            // just conjugate without transpose is not possible with standard BLAS/LAPACK
            if ( is_complex_type< value_t >::value )
                HERROR( ERR_NOT_IMPL, "", "" );
            
            return BLAS_NORMAL;
        }// switch

        return BLAS_TRANSPOSED;
    }
    
    value_t    operator ()  ( const idx_t i, const idx_t j ) const noexcept { return _mat( j, i ); }
};

//
// trait for providing matrix properties
//
template < typename T >
struct matrix_trait< TransposeView< T > >
{
    using  value_t = typename T::value_t;
};

//
// signals, that T is of matrix type
//
template < typename T >
struct is_matrix< TransposeView< T > >
{
    static const bool  value = true;
};

//!
//! \ingroup  BLAS_Module
//! \brief  return transposed view object for matrix
//! \param  M  matrix to transpose
//!
template < typename T >
TransposeView< T >
transposed ( const T &  M )
{
    return TransposeView< T >( M );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup  BLAS_Module
//! \class    AdjoinView
//! \brief    Provide adjoint view, e.g. conjugate transposed of a given matrix.
//!
template < typename T_matrix >
class AdjoinView : public MatrixBase< AdjoinView< T_matrix > >
{
public:
    using  value_t = typename T_matrix::value_t;

private:
    const T_matrix &  _mat;

public:
    AdjoinView ( const T_matrix &  M ) noexcept
            : _mat( M )
    {}

    value_t *        data         () const noexcept { return _mat.data(); }

    size_t           nrows        () const noexcept { return _mat.ncols(); }
    size_t           ncols        () const noexcept { return _mat.nrows(); }

    size_t           blas_nrows   () const noexcept { return _mat.nrows(); }
    size_t           blas_ncols   () const noexcept { return _mat.ncols(); }

    size_t           row_stride   () const noexcept { return _mat.row_stride(); }
    size_t           col_stride   () const noexcept { return _mat.col_stride(); }

    blasview_t       blas_view    () const
    {
        switch ( _mat.blas_view() )
        {
        case BLAS_NORMAL     :
            return BLAS_ADJOINT;
            
        case BLAS_TRANSPOSED :
            // just conjugate without transpose is not possible with standard BLAS/LAPACK
            if ( is_complex_type< value_t >::value )
                HERROR( ERR_NOT_IMPL, "", "" );
            
            return BLAS_NORMAL;
            
        case BLAS_ADJOINT    :
            return BLAS_NORMAL;
        }// switch

        return BLAS_ADJOINT;
    }
    
    value_t    operator ()  ( const idx_t i, const idx_t j ) const noexcept { return HLIB::conj( _mat( j, i ) ); }
};

//
// trait for providing matrix properties
//
template < typename T >
struct matrix_trait< AdjoinView< T > >
{
    using  value_t = typename T::value_t;
};

//
// signals, that T is of matrix type
//
template < typename T >
struct is_matrix< AdjoinView< T > >
{
    static const bool  value = true;
};

//!
//! \ingroup  BLAS_Module
//! \brief  return adjoint view object for matrix
//! \param  M  matrix to conjugate transpose
//!
template < typename T >
AdjoinView< T >
adjoint ( const T &  M )
{
    return AdjoinView< T >( M );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

//!
//! \ingroup  BLAS_Module
//! \class    MatrixView
//! \brief    Provide generic view to a matrix, e.g. transposed or adjoint.
//!
template < typename T_matrix >
class MatrixView : public MatrixBase< MatrixView< T_matrix > >
{
public:
    using  value_t = typename T_matrix::value_t;

private:
    const matop_t     _op;
    const T_matrix &  _mat;

public:
    MatrixView ( const matop_t     op,
                 const T_matrix &  M ) noexcept
            : _op(op), _mat( M )
    {}

    value_t *        data         () const noexcept { return _mat.data(); }

    size_t           nrows        () const noexcept
    {
        switch ( _op )
        {
            case MATOP_NORM  : return _mat.nrows();
            case MATOP_TRANS :
            case MATOP_ADJ   : return _mat.ncols();
        }// switch

        return _mat.nrows();
    }

    size_t           ncols        () const noexcept
    {
        switch ( _op )
        {
            case MATOP_NORM  : return _mat.ncols();
            case MATOP_TRANS :
            case MATOP_ADJ   : return _mat.nrows();
        }// switch

        return _mat.ncols();
    }

    size_t           blas_nrows   () const noexcept { return _mat.nrows(); }
    size_t           blas_ncols   () const noexcept { return _mat.ncols(); }

    size_t           row_stride   () const noexcept { return _mat.row_stride(); }
    size_t           col_stride   () const noexcept { return _mat.col_stride(); }

    blasview_t       blas_view    () const
    {
        if ( _op == MATOP_NORM )
        {
            return _mat.blas_view();
        }// if
        else if ( _op == MATOP_TRANS )
        {
            switch ( _mat.blas_view() )
            {
                case BLAS_NORMAL     :
                    return BLAS_TRANSPOSED;
            
                case BLAS_TRANSPOSED :
                    return BLAS_NORMAL;
            
                case BLAS_ADJOINT    :
                    // just conjugate without transpose is not possible with standard BLAS/LAPACK
                    if ( is_complex_type< value_t >::value )
                        HERROR( ERR_NOT_IMPL, "", "" );
            
                    return BLAS_NORMAL;
            }// switch
        }// if
        else if ( _op == MATOP_ADJ )
        {
            switch ( _mat.blas_view() )
            {
                case BLAS_NORMAL     :
                    return BLAS_ADJOINT;
            
                case BLAS_TRANSPOSED :
                    // just conjugate without transpose is not possible with standard BLAS/LAPACK
                    if ( is_complex_type< value_t >::value )
                        HERROR( ERR_NOT_IMPL, "", "" );
            
                    return BLAS_NORMAL;
            
                case BLAS_ADJOINT    :
                    return BLAS_NORMAL;
            }// switch
        }// else
        
        return BLAS_NORMAL;
    }
    
    value_t    operator ()  ( const idx_t i, const idx_t j ) const noexcept
    {
        switch ( _op )
        {
            case MATOP_NORM  : return _mat( i, j );
            case MATOP_TRANS : return _mat( j, i );
            case MATOP_ADJ   : return HLIB::conj( _mat( j, i ) );
        }// switch

        return _mat( i, j );
    }
};

//
// trait for providing matrix properties
//
template < typename T >
struct matrix_trait< MatrixView< T > >
{
    using  value_t = typename T::value_t;
};

//
// signals, that T is of matrix type
//
template < typename T >
struct is_matrix< MatrixView< T > >
{
    static const bool  value = true;
};

//!
//! \ingroup  BLAS_Module
//! \brief  convert matop_t into view object
//! \param  op  matop_t value
//! \param  M   matrix to create view for
//!
template < typename T >
MatrixView< T >
mat_view ( const matop_t  op,
           const T &      M )
{
    return MatrixView< T >( op, M );
}

//! \}

}// namespace BLAS

}// namespace HLIB

#endif  // __HLIB_BLAS_MATRIX_VIEW_HH
