#ifndef __HLIB_BLAS_VECTOR_HH
#define __HLIB_BLAS_VECTOR_HH
//
// Project     : HLib
// File        : Vector.hh
// Description : implements dense vector class for BLAS operations
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <ostream>

#include "hpro/base/error.hh"

#include "hpro/blas/VectorBase.hh"
#include "hpro/blas/MatrixBase.hh"
#include "hpro/blas/MemBlock.hh"

namespace HLIB
{

namespace BLAS
{

template < typename T > class Matrix;

//!
//! \ingroup  BLAS_Module
//! \class    Vector
//! \brief    Standard vector in basic linear algebra, i.e. BLAS/LAPACK.
//!
//! \detail   Standard vector class in linear algebra routines. It holds a
//!           reference to a memory block, which is controlled by this vector,
//!           e.g. if the vector is gone or resized, the memory block is also gone.
//!
template < typename T_value >
class Vector : public VectorBase< Vector< T_value > >, public MemBlock< T_value >
{
public:
    //! internal value type
    using   value_t = T_value;

    //! super class type
    using   super_t = MemBlock< value_t >;
    
private:
    //! length of vector
    size_t   _length;

    //! stride of data in memory block
    size_t   _stride;
    
public:
    //
    // constructor and destructor
    //

    //! create zero sized vector
    Vector () noexcept
            : MemBlock< value_t >()
            , _length( 0 )
            , _stride( 0 )
    {}

    //! create vector of size \a n
    Vector ( const size_t  n )
            : MemBlock<value_t>( n )
            , _length( n )
            , _stride( 1 )
    {}

    //! copy ctor
    Vector ( const Vector &       v,
             const copy_policy_t  p = copy_reference )
            : MemBlock< value_t >()
            , _length( 0 )
            , _stride( 0 )
    {
        switch ( p )
        {
            case copy_reference :
                (*this) = v;
                break;

            case copy_value :
                super_t::alloc_wo_value( v.length() );
                _length = v.length();
                _stride = 1;
                
                for ( idx_t i = 0; i < idx_t( _length ); i++ )
                    (*this)(i) = v(i);
                
                break;
        }// switch
    }
    
    //! move ctor
    Vector ( Vector &&  v ) noexcept
            : MemBlock< value_t >( std::move( v ) )
            , _length( 0 )
            , _stride( 0 )
    {
        v._data = nullptr;
        
        std::swap( _length, v._length );
        std::swap( _stride, v._stride );
    }

    // create copy of sub vector of \a v defined by range \a r
    Vector ( const Vector< value_t > & v,
             const Range &             ar,
             const copy_policy_t       p = copy_reference )
            : MemBlock< value_t >()
            , _length( 0 )
            , _stride( 0 )
    {
        const Range  r( ar == Range::all ? Range( 0, idx_t(v.length())-1 ) : ar );
        
        HASSERT( r.size() * r.stride() <= v.length(),
                 ERR_SIZE, "(Vector) ctor", "range size × range stride > vector length" );

        _length = r.size();

        switch ( p )
        {
            case copy_reference :
                _stride = r.stride() * v.stride();
                super_t::init( v.data() + r.first() * v.stride() );
                break;

            case copy_value :
                super_t::alloc_wo_value( _length );
                _stride = 1;
                for ( idx_t  i = 0; i < idx_t(_length); ++i )
                    (*this)(i) = v( r.first() + i * idx_t(r.stride()) );
                break;
        }// switch
    }

    //! create copy of column of matrix \a M defined by \a r and \a col
    Vector ( const Matrix< value_t > &  M,
             const Range &              ar,
             const idx_t                col,
             const copy_policy_t        p = copy_reference )
            : MemBlock< value_t >()
            , _length( 0 )
            , _stride( 0 )
    {
        const Range  r( ar == Range::all ? M.row_range() : ar );

        if ( r.size() * r.stride() > M.nrows() )
            HERROR( ERR_SIZE, "(Vector) ctor", "range size × range stride > matrix rows" );

        _length = r.size();

        switch ( p )
        {
            case copy_reference :
                _stride = r.stride() * M.row_stride();
                super_t::init( M.data() + r.first() * M.row_stride() + col * M.col_stride() );
                break;
                
            case copy_value :
                super_t::alloc_wo_value( _length );
                _stride = 1;
                for ( idx_t  i = 0; i < idx_t(_length); i++ )
                    (*this)(i) = M( r.first() + i * idx_t(r.stride()), col );
                break;
        }// switch
    }
    
    //! create copy/reference to row of matrix \a M defined by \a r and \a row
    Vector ( const Matrix< value_t > &  M,
             const idx_t                row,
             const Range &              ar,
             const copy_policy_t        p = copy_reference )
            : MemBlock< value_t >()
            , _length( 0 )
            , _stride( 0 )
    {
        const Range  r( ar == Range::all ? M.col_range() : ar );
        
        if ( r.size() * r.stride() > M.ncols() )
            HERROR( ERR_SIZE, "(Vector) ctor", "range size × range stride > matrix columns" );

        _length = r.size();

        switch ( p )
        {
            case copy_reference :
                _stride = r.stride() * M.col_stride();
                super_t::init( M.data() + r.first() * M.col_stride() + row * M.row_stride() );
                break;
                
            case copy_value :
                super_t::alloc_wo_value( _length );
                _stride = 1;
        
                for ( idx_t  i = 0; i < idx_t(_length); i++ )
                    (*this)(i) = M( row, r.first() + i * idx_t(r.stride()) );
                break;
        }// switch
    }

    //! copy operator (always copy reference! for real data copy, use copy function below!)
    Vector & operator = ( const Vector & v )
    {
        super_t::init( v.data(), false );
        _length = v.length();
        _stride = v.stride();
        
        return *this;
    }

    //! move operator (move ownership of data)
    Vector & operator = ( Vector && v ) noexcept
    {
        if ( this != & v ) // prohibit self-moving
        {
            super_t::init( v, v._is_owner );
            _length = v.length();
            _stride = v.stride();

            // reset data of v
            v._data   = nullptr;
            v._length = 0;
            v._stride = 0;
        }// if
        
        return *this;
    }
    
    //
    // data access
    //

    //! return length of vector
    size_t    length      () const noexcept { return _length; }

    //! return stride of index set
    size_t    stride      () const noexcept { return _stride; }
    
    //! return coefficient at position \a i
    value_t   operator () ( const idx_t  i ) const noexcept
    {
        HASSERT( i < idx_t(_length), ERR_ARR_BOUND, "(Vector) operator ()", "" );
        return super_t::_data[ i * _stride ];
    }
    
    //! return reference to coefficient at position \a i
    value_t & operator () ( const idx_t  i ) noexcept
    {
        HASSERT( i < idx_t(_length), ERR_ARR_BOUND, "(Vector) operator ()", "" );
        return super_t::_data[ i * _stride ];
    }

    //
    // construction operators
    //

    //! create real copy of matrix
    Vector< value_t >  copy () const
    {
        Vector< value_t >  v( *this, copy_value );

        return v;
    }
    
    //! create reference to this matrix
    Vector< value_t >  reference () const
    {
        Vector< value_t >  v( *this, copy_reference );

        return v;
    }
    
    //! return reference to sub vector defined by \a r
    Vector    operator () ( const Range & r ) const
    {
        return Vector( *this, r, copy_reference );
    }

    //! give access to internal data
    value_t * data     () const noexcept { return super_t::_data; }
};

//
// trait for giving access to vector properties
//
template < typename T_value >
struct vector_trait< Vector< T_value > >
{
    using  value_t = T_value;
};

template < typename T_value >
struct is_vector< Vector< T_value > >
{
    static const bool  value = true;
};

//
// return real copy of given vector
//
template < typename T >
Vector< T >
copy ( const Vector< T > &  v )
{
    return v.copy();
}

//
// create random vector
//
template < typename T >
Vector< T >
random ( const size_t  length );

}// namespace BLAS

//////////////////////////////////////////////////////////
//
// debug helpers
//

namespace DBG
{

//!
//! write vector to file
//!
template <typename T>
void
write ( const BLAS::Vector< T > &  v,
        const std::string &        filename,
        const std::string &        matname );

}// namespace DBG

}// namespace HLIB

#endif  // __HLIB_BLAS_VECTOR_HH
