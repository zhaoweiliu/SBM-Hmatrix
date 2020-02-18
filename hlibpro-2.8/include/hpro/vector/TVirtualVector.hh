#ifndef __HLIB_TVIRTUALVECTOR_HH
#define __HLIB_TVIRTUALVECTOR_HH
//
// Project     : HLib
// File        : TVirtualVector.hh
// Description : class for a vector which has no own data
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/vector/TScalarVector.hh"

namespace HLIB
{

//////////////////////////////////////////////////////////
//
// real valued virtual scalar vector
//

// local matrix type
DECLARE_TYPE( TVirtualVector );

//!
//! \ingroup Vector_Module
//! \class   TVirtualVector
//! \brief   A virtual vector gets his data from some real vector and behaves just like it,
//!          except memory-management.
//!
class TVirtualVector : public TScalarVector
{
public:
    ////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct virtual vector without data
    TVirtualVector ()
    {}

    //! copy constructor
    TVirtualVector ( const TVirtualVector & v )
            : TScalarVector()
    {
        set_vector( & v, v.is() );
    }
    
    //! reference data in \a vec with additional offset (real valued)
    TVirtualVector ( const BLAS::Vector< real > &  vec,
                     const TIndexSet &             ais )
    {
        set_vector( vec, ais );
    }

    //! reference data in \a vec with additional offset (complex valued)
    TVirtualVector ( const BLAS::Vector< complex > &  vec,
                     const TIndexSet &                ais )
    {
        set_vector( vec, ais );
    }

    //! reference part of vector vector \a x defined by index set \a ais
    TVirtualVector ( const TScalarVector *  x,
                     const TIndexSet &      ais )
    {
        set_vector( x, ais );
    }
        
    //! destructor
    ~TVirtualVector () {}

    /////////////////////////////////////////
    //
    // manipulate array
    //

    //! reference vector \a vec (real valued)
    virtual void set_vector ( const BLAS::Vector< real > &  vec )
    {
        _complex = false;
        _rvec    = vec;
        _size    = vec.length();
        set_ofs( 0 );
    }

    //! reference vector \a vec (complex valued)
    virtual void set_vector ( const BLAS::Vector< complex > &  vec )
    {
        _complex = true;
        _cvec    = vec;
        _size    = vec.length();
        set_ofs( 0 );
    }

    //! reference part of vector \a vec defined by index set \a ais (real valued)
    virtual void set_vector ( const BLAS::Vector< real > &  vec, const TIndexSet &  ais )
    {
        _complex = false;
        _rvec    = BLAS::Vector< real >( vec, BLAS::Range( 0, idx_t(ais.size()) - 1 ) );
        _size    = ais.size();
        set_ofs( ais.first() );
    }

    //! reference part of vector \a vec defined by index set \a ais (complex valued)
    virtual void set_vector ( const BLAS::Vector< complex > &  vec, const TIndexSet &  ais )
    {
        _complex = true;
        _cvec    = BLAS::Vector< complex >( vec, BLAS::Range( 0, idx_t(ais.size()) - 1 ) );
        _size    = ais.size();
        set_ofs( ais.first() );
    }

    //! reference part of vector \a vec defined by index set \a ais
    virtual void set_vector ( const TScalarVector *  vec, const TIndexSet &  ais );
    using TScalarVector::set_vector;

    //! set size of the vector (must not differ from current size in order to not reallocate memory)
    virtual void set_size ( const size_t  n )
    {
        if ( n != _size )
            HERROR( ERR_CONSISTENCY, "(TVirtualVector) set_size", "must not change size in virtual vector" );
        
        _size = n;
    }
    
    //! switch to complex representation
    virtual void set_complex ( const bool b )
    {
        if ( b == is_complex() )
            return;
        else
            HERROR( ERR_CONSISTENCY, "(TVirtualVector) set_complex", "can not change field type" );
    }
    
    //! copy operator (copy pointers !!!)
    virtual TVirtualVector & operator = ( const TVirtualVector &  vec )
    {
        _complex = vec.is_complex();
        _rvec    = vec.blas_rvec();
        _cvec    = vec.blas_cvec();
        _size    = vec.size();
        
        set_ofs( vec.ofs() );
        
        return *this;
    }

    //! compute this = α·x
    virtual void assign  ( const real      alpha, const TVector * x );

    //! compute this = α·x
    virtual void cassign ( const complex & alpha, const TVector * x );
    
    //////////////////////////////////////////////////
    //
    // misc. methods
    //

    //
    // virtual constructor
    //

    //! return copy of vector (real copy, not a virtual one)
    virtual auto  copy   () const -> std::unique_ptr< TVector > { return TScalarVector::copy(); }
    
    //! return object of same class
    virtual auto  create () const -> std::unique_ptr< TVector > { return std::make_unique< TVirtualVector >(); }

    HLIB_RTTI_DERIVED( TVirtualVector, TScalarVector );
};

}// namespace

#endif  // __HLIB_TVIRTUALVECTOR_HH
