#ifndef __HLIB_TUNIFORMVECTOR_HH
#define __HLIB_TUNIFORMVECTOR_HH
//
// Project     : HLib
// File        : TUniformVector.hh
// Description : class for uniform vector over a cluster basis
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/cluster/TClusterBasis.hh"
#include "hpro/vector/TVector.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TUniformVector );

//!
//! \ingroup Vector_Module
//! \class   TUniformVector 
//! \brief   Class for a uniform vector, e.g. represented as \f$x = V s\f$ with
//!          cluster basis \f$V\f$ and coefficients \f$s\f$.
//!
class TUniformVector : public TVector
{
private:
    // cluster basis
    const TClusterBasis< real > *     _rbasis;
    const TClusterBasis< complex > *  _cbasis;
    
    // vector coefficients
    BLAS::Vector< real >              _rcoeff;
    BLAS::Vector< complex >           _ccoeff;

    // size of vector
    size_t                            _size;
    
public:
    //////////////////////////////////////////////////
    //
    // constructors and destructor
    //

    //! construct vector over index set \a ais, cluster basis \a acb
    //! and with coefficients \a acoeff
    TUniformVector ( const TIndexSet &                 ais,
                     const TClusterBasis< real > *     acb,
                     const BLAS::Vector< real > &      acoeff )
            : TVector( ais.first(), real_valued ),
              _rbasis( acb ),
              _cbasis( nullptr ),
              _rcoeff( acoeff, copy_value ),
              _size( ais.size() )
    {}
    
    //! construct vector over index set \a ais, cluster basis \a acb
    //! and with coefficients \a acoeff
    TUniformVector ( const TIndexSet &                 ais,
                     const TClusterBasis< complex > *  acb,
                     const BLAS::Vector< complex > &   acoeff )
            : TVector( ais.first(), complex_valued ),
              _rbasis( nullptr ),
              _cbasis( acb ),
              _ccoeff( acoeff, copy_value ),
              _size( ais.size() )
    {}
    
    //! copy constructor
    TUniformVector ( const TUniformVector &  v )
            : TVector(),
              _size(0)
    {
        assign( 1.0, & v );
    }

    //! destructor
    virtual ~TUniformVector ()
    {}

    //////////////////////////////////////////////////
    //
    // access vector data
    //

    //! return size of vector
    virtual size_t  size        () const
    {
        return _size;
    }

    //! return rank of vector, e.g. rank of cluster basis
    virtual size_t  rank        () const
    {
        if ( is_complex() ) return _cbasis->rank();
        else                return _rbasis->rank();
    }

    //! access coefficent \a i (real valued)
    virtual real            entry       ( const idx_t ) const
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        return 0.0;
    }
    
    //! access coefficent \a i (complex valued)
    virtual const complex   centry      ( const idx_t ) const
    {
        HERROR( ERR_NOT_IMPL, "", "" );
        return 0.0;
    }
    
    //! return basis coefficients
    const TClusterBasis< real > *     rbasis () const { return _rbasis; }
    const TClusterBasis< complex > *  cbasis () const { return _cbasis; }
    
    //! return basis coefficients
    BLAS::Vector< real > &            rcoeff ()       { return _rcoeff; }
    BLAS::Vector< complex > &         ccoeff ()       { return _ccoeff; }
    const BLAS::Vector< real > &      rcoeff () const { return _rcoeff; }
    const BLAS::Vector< complex > &   ccoeff () const { return _ccoeff; }
    
    //
    // handle index set
    //
    
    //! set size of vector
    virtual void set_size    ( const size_t  n );
    
    //! define vector by cluster
    virtual void set_cluster ( const TCluster *  c )
    {
        if ( c != nullptr )
        {
            set_is( *c );
        }// if
    }

    //! define vector by indexset
    virtual void set_is ( const TIndexSet &  ais )
    {
        set_size( ais.size() );
        set_ofs( ais.first() );
    }

    // //
    // // copy methods
    // //
    
    // //! copy from vector \a v
    // virtual void copy_from ( const TUniformVector * )
    // {
    //     HERROR( ERR_NOT_IMPL, "", "" );
    // }

    // //! copy to vector \a v
    // virtual void copy_to ( TUniformVector * ) const
    // {
    //     HERROR( ERR_NOT_IMPL, "", "" );
    // }
    
protected:

    //
    // switch between real and complex
    //
    
    //! switch to real valued representation if possible
    virtual void to_real    ();

    //! switch to complex valued representation
    virtual void to_complex ();

public:
    //////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    //! fill vector with constant \a α
    virtual void   fill        ( const real alpha );

    //! fill vector with random numbers
    virtual void   fill_rand   ( const uint seed );

    //! set this ≔ α · this
    virtual void   scale       ( const real alpha );

    //! set this ≔ α · x
    virtual void   assign      ( const real alpha, const TVector * x );

    //! compute ‖·‖₂
    virtual real   norm2       () const;

    //! compute ‖·‖∞
    virtual real   norm_inf    () const;
    
    //! set this ≔ this + α · x
    virtual void   axpy        ( const real alpha, const TVector * x );

    //////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    //! conjugate coefficients
    virtual void     conjugate  ();
    
    //! fill vector with constant \a α
    virtual void     cfill      ( const complex &  alpha );

    //! set this ≔ α · this
    virtual void     cscale     ( const complex &  alpha );

    //! set this ≔ α · x
    virtual void     cassign    ( const complex &  alpha, const TVector * x );

    //! return inner product <this, x> = this^H x
    virtual complex  dot        ( const TVector * x ) const;

    //! return inner product <this, x> = this^T x
    virtual complex  dotu       ( const TVector * x ) const;

    //! set this ≔ this + α · x
    virtual void     caxpy      ( const complex & f, const TVector * x );

    //////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // virtual constructor
    //

    //! return copy of vector
    virtual auto  copy   () const -> std::unique_ptr< TVector > { return std::make_unique< TUniformVector >( *this ); }
    
    //! return object of same class
    virtual auto  create () const -> std::unique_ptr< TVector > { return nullptr; }
    
    //
    // restriction
    //

    //! return vector restricted to real part of coefficients
    virtual auto restrict_re () const -> std::unique_ptr< TVector >;

    //! return vector restricted to imaginary part of coefficients
    virtual auto restrict_im () const -> std::unique_ptr< TVector >;
    
    //
    // serialisation
    //

    //! read vector from stream
    virtual void read  ( TByteStream & s );

    //! write vector to stream
    virtual void write ( TByteStream & s ) const;

    //! returns size of object in bytestream
    virtual size_t bs_size () const;

    //
    // parallel methods
    //

    //! pointwise summation between all vectors in \a ps
    virtual void sum ( const TProcSet & ps );
    using TVector::sum;


    HLIB_RTTI_DERIVED( TUniformVector, TVector );
};

//
// wrapper to rbasis/cbasis
//
template <typename T>
const TClusterBasis< T > *
basis            ( const TUniformVector *  v );

template <>
inline
const TClusterBasis< real > *
basis< real >    ( const TUniformVector *  v )
{
    return v->rbasis();
}

template <>
inline
const TClusterBasis< complex > *
basis< complex > ( const TUniformVector *  v )
{
    return v->cbasis();
}

//
// wrapper to rcoeff/ccoeff
//
template <typename T>
const BLAS::Vector< T > &
blas_coeff            ( const TUniformVector *  v );

template <>
inline
const BLAS::Vector< real > &
blas_coeff< real >    ( const TUniformVector *  v )
{
    return v->rcoeff();
}

template <>
inline
const BLAS::Vector< complex > &
blas_coeff< complex > ( const TUniformVector *  v )
{
    return v->ccoeff();
}

}// namespace HLIB

#endif  // __HLIB_TUNIFORMVECTOR_HH
