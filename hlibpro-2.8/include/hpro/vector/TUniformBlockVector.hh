#ifndef __HLIB_TUNIFORMBLOCKVECTOR_HH
#define __HLIB_TUNIFORMBLOCKVECTOR_HH
//
// Project     : HLib
// File        : TUniformBlockVector.hh
// Description : class for uniform block vector over a cluster basis
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/vector/TUniformVector.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TUniformBlockVector );

//!
//! \ingroup Vector_Module
//! \class   TUniformBlockVector 
//! \brief   Class for a uniform block vector, e.g. of _uniform_ sub blocks.
//!
//! In contrast to a TBlockVector, a TUniformBlockVector
//! stores only uniform sub vectors and, furthermore, also contains local
//! data. Hence, the data of the whole vector is not stored in the leaf
//! vectors alone, as in the case of a TBlockVector, but is distributed
//! through the complete hierarchy of a TUniformBlockVector.
//!
class TUniformBlockVector : public TUniformVector
{
private:
    // individual vector blocks
    std::vector< TUniformVector * >  _blocks;

public:
    //////////////////////////////////////////////////
    //
    // constructors and destructor
    //

    //! construct vector over index set \a ais, cluster basis \a acb,
    //! with coefficients \a acoeff and sub vectors defined in \a subvec
    TUniformBlockVector ( const TIndexSet &                        ais,
                          const TClusterBasis< real > *            acb,
                          const BLAS::Vector< real > &             acoeff,
                          const std::vector< TUniformVector * > &  asubvec )
            : TUniformVector( ais, acb, acoeff )
            , _blocks( asubvec )
    {}
    
    //! construct vector over index set \a ais, cluster basis \a acb,
    //! with coefficients \a acoeff and sub vectors defined in \a subvec
    TUniformBlockVector ( const TIndexSet &                        ais,
                          const TClusterBasis< complex > *         acb,
                          const BLAS::Vector< complex > &          acoeff,
                          const std::vector< TUniformVector * > &  asubvec )
            : TUniformVector( ais, acb, acoeff )
            , _blocks( asubvec )
    {}
    
    //! copy constructor
    TUniformBlockVector ( const TUniformBlockVector &  v )
            : TUniformVector( v )
    {
        assign( 1.0, & v );
    }

    //! destructor
    virtual ~TUniformBlockVector ()
    {
    }

    //////////////////////////////////////////////////
    //
    // access vector data
    //

    //! return number of blocks
    virtual uint            n_blocks    () const { return uint(_blocks.size()); }

    //! access single vector block
    TUniformVector *        block       ( const uint i )       { return _blocks[i]; }

    //! access single vector block
    const TUniformVector *  block       ( const uint i ) const { return _blocks[i]; }

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
    
    // //
    // // copy methods
    // //
    
    // //! copy from vector \a v
    // virtual void copy_from ( const TUniformBlockVector * v )
    // {
    //     HERROR( ERR_NOT_IMPL, "", "" );
    // }

    // //! copy to vector \a v
    // virtual void copy_to ( TUniformBlockVector * v ) const
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
    // virtual void     conjugate  ();
    
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

    //
    // memory consumption
    //
    
    //! return size in bytes used by this object
    virtual size_t byte_size () const;

    //
    // virtual constructor
    //

    //! return copy of vector
    virtual auto  copy   () const -> std::unique_ptr< TVector > { return std::make_unique< TUniformBlockVector >( *this ); }
    
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


    HLIB_RTTI_DERIVED( TUniformBlockVector, TUniformVector );
};

}// namespace HLIB

#endif  // __HLIB_TUNIFORMBLOCKVECTOR_HH
