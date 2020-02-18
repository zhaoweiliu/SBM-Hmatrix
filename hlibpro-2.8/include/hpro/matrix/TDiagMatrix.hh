#ifndef __HLIB_TDIAGMATRIX_HH
#define __HLIB_TDIAGMATRIX_HH
//
// Project     : HLib
// File        : TDiagMatrix.hh
// Description : class for a diagonal-matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TScalarVector.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TDiagMatrix );

//
// class for a diagonal matrix
//
class TDiagMatrix : public TMatrix
{
private:
    // data and size-info (nxn format !!)
    TScalarVector   _diag;

public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    TDiagMatrix () {}
    TDiagMatrix ( const uint n ) : _diag( n, 0 ) {}
    TDiagMatrix ( const TBlockCluster * bct )
            : TMatrix(bct)
    { set_cluster( bct ); }

    virtual ~TDiagMatrix () {}
    
    ////////////////////////////////////////////////
    //
    // access internal data
    //

    virtual size_t  rows () const { return _diag.size(); }
    virtual size_t  cols () const { return _diag.size(); }
    
    void set_size ( const size_t n, const size_t m )
    {
        if ( n != m )
            HERROR( ERR_MAT_STRUCT, "(TDiagMatrix) set_size", "only support square matrices" );

        _diag.set_size( n );
    }
    
    // handle size (via cluster)
    virtual void set_cluster ( const TBlockCluster * c );

    // return pointer to diagonal (DO NOT DELETE !!!)
    BLAS::Vector< real > &           blas_rdiag  ()       { return _diag.blas_rvec(); }
    const BLAS::Vector< real > &     blas_rdiag  () const { return _diag.blas_rvec(); }
    BLAS::Vector< complex > &        blas_cdiag ()        { return _diag.blas_cvec(); }
    const BLAS::Vector< complex > &  blas_cdiag () const  { return _diag.blas_cvec(); }
    
    // switch between complex and real format
    virtual void to_real    () { _diag.set_complex( false ); }
    virtual void to_complex () { _diag.set_complex( true ); }
    
    ///////////////////////////////////////////////
    //
    // access coeff.
    //

    real          entry  ( const idx_t i, const idx_t j ) const { return ( i == j ? _diag.entry(i)  : real(0)); }
    const complex centry ( const idx_t i, const idx_t j ) const { return ( i == j ? _diag.centry(i) : real(0)); }

    void set_entry ( const idx_t i, const idx_t j, const real f )
    {
        if ( i != j )
            HERROR( ERR_DIAG_ENTRY, "(TDiagMatrix) entry", "" );
        
        return _diag.set_entry( i, f );
    }
    
    void set_centry ( const idx_t i, const idx_t j, const complex f )
    {
        if ( i != j )
            HERROR( ERR_DIAG_ENTRY, "(TDiagMatrix) centry", "" );
        
        return _diag.set_centry( i, f );
    }
    
    /////////////////////////////////////////////////
    //
    // BLAS-routines (real valued)
    //

    // scale matrix by constant factor
    virtual void scale ( const real f ) { _diag.scale( f ); }
    
    // matrix-vector-multiplication
    virtual void mul_vec ( const real      alpha,
                           const TVector * x,
                           const real      beta,
                           TVector       * y,
                           const matop_t   op = MATOP_NORM ) const;

    /////////////////////////////////////////////////
    //
    // BLAS-routines (complex valued)
    //

    // scale matrix by constant factor
    virtual void cscale ( const complex f ) { _diag.cscale( f ); }
    
    // matrix-vector-multiplication
    virtual void cmul_vec ( const complex   alpha,
                            const TVector * x,
                            const complex   beta,
                            TVector       * y,
                            const matop_t   op = MATOP_NORM ) const;

    /////////////////////////////////////////////////
    //
    // misc.
    //

    // transpose matrix
    virtual void transpose ()
    {}
    
    // conjugate matrix coefficients
    virtual void conjugate ()
    {
        _diag.conjugate();
    }
    
    // truncate matrix to given accuracy
    virtual void truncate ( const TTruncAcc & )
    {}
    
    // virtual constructors
    virtual auto  create  () const -> std::unique_ptr< TMatrix > { return std::make_unique< TDiagMatrix >(); }
    virtual auto  copy    () const -> std::unique_ptr< TMatrix >;
    virtual void  copy_to ( TMatrix * A ) const;
    using TMatrix::copy;
    using TMatrix::copy_to;

    // return size in bytes used by this object
    virtual size_t byte_size () const { return TMatrix::byte_size() + _diag.byte_size(); }

    //
    // RTTI
    //

    HLIB_RTTI_DERIVED( TDiagMatrix, TMatrix )

    //
    // serialisation
    //

    virtual void read  ( TByteStream & s );
    virtual void build ( TByteStream & s );
    virtual void write ( TByteStream & s ) const;

    // returns size of object in bytestream
    virtual size_t bs_size () const;
};

}// namespace

#endif  // __HLIB_TDIAGMATRIX_HH
