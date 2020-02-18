#ifndef __HLIB_TCLUSTERBASIS_HH
#define __HLIB_TCLUSTERBASIS_HH
//
// Project     : HLib
// File        : TClusterBasis.hh
// Description : classes for representing a cluster basis
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <vector>
#include <utility>

#include "hpro/base/TTypeInfo.hh"
#include "hpro/cluster/TCluster.hh"
#include "hpro/cluster/TBlockCluster.hh"
#include "hpro/blas/Matrix.hh"
#include "hpro/blas/Algebra.hh"
#include "hpro/base/TTruncAcc.hh"
#include "hpro/vector/TVector.hh"

namespace HLIB
{

// forward decl.
class TUniformVector;

// local type
DECLARE_TYPE( TClusterBasis );

//!
//! \ingroup  Cluster_Module
//! \class    TClusterBasis
//! \brief    class representing a nested cluster basis
//!
template < typename T >
class TClusterBasis : public TIndexSet, public TTypeInfo
{
public:
    //! value type of cluster basis
    using  value_t    = T;

    //! iterator for sons
    using  son_iter_t = typename std::vector< TClusterBasis * >::iterator;
    
private:
    //!@cond
    
    //! local basis (usually in case of a leaf)
    BLAS::Matrix< value_t >                 _V;

    //! transfer matrices from the sons
    std::vector< BLAS::Matrix< value_t > >  _E;

    //! rank of the basis
    uint                                    _rank;
    
    //! the sons
    std::vector< TClusterBasis * >          _sons;

    //!@endcond
    
public:
    ///////////////////////////////////////////////////////
    //
    // ctor and dtor
    //

    //! construct cluster basis corresponding to cluster \a cl
    //! with basis defined by \a V
    TClusterBasis ( const TIndexSet &                            is,
                    const BLAS::Matrix< T > &                    V );

    //! construct cluster basis corresponding to cluster \a cl
    //! with son bases \a asons and transfer matrices \a aE
    TClusterBasis ( const TIndexSet &                            is,
                    const std::vector< TClusterBasis< T > * > &  asons,
                    const std::vector< BLAS::Matrix< T > > &     aE );

    //! construct cluster basis corresponding to cluster \a cl
    //! with son bases \a asons (no transfer matrices)
    TClusterBasis ( const TIndexSet &                            is,
                    const std::vector< TClusterBasis< T > * > &  asons );

    //! dtor: delete all sons
    virtual ~TClusterBasis ();

    ///////////////////////////////////////////////////////
    //
    // tree management
    //

    //! return no of nodes
    virtual size_t  nnodes () const;

    //! return depth of tree
    virtual size_t  depth () const;

    ///////////////////////////////////////////////////////
    //
    // basis management
    //

    //! return basis rank
    uint                   rank    () const                  { return _rank; }

    //! return number of sons
    size_t                 nsons   () const                  { return _sons.size(); }

    //! return \a i'th son basis
    TClusterBasis *        son     ( const size_t  i )       { return _sons[i]; }
    const TClusterBasis *  son     ( const size_t  i ) const { return _sons[i]; }

    //! access to explicit local basis
    const BLAS::Matrix< T > &  basis        ()                const { return _V; }
    
    //! access to transfer matrices
    const BLAS::Matrix< T > &  transfer_mat ( const uint  i ) const { return _E[i]; }
    
    ///////////////////////////////////////////////////////
    //
    // basis functionality
    //

    //! truncate basis as defined by accuracy \a acc
    void               truncate        ( const TTruncAcc &  acc );

    //! transfer given data to basis of \a i'th son and return result
    BLAS::Vector< T >  transfer_to_son ( const uint                 i,
                                         const BLAS::Vector< T > &  s ) const
    {
        BLAS::Vector< T >  s_i;

        if ( i >= nsons() )
            HERROR( ERR_ARR_BOUND, "(TClusterBasis) transfer_to_son", "" );
        
        if ( s.length() > 0 )
        {
            s_i = BLAS::Vector< T >( _E[i].nrows() );

            BLAS::mulvec( T(1), _E[i], s, T(0), s_i );
        }// if
        
        return s_i;
    }
    
    ///////////////////////////////////////////////////////
    //
    // vector management
    //

    //! construct uniform vector corresponding to cluster basis
    TUniformVector *  build_vec () const;
    
    ///////////////////////////////////////////////////////
    //
    // vector and matrix transformation
    //

    //
    // forward transformation
    //
    
    //! forward transformation: s ≔ V^H · v
    void              transform_forward   ( const BLAS::Vector< T > &          v,
                                            BLAS::Vector< T > &                s ) const;
    
    //! forward transformation: return uniform vector
    TUniformVector *  transform_forward   ( const TVector *                    v ) const;

    //! block forward transformation: S ≔ V^H · M
    template < typename T_mat >
    void              transform_forward   ( const BLAS::MatrixBase< T_mat > &  M,
                                            BLAS::Matrix< T > &                S ) const;
    //
    // backward transformation
    //
    
    //! backward transformation: v ≔ V · s
    void              transform_backward  ( const BLAS::Vector< T > &          s,
                                            BLAS::Vector< T > &                v ) const;

    //! backward transformation: return transformed vector
    TVector *         transform_backward  ( const TUniformVector *             v ) const;

    //! block backward transformation: M ≔ V · S
    template < typename T_mat >
    void              transform_backward  ( const BLAS::MatrixBase< T_mat > &  S,
                                            BLAS::Matrix< T > &                M ) const;

    ////////////////////////////////////////////////////////
    //
    // misc. methods
    //

    HLIB_RTTI_BASE( TClusterBasis );

    //! return copy of basis
    virtual auto    copy       () const -> std::unique_ptr< TClusterBasis< T > >;
    
    //! return size in bytes used by this object
    virtual size_t  byte_size  () const;

    //! print basis to terminal
    virtual void    print      ( const uint  ofs = 0 ) const;
};

////////////////////////////////////////////////////////
//
// non instantiable methods
//

template < typename T >
template < typename T_derived >
inline
void
TClusterBasis<T>::transform_forward  ( const BLAS::MatrixBase< T_derived > &  M,
                                       BLAS::Matrix< T > &                    S ) const
{
    S.resize( rank(), M.ncols() );
    
    if ( nsons() == 0 )
    {
        prod( T(1), adjoint( _V ), M, T(0), S );
    }// if
    else
    {
        //
        // S ≔ V^H·M = (∑_i V_i E_i)^H M = ∑_i E_i^H V_i^H M
        //
        
        idx_t  son_ofs = 0;
        
        BLAS::fill( T(0), S );
        
        for ( size_t  i = 0; i < nsons(); ++i )
        {
            TClusterBasis< T > *  son_i = _sons[i];
            BLAS::Matrix< T >     Si;

            // restrict M to son cluster
            // TO DO: do not copy, just map
            const size_t       nrows_i = son_i->size();
            BLAS::Matrix< T >  Mi( nrows_i, M.ncols() );

            for ( idx_t  jj = 0; jj < idx_t(M.ncols()); ++jj )
                for ( idx_t  ii = 0; ii < idx_t(nrows_i); ++ii )
                    Mi( ii, jj ) = M( ii + son_ofs, jj );
                      
            son_i->transform_forward( Mi, Si );
            BLAS::prod( T(1), adjoint(_E[i]), Si, T(1), S );

            son_ofs += idx_t( nrows_i );
        }// for
    }// else
}

template < typename T >
template < typename T_derived >
inline
void
TClusterBasis<T>::transform_backward  ( const BLAS::MatrixBase< T_derived > &  S,
                                        BLAS::Matrix< T > &                    M ) const
{
    M.resize( size(), S.ncols() );
    
    if ( nsons() == 0 )
    {
        BLAS::prod( T(1), _V, S, T(0), M );
    }// if
    else
    {
        //
        // M ≔ V·S = ∑_i V_i·E_i·S = ∑_i transform_backward( son_i, E_i·S )
        //
        
        idx_t  son_ofs = 0;
        
        for ( size_t  i = 0; i < nsons(); ++i )
        {
            TClusterBasis< T > *  son_i = _sons[i];
            const size_t          nrows_i = son_i->size();
            BLAS::Matrix< T >     S_i( _E[i].nrows(), S.ncols() );

            BLAS::prod( T(1), _E[i], S, T(0), S_i );

            BLAS::Matrix< T >  M_i( M,
                                    BLAS::Range( son_ofs, son_ofs + idx_t( nrows_i ) - 1 ),
                                    BLAS::Range( 0, idx_t( S_i.ncols() )-1 ) );

            son_i->transform_backward( S_i, M_i );
            
            son_ofs += idx_t( nrows_i );
        }// for
    }// else
}

}// namespace HLIB

#endif  // __HLIB_TCLUSTERBASIS_HH
