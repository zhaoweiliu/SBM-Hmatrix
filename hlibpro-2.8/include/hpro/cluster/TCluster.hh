#ifndef __HLIB_TCLUSTER_HH
#define __HLIB_TCLUSTER_HH
//
// Project     : HLib
// File        : TCluster.hh
// Description : baseclass for a node in a cluster tree
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <vector>
#include <list>

#include "hpro/base/System.hh"
#include "hpro/base/TTypeInfo.hh"
#include "hpro/cluster/TIndexSet.hh"

namespace HLIB
{

// local type
DECLARE_TYPE( TCluster );

//!
//! \ingroup  Cluster_Module
//! \class    TCluster
//! \brief    Represents a node in a cluster tree with an arbitrary number of sons
//!           - also for nested dissection case with interface clusters
//!
class TCluster : public TIndexSet, public TTypeInfo
{
protected:
    //!@cond
    
    //! son-clusters
    std::vector< TCluster * >  _sons;
    
    //! true if cluster is domain cluster
    bool                       _is_domain;

    //!@endcond
    
public:
    ////////////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct node with empty index set
    TCluster ()
            : _is_domain(false)
    {}

    //! construct node with index set \a is
    TCluster ( const TIndexSet &  is )
            : TIndexSet( is )
            , _is_domain(false)
    {}

    //! construct node with index set [\a first_idx, \a last_idx]
    TCluster ( const idx_t  first_idx,
               const idx_t  last_idx )
            : TIndexSet(first_idx,last_idx)
            , _is_domain(false)
    {}

    //! dtor
    virtual ~TCluster ();
    
    ////////////////////////////////////////////////////////
    //
    // access internal data
    //

    //! return true if node is domain cluster
    bool is_domain  () const         { return _is_domain; }

    //! set domain status of node
    void set_domain ( const bool b ) { _is_domain = b; }

    ////////////////////////////////////////////////
    //
    // manage tree
    //

    //! return number of sons
    virtual uint             nsons       () const { return uint(_sons.size()); }

    //! set number of sons
    virtual void             set_nsons   ( const uint  n );
    
    //! return i'th son
    virtual TCluster *       son         ( const uint  i )       { return _sons[i]; }

    //! return i'th son
    virtual const TCluster * son         ( const uint  i ) const { return _sons[i]; }

    //! set i'th son. If \a del is true, former son_i is deleted
    virtual void             set_son     ( const uint  i,
                                           TCluster *    son,
                                           const bool    del = true );

    //! add a son (gets first unused slot); if \a inc_nsons is true,
    //! the number of sons will be increased, if no free slot is available
    virtual void             add_son     ( TCluster *  son,
                                           const bool  inc_nsons = false );

    //! remove references to sons, no deletion
    virtual void             clear_sons  () { set_nsons( 0 ); }
    
    //! return true if node is leaf
    virtual bool             is_leaf     () const
    {
        if ( nsons() == 0 )
            return true;
        
        for ( uint  i = 0; i < nsons(); i++ )
            if ( _sons[i] != NULL )
                return false;
        
        return true;
    }

    //! return no of nodes
    virtual uint  nnodes  () const;

    //! return depth of tree
    virtual uint  depth   () const;

    ////////////////////////////////////////////////////////
    //
    // misc. methods
    //

    //! collect leaves (or nodes with depth \a depth) in \a list
    virtual void       collect_leaves ( std::list< TCluster * > &  leaves,
                                        const int                  depth = -1,
                                        const int                  level = 0 ) const;
    
    //! return object of same type
    virtual TCluster * create () const { return new TCluster(); }

    //! return copy of node/subtree
    virtual TCluster * copy  () const;

    //! stream output
    virtual void       print ( const uint ofs = 0 ) const;

    //! return size in bytes used by this object
    virtual size_t     byte_size () const;

    ///////////////////////////////////////////////////////////
    //
    // RTTI
    //

    HLIB_RTTI_BASE( TCluster );
};

//!
//! \ingroup  Cluster_Module
//! \brief    flatten hierarchy of cluster tree
//! \detail   Flatten hierarchy of cluster tree starting at \a root 
//!
void flatten ( TCluster * root ); 

}// namespace HLIB

#endif  // __HLIB_TCLUSTER_HH
