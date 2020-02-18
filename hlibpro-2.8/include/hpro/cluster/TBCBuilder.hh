#ifndef __HLIB_TBCBUILDER_HH
#define __HLIB_TBCBUILDER_HH
//
// Project     : HLib
// File        : TBCBuilder.hh
// Description : builds a block-cluster-tree
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/cluster/types.hh"
#include "hpro/cluster/TBlockCluster.hh"
#include "hpro/cluster/TBlockClusterTree.hh"
#include "hpro/cluster/TClusterTree.hh"
#include "hpro/cluster/TAdmCondition.hh"
#include "hpro/parallel/TProcSet.hh"

namespace HLIB
{

//!
//! \ingroup  Cluster_Module
//! \class    TBCBuilder
//! \brief    Recursively build block cluster tree with supplied admissibility condition.
//!
class TBCBuilder
{
protected:
    //! @cond
    
    //! minimal level to accept admissible clusters
    const uint  _min_lvl;

    //! disallow/allow different cluster levels in block cluster
    const bool  _same_cluster_level;
    
    //! @endcond
    
public:
    ////////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct block cluster tree builder 
    TBCBuilder ( const uint                  min_lvl     = 0,
                 const cluster_level_mode_t  cl_lvl_mode = CFG::Cluster::cluster_level_mode );

    //! dtor
    virtual ~TBCBuilder () {}

    ////////////////////////////////////////////////
    //
    // build a block-clustertree
    //

    //! build block cluster tree for \a rowct × \a colct with
    //! admissibility condition \a ac
    virtual std::unique_ptr< TBlockClusterTree >  build ( const TClusterTree   * rowct,
                                                          const TClusterTree   * colct,
                                                          const TAdmCondition  * ac ) const;

    //! refine given block cluster
    virtual void refine     ( TBlockCluster *  bc ) const;
    
protected:
    //! recusivly build a block cluster tree for \a rowcl × \a colcl
    virtual void  rec_build ( TBlockCluster *        bc,
                              const TAdmCondition *  ac,
                              const uint             level ) const;

    //! create new node in tree with father node \a parent
    virtual std::unique_ptr< TBlockCluster >  create_bc ( TBlockCluster *  parent ) const;

    //! create new node in tree with father node \a parent and given row/col clusters
    virtual std::unique_ptr< TBlockCluster >  create_bc ( TBlockCluster *  parent,
                                                          TCluster *       rowcl,
                                                          TCluster *       colcl ) const;
};

}// namespace HLIB

#endif  // __HLIB_TBCBUILDER_HH
