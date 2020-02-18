#ifndef __HLIB_TBSPCTBUILDER_HH
#define __HLIB_TBSPCTBUILDER_HH
//
// Project     : HLib
// File        : TBSPCTBuilder.hh
// Description : build clustertrees for coordinate based indexsets via BSP
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <vector>

#include "hpro/base/config.hh"

#include "hpro/cluster/types.hh"
#include "hpro/cluster/TNodeSet.hh"
#include "hpro/cluster/TGeomCluster.hh"
#include "hpro/cluster/TClusterTree.hh"
#include "hpro/cluster/TBSPPartStrat.hh"
#include "hpro/cluster/TCoordinate.hh"
#include "hpro/cluster/TPermutation.hh"

#include "hpro/matrix/TSparseMatrix.hh"

namespace HLIB
{

//!
//! \ingroup  Cluster_Module
//! \class    TGeomCTBuilder
//! \brief    Base class for all cluster tree constructors based on geometry data.
//!
class TGeomCTBuilder
{
protected:
    //!
    //! \struct  data_t
    //! \brief   Datatype for internal argument transfer
    //!
    struct data_t
    {
        const TCoordinate *  coord;        // coordinate information of indices
        TPermutation *       perm;         // permutation to built
        const uint           nmin;         // minimale leaf size
        const uint           min_leaf_lvl; // minimal level where leaves might appear
        const uint           max_lvl;      // maximal level to reach in clustering
    };

    //!
    //! \class  TOptClusterSize
    //! \brief  Controls optimal cluster size per tree level.
    //!
    class TOptClusterSize
    {
    private:
        //! optimal size of cluster 
        size_t   _opt_size;

        //! optimal reduction of indices per level
        double   _reduction;

    public:
        //! construct size control
        TOptClusterSize ( const size_t  opt_size  = 0,
                          const double  reduction = 0.0 )
                : _opt_size(opt_size), _reduction( reduction )
        {}

        //! return true, if \a n is of optimal size
        bool  is_optimal ( const size_t  n ) const { return n >= _opt_size; }

        //! return size control for next level in tree
        TOptClusterSize  recurse () const
        {
            return TOptClusterSize( size_t( double(_opt_size) * _reduction ),
                                    _reduction );
        }
    };

protected:
    //! minimal size of a cluster, i.e. not smaller than this
    uint                   _n_min;
 
    //! minimal level on which leaves may occur
    uint                   _min_leaf_lvl;

    //! flag for adjusting bounding boxes of nodes
    bool                   _adjust_bb;

    //! flag for sorting sub clusters w.r.t. size
    bool                   _sort_wrt_size;
    

public:
    //! construct cluster tree builder
    TGeomCTBuilder ( const uint  n_min        = CFG::Cluster::nmin,
                     const uint  min_leaf_lvl = 0 );

    // dtor
    virtual ~TGeomCTBuilder () {}
    
    //!
    //! build cluster tree out of given coordinate set
    //! \param   coord   : geometry information for each index
    //! \param   idx_ofs : start renumbering indices from \a idx_ofs
    //!
    virtual std::unique_ptr< TClusterTree >  build  ( const TCoordinate *      coord,
                                                      const idx_t              idx_ofs = 0 ) const;

    //! recursively build cluster tree for indices in \a dofs
    virtual std::unique_ptr< TGeomCluster >  divide ( const TNodeSet &         dofs,
                                                      const uint               lvl,
                                                      const TBBox &            bbox,
                                                      const TOptClusterSize &  csize,
                                                      const idx_t              index_ofs,
                                                      data_t &                 data ) const = 0;

    //
    // give access to local parameters
    //

    uint  n_min          () const { return _n_min;         }
    uint  min_leaf_lvl   () const { return _min_leaf_lvl;  }
    bool  adjust_bb      () const { return _adjust_bb;     }
    bool  sort_wrt_size  () const { return _sort_wrt_size; }
    
    //! set flag for adjusting bounding box
    void  adjust_bb      ( const bool  b ) { _adjust_bb = b; }

    //! set flag for sorting son cluster wrt. size
    void  sort_wrt_size  ( const bool  b ) { _sort_wrt_size = b; }
    
protected:
    //! create a leaf in a clustertree containing indices in \a dofs
    virtual std::unique_ptr< TGeomCluster >  build_leaf ( const TNodeSet & dofs,
                                                          const uint       lvl,
                                                          const idx_t      index_ofs,
                                                          const TBBox &    bbox,
                                                          data_t &         data ) const;

    //! compute bounding box of index set defined by \a dofs
    virtual TBBox  compute_bb   ( const TNodeSet & dofs,
                                  const data_t &   data ) const;

    //! update bounding box of index set defined by \a dofs
    virtual void   update_bb    ( const TNodeSet & dofs,
                                  TBBox &          bbox,
                                  const data_t &   data ) const;

    //! check and update bbox in case of degenerate axis, e.g. very small length
    virtual void   check_bb     ( TBBox &          bbox,
                                  const data_t &   data ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TBSPCTBuilder
//! \brief    Base class for all cluster tree constructors based on BSP.
//!
class TBSPCTBuilder : public TGeomCTBuilder
{
protected:
    //! type of partitioning strategy
    const TBSPPartStrat *  _part_strat;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct BSP cluster tree builder with partitioning strategy \a part_strat
    TBSPCTBuilder ( const TBSPPartStrat *  part_strat,
                    const uint             n_min        = CFG::Cluster::nmin,
                    const uint             min_leaf_lvl = 0 );

    //! dtor
    virtual ~TBSPCTBuilder ();

    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //! recursively build cluster tree for indices in \a dofs
    virtual std::unique_ptr< TGeomCluster >  divide ( const TNodeSet &         dofs,
                                                      const uint               lvl,
                                                      const TBBox &            bbox,
                                                      const TOptClusterSize &  csize,
                                                      const idx_t              index_ofs,
                                                      data_t &                 data ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TBSPNDCTBuilder
//! \brief    Combines binary space partitioning with nested dissection
//!           based on connectivity defined by a sparse matrix.
//!
class TBSPNDCTBuilder : public TBSPCTBuilder
{
private:
    //! sparse matrix for connectivity between indices
    const TSparseMatrix *  _sparse_mat;

    // mode for handling interface cluster tree depth
    bool                   _sync_interface_depth;

public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct cluster tree with partition strategy defined by \a part_strat
    //! and connectivity defined by \a S
    TBSPNDCTBuilder ( const TSparseMatrix *  S,
                      const TBSPPartStrat *  part_strat,
                      const uint             n_min        = CFG::Cluster::nmin,
                      const uint             min_leaf_lvl = 0 );

    //! construct cluster tree with partition strategy defined by \a part_strat
    //! (must use "build( coord, S )")
    TBSPNDCTBuilder ( const TBSPPartStrat *  part_strat,
                      const uint             n_min        = CFG::Cluster::nmin,
                      const uint             min_leaf_lvl = 0 );

    virtual ~TBSPNDCTBuilder ();

    //////////////////////////////////////////////
    //
    // options
    //

    //! set mode for handling interface cluster tree depth, e.g. synchronise with domain clusters
    void  sync_interface_depth ( const bool  b )
    {
        _sync_interface_depth = b;
    }
    
    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //! build nested dissection cluster tree out of coordinate set \a coord and
    //! additional connectivity defined by local sparse matrix
    virtual std::unique_ptr< TClusterTree >  build     ( const TCoordinate *      coord,
                                                         const idx_t              idx_ofs = 0 ) const;

    //! build nested dissection cluster tree out of coordinate set \a coord and
    //! additional connectivity defined by sparse matrix \a S
    virtual std::unique_ptr< TClusterTree >  build     ( const TCoordinate *      coord,
                                                         const TSparseMatrix *    S,
                                                         const idx_t              idx_ofs = 0 ) const;
    
    //! recursively build cluster tree
    virtual std::unique_ptr< TGeomCluster >  divide    ( const TNodeSet &         dofs,
                                                         const uint               lvl,
                                                         const TBBox &            bbox,
                                                         const TOptClusterSize &  csize,
                                                         const idx_t              index_ofs,
                                                         data_t &                 data ) const;

    //! recursively build cluster tree for interfaces clusters
    virtual std::unique_ptr< TGeomCluster >  divide_if ( const TNodeSet &         dofs,
                                                         const uint               lvl,
                                                         const uint               max_lvl,
                                                         const TBBox &            bbox,
                                                         const TOptClusterSize &  csize,
                                                         const idx_t              index_ofs,
                                                         data_t &                 data ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TGeomPartCTBuilder
//! \brief    Enhances other geometrical ct builders by allowing the user to
//!           define the first level of index partitioning by a given vector.
//!
class TGeomPartCTBuilder : public TGeomCTBuilder
{
private:
    //! array defining first index partition:
    //! - _first_lvl_part[i] = j  ⇔  index i belongs to son j
    //! - partition ids have to start at 0 and must be consecutively numbered
    //! - root of cluster tree will have <max partition id> sons
    const std::vector< idx_t >  _first_lvl_part;
        
    //! controls adjustment of path lengths in subtrees, e.g.
    //! all to have similar depths in sub trees of all sons
    const bool                  _adjust_depth;
    
    //! base cluster tree builder, used for constructing all levels except 0
    const TGeomCTBuilder *      _base_builder;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct cluster tree with first partition defined by \a first_lvl_part
    //! and base ct builder \a base_builder
    TGeomPartCTBuilder ( const std::vector< idx_t > &  first_lvl_part,
                         const TGeomCTBuilder *        base_builder,
                         const adjust_depth_mode_t     adjust_depth );

    //! dtor
    virtual ~TGeomPartCTBuilder ();

    //! use user defined partition for first level and partitioning strategy for rest
    virtual std::unique_ptr< TGeomCluster >  divide ( const TNodeSet &         dofs,
                                                      const uint               lvl,
                                                      const TBBox &            bbox,
                                                      const TOptClusterSize &  csize,
                                                      const idx_t              index_ofs,
                                                      data_t &                 data ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TGeomGroupCTBuilder
//! \brief    Enhances geometrical ct builder by allowing to group indices, e.g.
//!           the groups are clustered and later expanded, ensuring
//!           that all indices in a group are in the same cluster
//!           NOTE: bounding boxes per index are not yet supported (bb_min/bb_max)
//!
class TGeomGroupCTBuilder : public TGeomCTBuilder
{
private:
    //! array defining index groups: _groups[i] = group of index i
    //! - group ids have to start at 0 and must be consecutively numbered
    const std::vector< uint > &  _groups;

    //! base cluster tree builder, used for constructing tree of grouped indices
    const TBSPCTBuilder *        _base_builder;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct cluster tree with index grouping defined by \a groups
    //! and base ct builder \a base_builder
    TGeomGroupCTBuilder ( const std::vector< uint > &  groups,
                          const TBSPCTBuilder *        base_builder );

    virtual ~TGeomGroupCTBuilder ();

    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //!
    //! build cluster tree out of coordinate set \coord but
    //! with respect to additional index groups defined by \see _groups
    //!
    virtual std::unique_ptr< TClusterTree >  build  ( const TCoordinate *      coord,
                                                      const idx_t              idx_ofs = 0 ) const;

    //! recursively build cluster tree for indices in \a dofs
    virtual std::unique_ptr< TGeomCluster >  divide ( const TNodeSet &         dofs,
                                                      const uint               lvl,
                                                      const TBBox &            bbox,
                                                      const TOptClusterSize &  csize,
                                                      const idx_t              index_ofs,
                                                      data_t &                 data ) const;
protected:

    //!
    //! recursively expand index groups and update the corresponding
    //! cluster, e.g. size, bounding box, etc.
    //!
    std::unique_ptr< TGeomCluster >  expand ( TCluster *                       cluster,
                                              const std::vector< std::vector< uint > > & group_to_idx,
                                              const TCoordinate *              coord,
                                              const TPermutation &             group_perm,
                                              TPermutation &                   perm,
                                              const idx_t                      idx_ofs ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TMBLRBSPCTBuilder
//! \brief    implement MBLR clustering
//!
class TMBLRCTBuilder : public TGeomCTBuilder
{
protected:
    //! type of partitioning strategy
    const TBSPPartStrat *  _part_strat;

    //! number of levels in MBLR clustering
    const size_t           _nlevel;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct BSP cluster tree builder with partitioning strategy \a part_strat
    TMBLRCTBuilder ( const size_t           nlevel,
                     const TBSPPartStrat *  part_strat,
                     const uint             n_min = CFG::Cluster::nmin );

    //! dtor
    virtual ~TMBLRCTBuilder ();

    //////////////////////////////////////////////
    //
    // clustertree creation
    //

    //! recursively build cluster tree for indices in \a dofs
    virtual std::unique_ptr< TGeomCluster >  divide ( const TNodeSet &         dofs,
                                                      const uint               lvl,
                                                      const TBBox &            bbox,
                                                      const TOptClusterSize &  csize,
                                                      const idx_t              index_ofs,
                                                      data_t &                 data ) const;
};

}// namespace HLIB

#endif  // __HLIB_TBSPCTBUILDER_HH
