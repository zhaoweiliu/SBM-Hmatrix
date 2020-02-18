#ifndef __HLIB_TDISTRBC_HH
#define __HLIB_TDISTRBC_HH
//
// Project     : HLib
// File        : TDistrBC.hh
// Description : distribute given block-clustertree over processors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <map>
#include <list>
#include <vector>

#include "hpro/base/types.hh"

#include "hpro/cluster/TBlockCluster.hh"
#include "hpro/cluster/TBlockClusterTree.hh"

namespace HLIB
{

///////////////////////////////////////////////////////////////////////////
//!
//! \class TDistrBC
//! \brief Base class for all block cluster distribution methods.
//!
class TDistrBC
{
public:
    
    //!
    //! \class  TCostFunc
    //! \brief  Cost function for block clusters in load balancing.
    //!
    //!         The cost function has to return costs associated with
    //!         given cluster only, i.e. not for son clusters in case
    //!         of inner nodes.
    //!
    class TCostFunc
    {
    public:
        // ctor
        TCostFunc () {}

        // dtor
        virtual ~TCostFunc () {}

        //!
        //! return cost of cluster \a cl \b not including costs for son clusters
        //!
        virtual double  eval      ( const TBlockCluster * cl ) const = 0;

        //!
        //! return cost of cluster \a cl including costs for son clusters
        //!
        virtual double  eval_rec  ( const TBlockCluster * cl ) const
        {
            double  costs = eval( cl );

            for ( uint  i = 0; i < cl->nsons(); ++i )
            {
                if ( cl->son( i ) != NULL )
                    costs += eval_rec( cl->son( i ) );
            }// for

            return costs;
        }
    };

public:
    ///////////////////////////////////////////////
    //
    // ctors and dtor
    //

    TDistrBC () {}

    virtual ~TDistrBC () {}

    ///////////////////////////////////////////////
    //
    // distribution techniques
    //

    //! distribute block cluster tree
    /*!
      distribute given block cluster \a tree onto \a p processors with
      costs for block clusters provided by cost-function \a cf
    */
    virtual void distribute ( const uint                   p,
                              TBlockCluster *              tree,
                              const TDistrBC::TCostFunc *  cf ) const = 0;

    //! distribute block cluster tree (tree version)
    virtual void distribute ( const uint                   p,
                              TBlockClusterTree *          tree,
                              const TDistrBC::TCostFunc *  cf ) const
    {
        distribute( p, const_cast< TBlockCluster * >( tree->root() ), cf );
    }
};

///////////////////////////////////////////////////////////////////////////
//!
//! \class TBlockDistrBC
//! \brief class for block-wise block cluster tree distribution
//!
//!        Distribute block cluster trees by traversing tree using BFS and
//!        collect at least p nodes, where p is the number of processors,
//!        on a level. The collected nodes are then distributed onto the p
//!        processors. This ensures, that all nodes below that level are
//!        assigned to a unique processor. All nodes above belong to all
//!        processors. This is necessary for the factorisation algorithms.
//!
//!        The number of nodes to schedule can be further controlled by a
//!        minimal depth in the BFS, e.g. even if enough nodes were found
//!        earlier, the algorithm will proceed until level \c _min_lvl.
//!
//!        In case of symmetric/hermitian matrices, the scheduling can be
//!        performed for the lower left part of the block cluster tree.
//!        The upper right part will be assigned by mirroring the distribution
//!        pattern.
//!
//!        Make sure, that leaves in the tree will apear at most on the
//!        level of the scheduled nodes and not earlier.
//!
class TBlockDistrBC : public TDistrBC
{
private:
    // minimal level to traverse tree
    const uint  _min_lvl;
    
    // schedule only lower left half for symmetric matrices
    // - upper right is transposed
    const bool  _symmetric;
    
public:
    ///////////////////////////////////////////////
    //
    // ctors and dtor
    //

    TBlockDistrBC ( const uint  min_lvl = 0,
                    const bool  sym     = false )
            : _min_lvl( min_lvl ),
              _symmetric( sym )
    {}

    virtual ~TBlockDistrBC () {}

    ///////////////////////////////////////////////
    //
    // distribution techniques
    //

    //! distribute block cluster tree
    virtual void distribute ( const uint                   p,
                              TBlockCluster *              tree,
                              const TDistrBC::TCostFunc *  cf ) const;

    using TDistrBC::distribute;
};

///////////////////////////////////////////////////////////////////////////

class TBlockCyclicDistrBC : public TDistrBC
{
public:
    ///////////////////////////////////////////////
    //
    // ctors and dtor
    //

    TBlockCyclicDistrBC ();

    virtual ~TBlockCyclicDistrBC () {}

    ///////////////////////////////////////////////
    //
    // distribution techniques
    //

    //! distribute block cluster tree
    virtual void distribute ( const uint                   p,
                              TBlockCluster *              tree,
                              const TDistrBC::TCostFunc *  cf ) const;

    using TDistrBC::distribute;
};

///////////////////////////////////////////////////////////////////////////
//!
//! \class TSFCDistrBC
//! \brief Class for distributing block cluster trees using space
//!        filling curves.
//!
//!        The partitioning of the block index set as defined by the
//!        nodes of a given level of the block cluster tree is mapped
//!        to a one dimensional sequence using space filling
//!        curves. This sequence is afterwards partitioned into p
//!        intervals by a sequence partitioning scheduler.
//!
//!        This procedure is performed for each level of the block cluster
//!        from root to leaves, as long as no load balance was achieved.
//!        The difference between the individual loads in all intervals
//!        is defined by ε, e.g. [min load, max load] ⊆ [-ε,+ε]·avg load.
//!
class TSFCDistrBC : public TDistrBC
{
public:
    //! different types of space-filling-curves
    enum curve_t { Z_CURVE,
                   LEBESGUE_CURVE,
                   HILBERT_CURVE,
                   MOORE_CURVE };

private:

    //! adjust for coarsening
    const bool     _adjust_coarsen;
    
    //! type of space filling curve to use
    const curve_t  _curve_type;
    
    //! approximation ration for scheduling
    const double   _epsilon;

public:
    ///////////////////////////////////////////////
    //
    // ctors and dtor
    //

    //! construct block cluster distributor
    TSFCDistrBC ( const bool     adj_coarsen = false,
                  const curve_t  curve_type  = HILBERT_CURVE,
                  const double   eps         = 0.05 );

    virtual ~TSFCDistrBC () {}

    ///////////////////////////////////////////////
    //
    // distribution techniques
    //

    //! distribute \a tree, obtain costs from cost-function \a cf
    virtual void distribute ( const uint                   p,
                              TBlockCluster *              tree,
                              const TDistrBC::TCostFunc *  cf ) const;

    // for TBlockClusterTree version
    using TDistrBC::distribute;

protected:

    DISABLE_COPY_OP( TSFCDistrBC );
};

///////////////////////////////////////////////////////////////////////////
//!
//! \class TNDDistrBC
//! \brief class for block cluster tree distribution for nested dissection
//!
//!        The block cluster tree is distributed recursively starting with
//!        the root, which is assigned to all processors. On each level, 
//!        the processor set is partitioned according to the number of 
//!        domain clusters (diagonal blocks, except interface).
//!
//!        Offdiagonal blocks are assigned to the corresponding domain
//!        processor sets.
//!
//!        The cost function is not used in the current implementation
//!        (assuming, that the cluster partitioning results in equal load).
//!
class TNDDistrBC : public TDistrBC
{
public:
    ///////////////////////////////////////////////
    //
    // ctors and dtor
    //

    TNDDistrBC ();

    virtual ~TNDDistrBC () {}

    ///////////////////////////////////////////////
    //
    // distribution techniques
    //

    //! distribute block cluster tree
    virtual void distribute ( const uint                   p,
                              TBlockCluster *              tree,
                              const TDistrBC::TCostFunc *  cf ) const;

    using TDistrBC::distribute;
};

}// namespace HLIB

#endif  // __TDISTRBC_HH
