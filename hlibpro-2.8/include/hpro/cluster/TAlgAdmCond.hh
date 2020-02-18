#ifndef __HLIB_TALGADMCOND_HH
#define __HLIB_TALGADMCOND_HH
//
// Project     : HLib
// File        : TAlgAdmCond.hh
// Description : algebraic admissibility condition for sparse matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <vector>

#include "hpro/cluster/TPermutation.hh"
#include "hpro/cluster/TNodeSet.hh"
#include "hpro/cluster/TAdmCondition.hh"

#include "hpro/matrix/TSparseMatrix.hh"

namespace HLIB
{

//!
//! \ingroup  Cluster_Module
//! \class    TAlgAdmCond
//! \brief    base class for algebraic admissibility conditions
//!
class TAlgAdmCond : public TAdmCondition
{
protected:
    //! @cond
    
    // sparse matrix defining the matrix graph
    const TSparseMatrix *  _mat;

    // mapping of index-names from external (in sparse matrix)
    // to internal numbering (in cluster tree)
    const TPermutation *   _row_perm_i2e, * _col_perm_i2e;
    TPermutation *         _row_perm_e2i, * _col_perm_e2i;
    
    //! @endcond
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor with graph defined by \a S and mapping of internal to external
    //! indices defined by \a perm_i2e (row and column mappings identical)
    TAlgAdmCond ( const TSparseMatrix *  S,
                  const TPermutation *   perm_i2e );

    //! ctor with graph defined by \a S and mapping of internal to external
    //! indices defined by \a row_perm_i2e and \a col_perm_i2e
    TAlgAdmCond ( const TSparseMatrix *  S,
                  const TPermutation *   row_perm_i2e,
                  const TPermutation *   col_perm_i2e );

    //! dtor
    virtual ~TAlgAdmCond ();

protected:

    DISABLE_COPY_OP( TAlgAdmCond );
};

//!
//! \ingroup  Cluster_Module
//! \class    TStdAlgAdmCond
//! \brief    Standard admissibility condition based on matrix graph criteria.
//!
class TStdAlgAdmCond : public TAlgAdmCond
{
protected:
    //! @cond
    
    // admissibility parameter
    const real                   _eta;
    
    // mark visited nodes (not mt-safe !!!)
    mutable std::vector< bool >  _visited;
    
    //! @endcond
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! ctor
    TStdAlgAdmCond ( const real             eta,
                     const TSparseMatrix *  S,
                     const TPermutation *   perm_i2e );

    //! ctor
    TStdAlgAdmCond ( const real             eta,
                     const TSparseMatrix *  S,
                     const TPermutation *   row_perm_i2e,
                     const TPermutation *   col_perm_i2e );

    //! dtor
    virtual ~TStdAlgAdmCond ();

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is admissible
    virtual bool  is_adm    ( const TBlockCluster * cl ) const;

protected:
    //! determine diameter of cluster \a cl
    virtual uint  diameter  ( const TCluster *      cl,
                              const TPermutation *  perm_i2e,
                              const TPermutation *  perm_e2i ) const;

    //! Perform a BFS from set \a start in matrix and store last visited nodes
    //! in \a last. Stop BFS if all nodes in \a tau have been visited. Return the depth
    //! of the BFS iteration.
    virtual uint  bfs       ( TNodeSet &            start,
                              TNodeSet &            last,
                              const TCluster *      tau,
                              const TPermutation *  perm_i2e,
                              const TPermutation *  perm_e2i ) const;

    //! return true, if distance between \a tau and \a sigma is bigger than \a min_dist
    virtual bool  cmp_dist  ( const TCluster *      tau,
                              const TCluster *      sigma,
                              const uint            min_dist ) const;

    //! return true if \a node is local to cluster tree \a cl
    bool          is_local  ( const TCluster *      cl,
                              const node_t          node,
                              const TPermutation *  perm_e2i ) const
    {
        idx_t  idx;
        
        if ( perm_e2i != NULL ) idx = perm_e2i->permute( node );
        else                    idx = node;
        
        return (( idx >= cl->first() ) && ( idx <= cl->last() ));
    }

    DISABLE_COPY_OP( TStdAlgAdmCond );
};

//!
//! \ingroup  Cluster_Module
//! \class    TStdAlgAdmCond
//! \brief    Weak admissibility condition based on matrix graph criteria.
//!
class TWeakAlgAdmCond : public TAlgAdmCond
{
private:
    // distance (in no. of edges) to test for between clusters
    const uint  _distance;
    
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct object for algebraic weak admissibility based
    //! on connectivity in \a S with optional permutation \a perm_i2e
    //! (from internal to external ordering; identical for rows and columns)
    TWeakAlgAdmCond ( const TSparseMatrix *  S,
                      const TPermutation *   perm_i2e,
                      const uint             distance = 1 );

    //! construct object for algebraic weak admissibility based
    //! on connectivity in \a S with optional row and column permutations
    //! \a row_perm_i2e and \a col_perm_i2e (from internal to external
    //! ordering)
    TWeakAlgAdmCond ( const TSparseMatrix *  S,
                      const TPermutation *   row_perm_i2e,
                      const TPermutation *   col_perm_i2e,
                      const uint             distance = 1 );

    //! dtor
    virtual ~TWeakAlgAdmCond ();

    ///////////////////////////////////////////
    //
    // check block-cluster if admissible
    //

    //! return true if \a cl is weakly admissible
    virtual bool is_adm ( const TBlockCluster *  c ) const;
};

}// namespace

#endif  // __HLIB_TALGADMCOND_HH
