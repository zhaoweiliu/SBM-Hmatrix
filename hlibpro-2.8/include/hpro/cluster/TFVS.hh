#ifndef __HLIB_TFVS_HH
#define __HLIB_TFVS_HH
//
// Project     : HLib
// File        : TFVS.hh
// Description : class for computing a feedback-vertex-set of a matrix graph
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <vector>
#include <list>
#include <string>

#include "hpro/matrix/TSparseMatrix.hh"

namespace HLIB
{

//!
//! \ingroup  Cluster_Module
//! \class    TFVS
//! \brief    Uses a heuristic algorithm to compute feedback vertex set of a
//!           directed graph represented by a sparse matrix
//!
class TFVS
{
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    TFVS ()
    {}
    
    //////////////////////////////////////////////
    //
    // compute FVS
    //

    //! build FVS for \a S
    void build_fvs ( const TSparseMatrix *        S,
                     std::list< node_t > &        fvs ) const;

    //! build FVS for \a S with nodes in \a del_nodes already deleted
    void build_fvs ( const TSparseMatrix *        S,
                     const std::list< node_t > &  del_nodes,
                     std::list< node_t > &        fvs ) const;

protected:
    //
    // compute FVS by transforming graph
    //
    void build_fvs ( std::vector< std::list< node_t > > &  in_edges,
                     std::vector< std::list< node_t > > &  out_edges,
                     std::vector< bool > &                 node_del,
                     std::list< node_t > &                 fvs ) const;
    
    //
    // various graph transformations
    //
    
    bool trans1 ( const node_t                          node,
                  std::vector< std::list< node_t > > &  in_edges,
                  std::vector< std::list< node_t > > &  out_edges,
                  std::vector< bool > &                 node_del,
                  std::list< node_t > &                 fvs,
                  bool                                  check = true ) const;
    bool trans2 ( const node_t                          node,
                  std::vector< std::list< node_t > > &  in_edges,
                  std::vector< std::list< node_t > > &  out_edges,
                  std::vector< bool > &                 node_del,
                  std::list< node_t > &                 fvs ) const;
    bool trans3 ( const node_t                          node,
                  std::vector< std::list< node_t > > &  in_edges,
                  std::vector< std::list< node_t > > &  out_edges,
                  std::vector< bool > &                 node_del,
                  std::list< node_t > &                 fvs ) const;
    bool trans4 ( const node_t                          node,
                  std::vector< std::list< node_t > > &  in_edges,
                  std::vector< std::list< node_t > > &  out_edges,
                  std::vector< bool > &                 node_del,
                  std::list< node_t > &                 fvs ) const;
    bool trans5 ( const node_t                          node,
                  std::vector< std::list< node_t > > &  in_edges,
                  std::vector< std::list< node_t > > &  out_edges,
                  std::vector< bool > &                 node_del,
                  std::list< node_t > &                 fvs ) const;

    //
    // add node to FVS
    //
    
    void add_to_FVS ( const node_t                          node,
                      std::vector< std::list< node_t > > &  in_edges,
                      std::vector< std::list< node_t > > &  out_edges,
                      std::vector< bool > &                 node_del,
                      std::list< node_t > &                 fvs ) const;

    //
    // print current graph (for debugging)
    //
    void print ( const std::string &                   filename,
                 std::vector< std::list< node_t > > &  in_edges,
                 std::vector< std::list< node_t > > &  out_edges,
                 std::vector< bool > &                 node_del ) const;
};
    
}// namespace

#endif // __HLIB_TFVS_HH
