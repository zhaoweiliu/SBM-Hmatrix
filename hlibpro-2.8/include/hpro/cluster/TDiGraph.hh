#ifndef __HLIB_TDIGRAPH_HH
#define __HLIB_TDIGRAPH_HH
//
// Project     : HLib
// File        : TDiGraph.hh
// Description : classes for representing directed graphs
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <ostream>

#include "hpro/cluster/TGraph.hh"

namespace HLIB
{

//!
//! \ingroup  Cluster_Module
//! \class    TDiGraph
//! \brief    Class for directed graph represented by adjacency matrix
//!           in sparse format (assuming sparse graph!)
//!
class TDiGraph
{
public:
    //////////////////////////////////////////////
    //!
    //! \class TIterator
    //! \brief iterator to predecessor/successor lists
    //!
    class TIterator
    {
    private:
        //! adjacency list to reference
        const std::vector< idx_t > *  _adj_list;
        
        //! lower bound for adjacency list
        size_t                        _lb;

        //! upper bound for adjacency list
        size_t                        _ub;

        //! current pointer to adjacent node
        size_t                        _idx;

    private:
        // forbid default constructor
        TIterator ();
        
    public:
        //
        // ctors and dtor
        //

        TIterator ( const std::vector< idx_t > & adj_list,
                    const size_t                 lb,
                    const size_t                 ub )
                : _adj_list( & adj_list ), _lb( lb ), _ub( ub ), _idx( lb )
        {}

        TIterator ( const TIterator & iter )
        {
            *this = iter;
        }

        // get status (end-of-list)
        bool        eol  () const { return (_idx >= _ub); }
        
        // increment
        TIterator & next () { return ++*this; }
        
        // convert to current node
        operator idx_t () const { return (*_adj_list)[ _idx ]; }

        // access nodes
        idx_t  node        () const { return (*_adj_list)[ _idx ]; }
        idx_t  operator () () const { return (*_adj_list)[ _idx ]; }
        
        //
        // usual operators
        //

        // copy
        TIterator & operator = ( const TIterator & i )
        {
            _adj_list = i._adj_list;
            _lb       = i._lb;
            _ub       = i._ub;
            _idx      = i._idx;
            return *this;
        }

        // iterate (prefix/postfix)
        TIterator & operator ++ ()    { _idx++; return *this; }
        TIterator   operator ++ (int) { TIterator tmp(*this); ++*this; return tmp; }

        // compare
        bool        operator == ( const TIterator & i ) const
        {
            return ( ( _adj_list == i._adj_list ) &&
                     ( _lb == i._lb) && ( _ub == i._ub ) &&
                     ( _idx == i._idx ) );
        }
                
        bool        operator != ( const TIterator & i ) const
        {
            return ! (*this == i);
        }
    };
        
private:
    //! contains pointers into predecessor list for each node
    std::vector< idx_t >   _pred_list_ptr;
    
    //! contains list of predecessor nodes for each node
    std::vector< node_t >  _pred_nodes;

    //! contains pointers into successor list for each node
    std::vector< idx_t >   _succ_list_ptr;
    
    //! contains list of successor nodes for each node
    std::vector< node_t >  _succ_nodes;

public:
    //////////////////////////////////////////////
    //
    // constructor
    //

    //! construct empty digraph
    TDiGraph () {}

    //! construct digraph based on pattern in S using coefficients a_ij
    //! with |a_ij| > ε
    TDiGraph ( const TSparseMatrix * S,
               const real            eps = 0.0 );

    virtual ~TDiGraph () {}
    
    //! return number of nodes in graph
    size_t n_nodes () const
    {
        return (_pred_list_ptr.size() > 0 ? _pred_list_ptr.size() - 1 : 0);
    }
    
    //! return number of edges in graph
    size_t n_edges () const
    {
        return _pred_nodes.size() + _succ_nodes.size();
    }

    //////////////////////////////////////////////
    //
    // access graph data
    //

    //! return iterator to predecessors of node
    TIterator  predecessors ( const node_t  node ) const
    {
        return TIterator( _pred_nodes, _pred_list_ptr[ node ], _pred_list_ptr[ node+1 ] );
    }

    //! return iterator to successors of node
    TIterator  successors   ( const node_t  node ) const
    {
        return TIterator( _succ_nodes, _succ_list_ptr[ node ], _succ_list_ptr[ node+1 ] );
    }

    //////////////////////////////////////////////
    //
    // output
    //

    //! print graph
    void print ( std::ostream & os ) const;

    //! print graph with labels
    void print ( std::ostream &               os,
                 const std::vector< uint > &  label ) const;

    //! write in Chaco/Jostle/Metis file format
    void write ( std::ostream & os ) const;
    
    //////////////////////////////////////////////
    //
    // misc.
    //

    //! copy operator
    TDiGraph & operator = ( const TDiGraph & graph );
    
protected:
    //////////////////////////////////////////////
    //
    // internal methods
    //

    //! initialise graph for n nodes, m1 pred. edges and m2 succ. edges
    virtual void init ( const size_t  nnodes,
                        const size_t  npred,
                        const size_t  nsucc )
    {
        _pred_nodes.resize( npred );
        _succ_nodes.resize( nsucc );
        _pred_list_ptr.resize( nnodes + 1 );
        _succ_list_ptr.resize( nnodes + 1 );
    }
};

}// namespace HLIB

#endif  // __HLIB_TDIGRAPH_HH
