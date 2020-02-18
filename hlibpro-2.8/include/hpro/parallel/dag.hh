#ifndef __HLIB_DAG_HH
#define __HLIB_DAG_HH
//
// Project     : HLib
// File        : dag.hh
// Description : classes and functions for handling DAGs
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <list>
#include <deque>
#include <atomic>

#include <tbb/task.h>
#include <tbb/spin_mutex.h>

#include "hlib-config.h"

#include "hpro/base/error.hh"
#include "hpro/base/TTruncAcc.hh"
#include "hpro/misc/TProgressBar.hh"
#include "hpro/parallel/TMutex.hh"

namespace HLIB
{

//!
//! \class  TDAGNode
//! \brief  base class for DAG nodes defining interface and basic functionality
//!
class TDAGNode : public tbb::task
{
private:
    //! @cond
    
    // list of successors
    std::list< TDAGNode * >  _successors;

    #if HLIB_DEBUG == 1
    // DEBUG: list of predecessors
    std::list< TDAGNode * >  _predecessors;
    #endif
    
    //! @endcond
    
public:

    //! default ctor
    TDAGNode ()
    {}

    virtual  ~TDAGNode () {}

    //! overwritten TBB method executed by thread
    tbb::task *  execute ();

    //! do actual work associated with node; return true if work was finished
    virtual bool  compute  () = 0;
    
    //! add dependency, e.g. this node requires \a dep to finish
    void  add_dep ( TDAGNode *  dep )
    {
        if ( dep == NULL )
            HERROR( ERR_ARG, "(TDAGNode) add_dep", "node argument is NULL" );

        #if HLIB_DEBUG == 1
        // DEBUG
        _predecessors.push_back( dep );
        #endif
        
        // add this node as successor of \a dep
        dep->add_succ( this );
    }

    //! signals, if node has successors
    bool    has_succ      () const { return ! _successors.empty(); }
    
    //! give access to successor list
    const std::list< TDAGNode * > &  successors   () const { return _successors; }

    #if HLIB_DEBUG == 1
    //! give access to predecessor list
    const std::list< TDAGNode * > &  predecessors () const { return _predecessors; }
    #endif
    
    // return textual description of node
    virtual std::string  to_string () const = 0;
    
protected:

    //! spawn given node (for debugging)
    void  spawn_node ( TDAGNode * node );
    
    //! spawn ready successors
    void  spawn_succ ();
    
    //! add successor; also increasing dependency counter of \a succ
    void  add_succ ( TDAGNode *  succ )
    {
        if ( succ == NULL )
            HERROR( ERR_ARG, "(TDAGNode) add_succ", "node argument is NULL" );
        
        _successors.push_back( succ );
        succ->increment_ref_count();
    }

    // node allocation
    template < typename T,
               typename ... Args >
    T *
    alloc_child ( Args && ... args )
    {
        return new( allocate_child() ) T( std::forward< Args >( args ) ... );
    }

};

//!
//! \brief  Execute nodes in DAG reachable from nodes in \a start.
//!
void
run_dag ( std::list< TDAGNode * > &  start,
          tbb::task *                final );

//!
//! \brief  Test reachability in DAG
//!
void
test_dag ( std::list< TDAGNode * > &   start,
           tbb::task *                 final,
           std::deque< TDAGNode * > &  all_nodes );

//!
//! \brief  print DAG nodes (accessible from start nodes only)
//!
void
print_dag ( std::list< TDAGNode * > &  start_nodes );

//!
//! \brief  print statistical data about DAG
//!
void
print_dag_stat ( std::list< TDAGNode * > &  start_nodes );

//!
//! \brief  print given set of nodes
//!
void
print_nodes ( std::deque< TDAGNode * > &  nodes );

//!
//! \brief  print given set of nodes in DOT format
//!
void
print_nodes ( std::deque< TDAGNode * > &  nodes,
              const std::string &         filename );

//
// node allocation
//
template < typename T,
           typename ... Args >
T *
alloc_root ( Args && ... args )
{
    return new( tbb::task::allocate_root() ) T( std::forward< Args >( args ) ... );
}

//////////////////////////////////////////////////////////////////////////////////////////////
//
// alternative DAG interface
//
//////////////////////////////////////////////////////////////////////////////////////////////

namespace DAG
{

//
// forward declarations
//
class  Node;
class  RuntimeTask;

//
// defines a memory block by an Id and a block index set
// (memory block as used within a matrix)
//
struct mem_block_t
{
    id_t            id;
    TBlockIndexSet  is;
};

// set of memory blocks
using  block_list_t = std::vector< mem_block_t >;

// list/vector of nodes
using  node_vec_t   = std::vector< Node * >;
using  node_list_t  = std::list< Node * >;

//!
//! class for local sub graph during refinement
//!
class LocalGraph : public node_vec_t
{
protected:
    // signals finished graph, e.g. no update of edges needed
    bool  _finished;
    
public:
    // ctor
    LocalGraph ()
            : _finished( false )
    {}
    
    // add node and apply dependencies based on existing nodes
    void add_node_and_dependencies ( Node *  node );

    // only add node to graph
    void add_node                  ( Node *  node )
    {
        push_back( node );
    }

    // set dependencies between all nodes in graph based on
    // in/out data blocks of nodes
    void  set_dependencies ();

    // return finish status of graph
    bool  is_finalized () const { return _finished; }
    
    // signal finished graph
    void  finalize     ()       { _finished = true; }
    
    //
    // wrapper to simultaneously allocate node and put into list local graph
    //
    template < typename T,
               typename ... Args >
    T *
    alloc_node ( Args && ...    args )
    {
        auto  node = new T( std::forward< Args >( args ) ... );

        push_back( node );

        return node;
    }
};

//!
//! @class Node
//!
//! @brief Represents a node in a DAG
//!
//! A node in a DAG with incoming and outgoing edges (dependencies)
//! - data dependencies are described in the form of mem_block_t lists (in/out)
//! - also implements actual actions to perform per node
//! - user defined functions are "run_", "in_blocks_" and "out_blocks_"
//!
class Node
{
private:
    // successor nodes in DAG
    node_vec_t             _successors;

    // number of dependencies (incoming edges)
    int                    _ndeps;

    // dependency counter (#incoming edges)
    std::atomic< int >     _dep_cnt;

    // block index set dependencies for automatic dependency refinement
    block_list_t           _in_blk_deps;
    block_list_t           _out_blk_deps;

    // set of sub nodes
    node_vec_t             _sub_nodes;

    // mutex to handle concurrent access to internal data
    std::mutex             _mutex;
    
public:
    // ctor
    Node ();

    // dtor
    virtual ~Node () {}

    // per node initialization (e.g., data dependencies)
    void  init  ();
    
    // handles execution of node code and spawning of successor nodes
    void  run   ( const TTruncAcc & acc );

    // givess access to successor nodes
    node_vec_t &        successors ()       { return _successors; }
    const node_vec_t &  successors () const { return _successors; }

    //
    // task dependencies
    //
    
    // run <this> before <t>
    void  before   ( Node *  t )
    {
        if ( t == nullptr )
            HERROR( ERR_ARG, "Node::before", "node is null" );
        
        _successors.push_back( t );
    }
    
    // run <this> after <t>
    void  after    ( Node *  t )
    {
        if ( t == nullptr )
            HERROR( ERR_ARG, "Node::after", "node is null" );
        
        t->_successors.push_back( this );
    }

    // return dependency counter
    int   dep_cnt      () const   { return _dep_cnt; }
    
    // increase task dependency counter
    int   inc_dep_cnt  ()         { return ++_dep_cnt; }
    
    // decrease task dependency counter
    int   dec_dep_cnt  ()         { return --_dep_cnt; }
    
    // set task dependency counter
    void  set_dep_cnt  ( int  d ) { _dep_cnt = d; }
    
    // reset dependency counter
    void  reset_dep_cnt ()        { _dep_cnt = _ndeps; }
    
    //
    // data dependencies
    //
    
    // return local list of block index sets for input dependencies
    const block_list_t &  in_blocks  () const { return _in_blk_deps; }
    
    // return local list of block index sets for output dependencies
    const block_list_t &  out_blocks () const { return _out_blk_deps; }

    //
    // refinement
    //
    
    // return true if node is refined
    bool  is_refined  () const { return ! _sub_nodes.empty(); }

    // give access to sub nodes
    auto  sub_nodes ()       -> decltype(_sub_nodes) { return _sub_nodes; };
    auto  sub_nodes () const -> decltype(_sub_nodes) { return _sub_nodes; };
    
    // split node into subnodes and update dependencies
    // if retval is empty, no refinement was done
    void  refine ();

    // refine dependencies of <node>, either local or of sub nodes
    // - return true, if <node> changed (modified dependencies)
    // - if <do_lock> == true, lock affected nodes during refinement
    bool  refine_deps  ( const bool  do_lock );

    // finalize node data (if internal data will not change)
    void  finalize ();

    //
    // mutex functions
    //

    void lock      () { _mutex.lock(); }
    void unlock    () { _mutex.unlock(); }
    bool try_lock  () { return _mutex.try_lock(); }
    
    //
    // output and visualization
    //
    
    // print node with full edge information
    void  print () const;
    
    // return text version of node
    virtual std::string  to_string () const { return "Node"; }

    // (optional) color for DAG visualization (format: RRGGBB)
    virtual std::string  color     () const { return "FFFFFF"; }

private:

    //
    // private functions used by above public wrappers
    //
    
    virtual void  run_ ( const TTruncAcc & )
    {
        HERROR( ERR_NOT_IMPL, "", "" );
    }

    virtual const block_list_t  in_blocks_ () const
    {
        return block_list_t();
    }
    
    virtual const block_list_t  out_blocks_ () const
    {
        return block_list_t();
    }

    virtual LocalGraph          refine_ ()
    {
        HERROR( ERR_NOT_IMPL, "", "" );
    }
};

//
// wrapper to simultaneously allocate node and put into list of subnodes
//
template < typename T,
           typename ... Args >
T *
alloc_node ( LocalGraph &  g, Args && ... args )
{
    auto  node = new T( std::forward< Args >( args ) ... );

    g.add_node_and_dependencies( node );

    return node;
}

template < typename T,
           typename ... Args >
T *
alloc_node_wo_dep ( LocalGraph &  g, Args && ... args )
{
    auto  node = new T( std::forward< Args >( args ) ... );

    g.add_node( node );

    return node;
}

template < typename T,
           typename ... Args >
T *
alloc_node ( std::list< Node * > &  subnodes,
             Args && ...            args )
{
    auto  node = new T( std::forward< Args >( args ) ... );

    subnodes.push_back( node );

    return node;
}

//
// directed acyclic graph (DAG)
// - only holds list of nodes, start and end nodes
//
class Graph
{
private:
    node_list_t  _nodes;
    node_list_t  _start;
    node_list_t  _end;

public:
    // ctor
    Graph ( node_list_t &  nodes,
            node_list_t &  start,
            node_list_t &  end );

    Graph ( node_list_t &&  nodes,
            node_list_t &&  start,
            node_list_t &&  end );

    Graph ( Graph &&        g );

    // return number of nodes
    size_t  nnodes () const { return _nodes.size(); }

    // return number of (out) edges
    size_t  nedges () const;

    // return list of all/start/end nodes
    node_list_t &        nodes ()       { return _nodes; }
    const node_list_t &  nodes () const { return _nodes; }
    node_list_t &        start ()       { return _start; }
    const node_list_t &  start () const { return _start; }
    node_list_t &        end   ()       { return _end; }
    const node_list_t &  end   () const { return _end; }

    // add given set of nodes to DAG
    // (assumption: dependency counter already set)
    void    add_nodes ( node_list_t &  nodes );
    
    // execute given DAG
    void    run    ( const TTruncAcc &  acc,
                     TProgressBar *     progress = nullptr );
    
    // simulate execution of DAG and look if all nodes are handled and
    // all end nodes are reached
    void    test     ();

    // output DAG
    void    print () const;
    
    // output DAG in DOT format
    void    print_dot ( const std::string &  filename ) const;
};

//
// construct DAG based on refinement of given node
//
Graph
refine   ( Node *   node );

//
// concatenate dag1 and dag2
// - end nodes of dag1 are dependencies for start nodes of dag2
//
Graph
concat   ( Graph &  dag1,
           Graph &  dag2 );

//
// merge dag1 and dag2
//
Graph
merge    ( Graph &  dag1,
           Graph &  dag2 );

}// namespace DAG

}// namespace HLIB

#endif  // __HLIB_DAG_HH
