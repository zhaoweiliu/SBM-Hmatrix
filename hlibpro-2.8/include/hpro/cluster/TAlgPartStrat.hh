#ifndef __HLIB_TALGPARTSTRAT_HH
#define __HLIB_TALGPARTSTRAT_HH
//
// Project     : HLib
// File        : TAlgPartStrat.hh
// Description : partitioning strategies for algebraic clustering
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/cluster/TGraph.hh"
#include "hpro/cluster/TNodeSet.hh"

namespace HLIB
{

//!
//! \ingroup  Cluster_Module
//! \class    TAlgPartStrat
//! \brief    Base class for partitioning strategies for algebraic clustering.
//!
class TAlgPartStrat
{
public:

    // dtor
    virtual ~TAlgPartStrat () {}
    
    //!
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    //!
    virtual void partition ( const TGraph &  graph,
                             TNodeSet &      left,
                             TNodeSet &      right ) const = 0;
};

//!
//! \ingroup  Cluster_Module
//! \class    TBFSAlgPartStrat
//! \brief    Graph partitioning using BFS algorithm and FM optimisation.
//!
class TBFSAlgPartStrat : public TAlgPartStrat
{
public:
    // dtor
    virtual ~TBFSAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph &  graph,
                             TNodeSet &      left,
                             TNodeSet &      right ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TMLAlgPartStrat
//! \brief    Multi level graph partitioning.
//!
class TMLAlgPartStrat : public TAlgPartStrat
{
public:
    // dtor
    virtual ~TMLAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph &  graph,
                             TNodeSet &      left,
                             TNodeSet &      right ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TMETISAlgPartStrat
//! \brief    Graph partitioning using METIS
//!
class TMETISAlgPartStrat : public TAlgPartStrat
{
public:

    // ctor
    TMETISAlgPartStrat ( const bool use_random = CFG::Cluster::METIS_random );
    
    // dtor
    virtual ~TMETISAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph &  graph,
                             TNodeSet &      left,
                             TNodeSet &      right ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TScotchAlgPartStrat
//! \brief    Graph partitioning using Scotch
//!
class TScotchAlgPartStrat : public TAlgPartStrat
{
public:
    // dtor
    virtual ~TScotchAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph & graph,
                             TNodeSet &     left,
                             TNodeSet &     right ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    TChacoAlgPartStrat
//! \brief    Graph partitioning using CHACO
//!
class TChacoAlgPartStrat : public TAlgPartStrat
{
public:
    // dtor
    virtual ~TChacoAlgPartStrat () {}
    
    //! compute graph bi-partitioning of \a graph and store result in \a left and \a right
    virtual void partition ( const TGraph & graph,
                             TNodeSet &     left,
                             TNodeSet &     right ) const;
};

}// namespace HLIB

#endif  // __HLIB_TALGPARTSTRAT_HH
