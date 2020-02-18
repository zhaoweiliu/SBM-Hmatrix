#ifndef __HLIB_TCLUSTERBASISBUILDER_HH
#define __HLIB_TCLUSTERBASISBUILDER_HH
//
// Project     : HLib
// File        : TClusterBasisBuilder.hh
// Description : classes for constructing cluster bases
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/cluster/TClusterBasis.hh"
#include "hpro/cluster/TBlockCluster.hh"

#include "hpro/base/TTruncAcc.hh"

#include "hpro/blas/Matrix.hh"

#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

//!
//! \ingroup  Cluster_Module
//! \class    TDenseClusterBasisBuilder
//! \brief    class for constructing cluster bases using dense matrices
//!
template < typename T >
class TDenseClusterBasisBuilder
{
public:
    ///////////////////////////////////////////////////////
    //
    // build cluster bases
    //

    //! Build row and column cluster basis for block cluster tree
    //! \a bct and dense matrix \a M. The accuracy of the cluster
    //! basis is determined by \a acc
    std::pair< std::unique_ptr< TClusterBasis< T > >,
               std::unique_ptr< TClusterBasis< T > > >
    build ( const TBlockCluster *      bct,
            const BLAS::Matrix< T > &  M,
            const TTruncAcc &          acc ) const;
};

//!
//! \ingroup  Cluster_Module
//! \class    THClusterBasisBuilder
//! \brief    class for constructing cluster bases using H-matrices
//!
template < typename T >
class THClusterBasisBuilder
{
public:
    ///////////////////////////////////////////////////////
    //
    // build cluster bases
    //

    //! Build row and column cluster basis for H matrix \a M.
    //! The accuracy of the cluster basis is determined by \a acc
    std::pair< std::unique_ptr< TClusterBasis< T > >,
               std::unique_ptr< TClusterBasis< T > > >
    build ( const TMatrix *        M,
            const TTruncAcc &      acc ) const;
};

}// namespace HLIB

#endif  // __HLIB_TCLUSTERBASISBUILDER_HH
