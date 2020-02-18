#ifndef __HLIB_VEC_CONV_HH
#define __HLIB_VEC_CONV_HH
//
// Project     : HLib
// File        : vec_conv.hh
// Description : vector conversion functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/vector/TVector.hh"
#include "hpro/cluster/TClusterBasis.hh"

namespace HLIB
{

//!
//! convert \a v into scalar vector
//!
TVector *
to_scalar  ( const TVector *             v );

//!
//! convert \a v into blocked vector as defined by
//! to cluster tree \a ct
//!
TVector *
to_blocked ( const TVector *             v,
             const TCluster *            ct );

//!
//! convert \a v into (blocked) uniform vector with basis
//! given in \a cb
//!
template <typename T>
TUniformVector *
to_uniform ( const TVector *             v,
             const TClusterBasis< T > *  cb );

}// namespace

#endif // __HLIB_VEC_CONV_HH
