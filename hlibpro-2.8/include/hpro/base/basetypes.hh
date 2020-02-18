#ifndef __HLIB_BASETYPES_HH
#define __HLIB_BASETYPES_HH
//
// Project     : HLib
// File        : basetypes.hh
// Description : basic type definitions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <cstdlib>

#include "hlib-config.h"

namespace HLIB
{

//
// abbr. for standard types
//
using  uchar = unsigned char;
using  uint  = unsigned int;
using  ulong = unsigned long;

//!
//! \typedef  real
//! \brief    default real valued type for all computations
//!
#if HLIB_SINGLE_PREC == 1
using  real = float;
#else
using  real = double;
#endif

//!
//! \typedef  idx_t
//! \brief    type for indices
//!
using  idx_t = long;

//!
//! \typedef  id_t
//! \brief    type for identifiers
//!
using  id_t  = ulong;

//!
//! type for type IDs
//!
using  typeid_t = uint;

}// namespace

#endif  // __HLIB_BASETYPES_HH
