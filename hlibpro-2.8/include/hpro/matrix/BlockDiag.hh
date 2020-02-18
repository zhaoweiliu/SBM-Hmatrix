#ifndef __HLIB_BLOCKDIAG_HH
#define __HLIB_BLOCKDIAG_HH
//
// Project     : HLib
// File        : BlockDiag.hh
// Description : restrict H-matrices to blockdiagonal form
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"

namespace HLIB
{

// determines level of leafs
const uint  lvl_leaf = uint(-1);

//
// restrict given H-Matrix to blockdiagonal part of first
// <lvl> levels or to diagonal blocks of size at most <blocksize>
//
void
restrict_blockdiag ( TMatrix *          A,
                     const uint         lvl       = lvl_leaf,
                     const size_t       blocksize = 0 );

//
// copy blockdiagonal part (first <lvl> levels or blocks of size
// at most <blocksize>) of given H-Matrix
//
std::unique_ptr< TMatrix >
blockdiag          ( const TMatrix *    A,
                     const uint         lvl       = lvl_leaf,
                     const size_t       blocksize = 0 );

//
// same as above, but copy up to given accuracy
//
std::unique_ptr< TMatrix >
blockdiag          ( const TMatrix *    A,
                     const TTruncAcc &  acc,
                     const uint         lvl       = lvl_leaf,
                     const size_t       blocksize = 0 );

}// namespace

#endif  // __HLIB_BLOCKDIAG_HH
