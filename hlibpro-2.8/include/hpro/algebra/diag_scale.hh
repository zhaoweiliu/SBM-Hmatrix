#ifndef __HLIB_DIAG_SCALE_HH
#define __HLIB_DIAG_SCALE_HH
//
// Project     : HLib
// File        : diag_scale.hh
// Description : compute and apply diagonal scaling for matrices
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/vector/TScalarVector.hh"

namespace HLIB
{

/////////////////////////////////////////////////////////////////////////////
//
// compute diagonal matrices D1 and D2 such that
// D1·A·D2 has (roughly) normalised row and column norms
// - the computed scaling factors are stored in the given vectors
//
/////////////////////////////////////////////////////////////////////////////

void
compute_diag_scale ( const TMatrix *        A,
                     TScalarVector &        D1,
                     TScalarVector &        D2 );

/////////////////////////////////////////////////////////////////////////////
//
// apply given diagonal scaling to matrix A
//
/////////////////////////////////////////////////////////////////////////////

void
apply_diag_scale   ( TMatrix *              A,
                     const TScalarVector &  D1,
                     const TScalarVector &  D2 );

}// namespace

#endif  // __HLIB_MUL_DIAG_HH
