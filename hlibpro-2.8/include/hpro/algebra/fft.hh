#ifndef __HLIB_FFT_HH
#define __HLIB_FFT_HH
//
// Project     : HLib
// File        : fft.hh
// Description : provides FFT for vectors
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

namespace HLIB
{

//!
//! \{
//! \name Fourier Transformation
//! Functions for the forward and backward Fourier transformation of vectors.
//!

//!
//! \ingroup Algebra_Module
//! \brief   perform FFT for vector \a v (inplace)
//!
void fft ( TVector * v );

//!
//! \ingroup Algebra_Module
//! \brief   perform inverse FFT for vector \a v (inplace)
//!
void ifft ( TVector * v );

//! \}

}// namespace

#endif  // __HLIB_FFT_HH
