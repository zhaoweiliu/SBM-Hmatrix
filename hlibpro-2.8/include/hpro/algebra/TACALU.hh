#ifndef __HLIB_TACALU_HH
#define __HLIB_TACALU_HH
//
// Project     : HLib
// File        : TACALU.hh
// Description : LU factorisation for H-matrices based on ACA
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/misc/TProgressBar.hh"
#include "hpro/algebra/mat_fac.hh"

namespace HLIB
{

//! @cond

//!
//! \ingroup Algebra_Module
//! \class  TACALU
//! \brief  Computes LU factorisation LU = A by means of ACA
//!
class TACALU
{
private:
    //! local options for factorisation
    fac_options_t  _options;
    
public:
    //////////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! standard constructor with optional LU settings
    TACALU  ( const fac_options_t  opts = fac_options_t() )
            : _options( opts )
    {}
    
    //! destructor
    ~TACALU () {}

    //! set options for factorisation
    void set_options ( const fac_options_t  opts )
    {
        _options = opts;
    }
    
    //
    // compute LU factorisation of given matrix
    // - A is overwritten with L,U
    //

    //! compute LU factorisation of given matrix
    //! \param  nthreads  number of threads to use
    //! \param  A         on input matrix to factorise; on output factors L and U 
    //! \param  acc       accuracy of factorisation
    void factorise ( const uint         nthreads,
                     TMatrix *          A,
                     const TTruncAcc &  acc ) const;

    //
    // misc.
    //

    //! return suitable matrix representation for evaluating factors L and U
    //! \param  A        LU factors to be represented
    TMatrix *  eval_matrix  ( TMatrix *  A ) const;
    
    //! return suitable inverse representation of factors L and U
    //! \param  A        Cholesky factor to be represented
    TMatrix *  inv_matrix  ( TMatrix *   A ) const;
    
    //! return number of factorisation steps for A for progress meter 
    size_t     pm_steps    ( const TMatrix * A ) const;

    DISABLE_COPY_OP( TACALU );
};

//! @endcond

}// namespace HLIB

#endif  // __HLIB_TACALU_HH
