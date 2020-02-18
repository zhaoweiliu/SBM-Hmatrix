#ifndef __HLIB_MAT_NORM_HH
#define __HLIB_MAT_NORM_HH
//
// Project     : HLib
// File        : mat_norm.hh
// Description : provides methods for matrix-norms
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <limits>

#include "hpro/matrix/TMatrix.hh"
#include "hpro/matrix/TFacInvMatrix.hh"

namespace HLIB
{

/////////////////////////////////////////////////////////////////////////
//!
//! \class  TMatrixNorm
//! \brief  baseclass for matrix norm computations
//!
class TMatrixNorm
{
public:
    //////////////////////////////////////////
    //
    // constructor and destructor
    //

    TMatrixNorm ();

    virtual ~TMatrixNorm () {}

    //////////////////////////////////////////
    //
    // matrix-norms
    //

    //! compute the norm of A
    virtual real norm       ( const TMatrix *  A ) const = 0;

    //! compute difference norm ‖A-B‖ or ‖A-B‖/‖A‖ (if \a rel == true)
    virtual real diff_norm  ( const TMatrix *  A,
                              const TMatrix *  B,
                              const bool       rel = true ) const = 0;
};

/////////////////////////////////////////////////////////////////////////
//!
//! \class  TFrobeniusNorm
//! \brief  computes Frobenius norm ‖·‖_F of a matrix
//!
class TFrobeniusNorm : public TMatrixNorm
{
public:
    //////////////////////////////////////////
    //
    // constructor and destructor
    //

    TFrobeniusNorm ();

    virtual ~TFrobeniusNorm () {}

    //////////////////////////////////////////
    //
    // matrix-norms
    //

    //! compute Frobenius norm of A
    virtual real norm       ( const TMatrix *  A ) const;

    //! compute difference norm ‖A-B‖_F or ‖A-B‖_F/‖A‖_F (if \a rel == true)
    virtual real diff_norm  ( const TMatrix *  A,
                              const TMatrix *  B,
                              const bool       rel = true ) const;
};

//!
//! \fn     norm_F
//! \brief  compute Frobenius norm of A (functional form of \see TFrobeniusNorm)
//!
inline real
norm_F ( const TMatrix *  A )
{
    TFrobeniusNorm  mnorm;

    return mnorm.norm( A );
}

//!
//! \fn     diff_norm_F
//! \brief  compute Frobenius norm of \a A - \a B (functional form of \see TFrobeniusNorm)
//!
inline real
diff_norm_F ( const TMatrix *  A,
              const TMatrix *  B,
              const bool       relative = true )
{
    TFrobeniusNorm  mnorm;

    return mnorm.diff_norm( A, B, relative );
}
inline real
diffnorm_F ( const TMatrix *  A,
             const TMatrix *  B,
             const bool       relative = true )
{
    return diff_norm_F( A, B, relative );
}

/////////////////////////////////////////////////////////////////////////
//!
//! \class  TSpectralNorm
//! \brief  computes spectral norm ‖·‖₂ of matrix (or linear operator)
//!
class TSpectralNorm : public TMatrixNorm
{
private:
    //! minimal number of steps for power iteration
    const uint  _min_it;

    //! maximal number of steps for power iteration
    const uint  _max_it;

    //! bound for relative error during power iteration
    const real  _rel_eps;

    //! bound for absolute error during power iteration
    const real  _abs_eps;
    
public:
    //////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! construct spectral norm object with given stop criterion
    TSpectralNorm ( const uint  min_it  = 10,
                    const uint  max_it  = 100,
                    const real  rel_eps = real(1e-5),
                    const real  abs_eps = real(100) * std::numeric_limits<real>::epsilon() );

    //! dtor
    virtual ~TSpectralNorm () {}

    //////////////////////////////////////////
    //
    // matrix norms
    //

    virtual real norm        ( const TMatrix *  A ) const
    {
        return norm( cptrcast( A, TLinearOperator ) );
    }

    virtual real diff_norm  ( const TMatrix *  A,
                              const TMatrix *  B,
                              const bool       rel = true ) const
    {
        return diff_norm( cptrcast( A, TLinearOperator ),
                          cptrcast( B, TLinearOperator ),
                          rel );
    }

    //////////////////////////////////////////
    //
    // linear operator norms
    //

    //! compute spectral norm of A
    virtual real norm        ( const TLinearOperator *  A ) const;
    
    //! compute difference norm ‖A-B‖₂ or ‖A-B‖₂/‖A‖₂ (if \a rel == true)
    virtual real diff_norm   ( const TLinearOperator *  A,
                               const TLinearOperator *  B,
                               const bool               rel = true ) const;
    
    //! compute spectral norm of A^-1
    virtual real inv_norm    ( const TLinearOperator *  A ) const;
    
    //! compute ‖ I - B A ‖₂
    virtual real inv_approx  ( const TLinearOperator *  A,
                               const TLinearOperator *  B ) const;
    
protected:
    //! return true if stop condition of power iteration is fulfilled
    bool converged ( const complex  new_val,
                     const complex  old_val,
                     const uint     step ) const;
};

//!
//! \fn     norm_2
//! \brief  compute spectral norm of A (functional form of \see TSpectralNorm)
//!
inline real
norm_2 ( const TLinearOperator * A )
{
    TSpectralNorm  mnorm;

    return mnorm.norm( A );
}

//!
//! \fn     diff_norm_2
//! \brief  compute spectral norm of \a A - \a B (functional form of \see TSpectralNorm)
//!
inline real
diff_norm_2 ( const TLinearOperator *  A,
              const TLinearOperator *  B,
              const bool               relative = true )
{
    TSpectralNorm  mnorm;

    return mnorm.diff_norm( A, B, relative );
}
inline real
diffnorm_2 ( const TLinearOperator *  A,
             const TLinearOperator *  B,
             const bool               relative = true )
{
    return diff_norm_2( A, B, relative );
}

//!
//! \fn     inv_approx_2
//! \brief  compute ‖ I - B A ‖₂ (functional form of \see TSpectralNorm)
inline real
inv_approx_2  ( const TLinearOperator *  A,
                const TLinearOperator *  B )
{
    TSpectralNorm  mnorm;

    return mnorm.inv_approx( A, B );
}
    

/////////////////////////////////////////////////////////////////////////
//!
//! \class  TInfinityNorm
//! \brief  computes infinity norm ‖·‖_∞ of a matrix
//!
class TInfinityNorm : public TMatrixNorm
{
public:
    //////////////////////////////////////////
    //
    // constructor and destructor
    //

    TInfinityNorm () {}

    virtual ~TInfinityNorm () {}

    //////////////////////////////////////////
    //
    // matrix-norms
    //

    //! compute Infinity norm of A
    virtual real norm       ( const TMatrix *  A ) const;

    //! compute difference norm ‖A-B‖_∞ or ‖A-B‖_∞/‖A‖_∞ (if \a rel == true)
    virtual real diff_norm  ( const TMatrix *  A,
                              const TMatrix *  B,
                              const bool       rel = true ) const;
};

//!
//! \fn     norm_inf
//! \brief  compute Infinity norm of A (functional form of \see TInfinityNorm)
//!
inline real
norm_inf ( const TMatrix *  A )
{
    TInfinityNorm  mnorm;

    return mnorm.norm( A );
}

//!
//! \fn     diff_norm_inf
//! \brief  compute Infinity norm of \a A - \a B (functional form of \see TInfinityNorm)
//!
inline real
diff_norm_inf ( const TMatrix *  A,
                const TMatrix *  B,
                const bool       relative = true )
{
    TInfinityNorm  mnorm;

    return mnorm.diff_norm( A, B, relative );
}
inline real
diffnorm_inf ( const TMatrix *  A,
               const TMatrix *  B,
               const bool       relative = true )
{
    return diff_norm_inf( A, B, relative );
}

/////////////////////////////////////////////////////////////////////////
//!
//! \class  TRowMatrixNorm
//! \brief  computes norm for each row of a matrix
//!
class TRowMatrixNorm
{
public:
    //////////////////////////////////////////
    //
    // constructor and destructor
    //

    TRowMatrixNorm () {}

    virtual ~TRowMatrixNorm () {}

    //////////////////////////////////////////
    //
    // matrix-norms
    //

    //! compute the row norms and return them in \a norm
    virtual void norm ( const TMatrix *  A,
                        TScalarVector &  norm ) const;
};

/////////////////////////////////////////////////////////////////////////
//!
//! \class  TColumnMatrixNorm
//! \brief  computes norm for each column of the matrix
//!
class TColumnMatrixNorm
{
public:
    //////////////////////////////////////////
    //
    // constructor and destructor
    //

    TColumnMatrixNorm () {}

    virtual ~TColumnMatrixNorm () {}

    //////////////////////////////////////////
    //
    // matrix-norms
    //

    //! compute the column norms and return them in \a norm
    virtual void norm ( const TMatrix *  A,
                        TScalarVector &  norm ) const;
};

/////////////////////////////////////////////////////////////////////////
//
// misc. other norm functions
//

//!
//! Assuming \a A is a LU (or Cholesky or LDL) factorisation matrix
//! (as TFacInvMatrix) the condest value ∥(LU)^-1 e∥∞, with e = (1,…,1),
//! is computed
//!
real
condest        ( const TLinearOperator *  LU );

//!
//! return norm of residual, i.e., |Ax-b|₂
//!
real
residual_norm  ( const TLinearOperator *  A,
                 const TVector *          x,
                 const TVector *          b );

//!
//! compute singular values of given matrix
//!
void
singular_values ( const TMatrix *         A,
                  BLAS::Vector< real > &  S );

BLAS::Vector< real >
singular_values ( const TMatrix *         A );

}// namespace HLIB

#endif  // __HLIB_MAT_NORM_HH
