#ifndef __HLIB_TSOR_HH
#define __HLIB_TSOR_HH
//
// Project     : HLib
// File        : TSOR.hh
// Description : classes representing the SOR and Gauss-Seidel iteration matrix
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/matrix/TMatrix.hh"
#include "hpro/algebra/solve_types.hh"

namespace HLIB
{

// local matrix type
DECLARE_TYPE( TSOR );

// type of SOR iteration
enum sor_type_t
{
    SOR_FORWARD  = 0,
    SOR_BACKWARD,
    SOR_SYMMETRIC
};

//!
//! \ingroup Matrix_Module
//! \class   TSOR
//! \brief   implements SOR preconditioner for sparse matrices
//!
//!          TSOR provides application of a SOR type preconditioner for a given
//!          sparse matrix \f$A = D - E - F\f$, with diagonal D, strictly lower triangular
//!          matrix E and strictly upper triangular matrix F.
//!
//!          For forward SOR, the applied operator is \f$ \omega ( D - \omega E )^{-1}\f$,
//!          for backward SOR \f$ \omega ( D - \omega F )^{-1}\f$ and the combination of both
//!          for the symmetric SOR.
//!
class TSOR : public TLinearOperator
{
private:

    //! @cond
    
    // matrix defining coefficients for SOR
    const TSparseMatrix *  _mat_A;

    // relaxation factor
    const real             _omega;

    // damping factor
    const real             _damping;

    // type of SOR step
    const sor_type_t       _sor_type;

    //! @endcond
    
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! constructs SOR preconditioner of type \a sor_type using coefficients
    //! as defined by sparse matrix \a A
    TSOR ( const TSparseMatrix *  A,
           const sor_type_t       sor_type,
           const real             omega   = 1.0,
           const real             damping = 1.0 );

    // dtor
    virtual ~TSOR ();
    
    /////////////////////////////////////////////////
    //
    // access internal data
    //

    //! return internal sparse matrix
    const TMatrix *  matrix          () const  { return _mat_A; }

    //! return SOR type
    sor_type_t       sor_type        () const  { return _sor_type; }

    //! return damping factor
    real             damping_factor  () const  { return _damping; }
                                       
    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    bool  is_complex      () const
    {
        return _mat_A->is_complex();
    }
    
    //! return true, of operator is self adjoint
    bool  is_self_adjoint () const
    {
        return ( _sor_type == SOR_SYMMETRIC ) && _mat_A->is_hermitian();
    }
    
    /////////////////////////////////////////////////
    //
    // linear operator mapping
    //

    //!
    //! mapping function of linear operator \f$A\f$, e.g. \f$ y := A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply       ( const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;

    //!
    //! mapping function with update: \f$ y := y + \alpha A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply_add   ( const real       alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;
    virtual void  capply_add  ( const complex    alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;

    virtual void  apply_add   ( const real       alpha,
                                const TMatrix *  X,
                                TMatrix *        Y,
                                const matop_t    op = apply_normal ) const;
    
    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const real                       alpha,
                                const BLAS::Vector< real > &     x,
                                BLAS::Vector< real > &           y,
                                const matop_t                    op = apply_normal ) const;
    virtual void  apply_add   ( const complex                    alpha,
                                const BLAS::Vector< complex > &  x,
                                BLAS::Vector< complex > &        y,
                                const matop_t                    op = apply_normal ) const;

    ///////////////////////////////////////////////////////////
    //
    // access to vector space elements
    //

    //! return dimension of domain
    virtual size_t  domain_dim () const { return _mat_A->range_dim(); }
    
    //! return dimension of range
    virtual size_t  range_dim  () const { return _mat_A->domain_dim(); }

    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector > { return _mat_A->row_vector(); }

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector > { return _mat_A->col_vector(); }


    HLIB_RTTI_DERIVED( TSOR, TLinearOperator );
};



// local matrix type
DECLARE_TYPE( TGaussSeidel );

// type of Gauss-Seidel iteration
enum gs_type_t
{
    GS_FORWARD  = 0,
    GS_BACKWARD,
    GS_SYMMETRIC
};

//!
//! \ingroup Matrix_Module
//! \class   TGaussSeidel
//! \brief   implements Gauss-Seidel preconditioner
//!
//!          TGaussSeidel provides application of a Gauss-Seidel type preconditioner 
//!          for a given matrix \f$A = D - E - F\f$, with diagonal D, strictly lower triangular
//!          matrix E and strictly upper triangular matrix F.
//!
//!          For forward GS, the applied operator is \f$ ( D - E )^{-1}\f$,
//!          for backward GS \f$ ( D - F )^{-1}\f$ and the combination of both
//!          for the symmetric GS, i.e., \f$ ( D - F )^{-1} D ( D - E )^{-1}\f$.
//!
//!          TGaussSeidel implements point-wise and block-wise GS steps. However, point-wise
//!          GS is only supported for sparse matrices while block-wise GS is only supported
//!          for H-matrices
//!
class TGaussSeidel : public TLinearOperator
{
private:

    //! @cond
    
    // matrix defining coefficients for GaussSeidel
    const TMatrix *             _mat_A;

    // damping factor
    const real                  _damping;

    // type of GS-step
    const gs_type_t             _gs_type;

    // inverse of diagonal (internal use)
    std::unique_ptr< TMatrix >  _blockdiag;
    
    //! @endcond
    
public:
    ///////////////////////////////////////////
    //
    // constructor and destructor
    //

    //! constructs Gauss-Seidel preconditioner of type \a gs_type
    TGaussSeidel ( const TMatrix *    A,
                   const gs_type_t    gs_type,
                   const real         damping = 1.0 );

    // dtor
    virtual ~TGaussSeidel ();
    
    /////////////////////////////////////////////////
    //
    // access internal data
    //

    //! return internal sparse matrix
    const TMatrix *  matrix          () const  { return _mat_A; }

    //! return GaussSeidel type
    gs_type_t        gs_type         () const  { return _gs_type; }

    //! return damping factor
    real             damping_factor  () const  { return _damping; }
                                       
    ///////////////////////////////////////////////////////////
    //
    // linear operator properties
    //

    //! return true, if field type is complex
    bool  is_complex      () const
    {
        return _mat_A->is_complex();
    }
    
    //! return true, of operator is self adjoint
    bool  is_self_adjoint () const
    {
        return ( _gs_type == GS_SYMMETRIC ) && _mat_A->is_hermitian();
    }
    
    /////////////////////////////////////////////////
    //
    // linear operator mapping
    //

    //!
    //! mapping function of linear operator \f$A\f$, e.g. \f$ y := A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply       ( const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;

    //!
    //! mapping function with update: \f$ y := y + \alpha A(x)\f$.
    //! Depending on \a op, either \f$A\f$, \f$A^T\f$ or \f$A^H\f$ is applied.
    //!
    virtual void  apply_add   ( const real       alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;
    virtual void  capply_add  ( const complex    alpha,
                                const TVector *  x,
                                TVector *        y,
                                const matop_t    op = apply_normal ) const;

    virtual void  apply_add   ( const real       alpha,
                                const TMatrix *  X,
                                TMatrix *        Y,
                                const matop_t    op = apply_normal ) const;
    
    //! same as above but only the dimension of the vector spaces is tested,
    //! not the corresponding index sets
    virtual void  apply_add   ( const real                       alpha,
                                const BLAS::Vector< real > &     x,
                                BLAS::Vector< real > &           y,
                                const matop_t                    op = apply_normal ) const;
    virtual void  apply_add   ( const complex                    alpha,
                                const BLAS::Vector< complex > &  x,
                                BLAS::Vector< complex > &        y,
                                const matop_t                    op = apply_normal ) const;

    ///////////////////////////////////////////////////////////
    //
    // access to vector space elements
    //

    //! return dimension of domain
    virtual size_t  domain_dim () const { return _mat_A->range_dim(); }
    
    //! return dimension of range
    virtual size_t  range_dim  () const { return _mat_A->domain_dim(); }

    //! return vector in domain space
    virtual auto    domain_vector  () const -> std::unique_ptr< TVector > { return _mat_A->row_vector(); }

    //! return vector in range space
    virtual auto    range_vector   () const -> std::unique_ptr< TVector > { return _mat_A->col_vector(); }


    HLIB_RTTI_DERIVED( TGaussSeidel, TLinearOperator );
};

}// namespace HLIB

#endif  // __HLIB_TSOR_HH
