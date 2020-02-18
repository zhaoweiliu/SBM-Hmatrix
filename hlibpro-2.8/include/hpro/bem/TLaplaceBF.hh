#ifndef __HLIB_TLAPLACEBF_HH
#define __HLIB_TLAPLACEBF_HH
//
// Project     : HLib
// File        : TLAPLACEBF.hh
// Description : bilinear forms for Laplace operator
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/algebra/TLowRankApx.hh"
#include "hpro/bem/TQuadBEMBF.hh"
#include "hpro/bem/TQuadHCAGenFn.hh"

namespace HLIB
{
    
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// standard bilinear forms for Laplace
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////
//
// single layer potential of Laplace operator
//
////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   TLaplaceSLPBF
//! \brief   Bilinear form for Laplace single layer potential
//!
//!          TLaplaceSLPBF implements the bilinear form for the Laplace
//!          single layer potential with the kernel function
//!          \f[ \frac{1}{\|x-y\|_2} \f]
//!          i.e. for the integral equation
//!          \f[ 4 \pi \int_{\Gamma} \frac{u(y)}{\|x-y\|_2} dy = f(x) \f]
//!
template < typename  T_ansatzsp,
           typename  T_testsp >
class TLaplaceSLPBF : public TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, real >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = real;

private:

    //
    // kernel function impl.
    //

    void  ( * _kernel_fn ) ( const TGrid::triangle_t &    tri0,
                             const TGrid::triangle_t &    tri1,
                             const tripair_quad_rule_t *  rule,
                             std::vector< real > &        values,
                             const ansatzsp_t *           ansatz_sp,
                             const testsp_t *             test_sp );

public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //

    TLaplaceSLPBF ( const ansatzsp_t *  aansatzsp,
                    const testsp_t *    atestsp,
                    const uint          quad_order = CFG::BEM::quad_order );

    virtual ~TLaplaceSLPBF () {}

    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;
            
protected:
    //
    // eval kernel function at quadrature points
    //
    virtual void  eval_kernel  ( const idx_t                  tri0idx,
                                 const idx_t                  tri1idx,
                                 const TGrid::triangle_t &    tri0,
                                 const TGrid::triangle_t &    tri1,
                                 const tripair_quad_rule_t *  quad_rule,
                                 std::vector< real > &        values ) const;

};

////////////////////////////////////////////////////
//
// double layer potential of Laplace operator
//
////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   TLaplaceDLPBF
//! \brief   Bilinear form for Laplace double layer potential
//!
//!          TLaplaceDLPBF implements the bilinear form for the Laplace
//!          double layer potential with the kernel function
//!          \f[ \frac{\langle n, y-x \rangle}{\|x-y\|_2^3} \f].
//!
template < typename  T_ansatzsp,
           typename  T_testsp >
class TLaplaceDLPBF : public TInvarBasisQuadBEMBF< T_ansatzsp, T_testsp, real >
{
public:
    //
    // template types as internal types
    //
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = real;
        
private:

    //
    // kernel function impl.
    //

    void  ( * _kernel_fn ) ( const TGrid::triangle_t &   tri0,
                             const TGrid::triangle_t &   tri1,
                             const tripair_quad_rule_t * rule,
                             std::vector< real > &       values,
                             const ansatzsp_t *          ansatz_sp,
                             const testsp_t *            test_sp,
                             const T3Point &             n );

public:
    //////////////////////////////////////
    //
    // constructor and destructor
    //
    
    TLaplaceDLPBF ( const ansatzsp_t *  aansatzsp,
                    const testsp_t *    atestsp,
                    const uint          quad_order = CFG::BEM::quad_order );

    virtual ~TLaplaceDLPBF () {}
    
    // return format of bilinear form, e.g. symmetric
    matform_t  format () const;
            
protected:
    //
    // eval kernel function at quadrature points
    //
    virtual void  eval_kernel  ( const idx_t                 tri0idx,
                                 const idx_t                 tri1idx,
                                 const TGrid::triangle_t &   tri0,
                                 const TGrid::triangle_t &   tri1,
                                 const tripair_quad_rule_t * quad_rule,
                                 std::vector< real > &       values ) const;
};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
//
// HCA functions for Laplace
//
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

//!
//! \ingroup BEM_Module
//! \class   TLaplaceSLPGenFn
//! \brief   kernel generator function for Laplace SLP
//!
template < typename T_ansatzsp,
           typename T_testsp >
class TLaplaceSLPGenFn : public TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, real >
{
public:
    //
    // template types as internal types
    //
    
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = real;
    
    //
    // inherit from base class
    //

    using  stat_t     = typename TQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::stat_t;

private:

    //
    // HCA function impl.
    //

    void ( * _eval_dxy_impl ) ( const tri_quad_rule_t &  quad_rule,
                                const T3Point            x[3],
                                const T3Point &          y,
                                std::vector< real > &    values );
    
public:
    TLaplaceSLPGenFn ( const ansatzsp_t *    ansatzsp,
                       const testsp_t *      testsp,
                       const TPermutation *  row_perm_i2e,
                       const TPermutation *  col_perm_i2e,
                       const uint            quad_order = CFG::BEM::quad_order );

    //
    // evaluate generator function at (\a x, \a y)
    //
    real
    eval ( const T3Point &  x,
           const T3Point &  y ) const
    {
        const real factor  = real( 1.0 / ( 4.0 * Math::pi< double >() ) );
        const real sq_dist = real( Math::square( x[0] - y[0] ) +
                                   Math::square( x[1] - y[1] ) +
                                   Math::square( x[2] - y[2] ) );
        
        return factor * Math::rsqrt( sq_dist );
    }

protected:
    //! Evaluate \f$ D_x \gamma(x, y) \f$ on with \f$x\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dx  ( const idx_t              tri_idx,
                             const T3Point &          y,
                             const tri_quad_rule_t &  quad_rule,
                             std::vector< real > &    values ) const
    {
        const TGrid *            grid = this->ansatz_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            x[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };

        _eval_dxy_impl( quad_rule, x, y, values );
    }
    
    //! Evaluate \f$ D_y \gamma(x, y) \f$ on with \f$y\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dy  ( const T3Point &          x,
                             const idx_t              tri_idx,
                             const tri_quad_rule_t &  quad_rule,
                             std::vector< real > &    values ) const
    {
        const TGrid *            grid = this->test_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            y[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
    
        _eval_dxy_impl( quad_rule, y, x, values );
    }
};

//!
//! \ingroup BEM_Module
//! \class   TLaplaceDLPGenFn
//! \brief   kernel generator function for Laplace DLP
//!
template < typename T_ansatzsp,
           typename T_testsp >
class TLaplaceDLPGenFn : public TInvarBasisQuadHCAGenFn< T_ansatzsp, T_testsp, real >
{
public:
    //
    // template types as internal types
    //
    
    using  ansatzsp_t = T_ansatzsp;
    using  testsp_t   = T_testsp;
    using  value_t    = real;
    
    //
    // inherit from base class
    //

    using  stat_t     = typename TQuadHCAGenFn< ansatzsp_t, testsp_t, value_t >::stat_t;
    
private:

    //
    // HCA function impl.
    //

    void ( * _eval_dx_impl ) ( const tri_quad_rule_t &  quad_rule,
                               const T3Point            x[3],
                               const T3Point &          y,
                               std::vector< real > &    values );
        
    void ( * _eval_dy_impl ) ( const tri_quad_rule_t &  quad_rule,
                               const T3Point &          x,
                               const T3Point            y[3],
                               const T3Point &          normal,
                               std::vector< real > &    values );
    
public:
    TLaplaceDLPGenFn ( const ansatzsp_t *    ansatzsp,
                       const testsp_t *      testsp,
                       const TPermutation *  row_perm_i2e,
                       const TPermutation *  col_perm_i2e,
                       const uint            quad_order = CFG::BEM::quad_order );

    //
    // evaluate generator function at (\a x, \a y)
    //
    real
    eval ( const T3Point & x,
           const T3Point & y ) const
    {
        const real  factor  = real( 1.0 / ( 4.0 * Math::pi< double >() ) );
        const real  sq_dist = real( Math::square( x[0] - y[0] ) +
                                    Math::square( x[1] - y[1] ) +
                                    Math::square( x[2] - y[2] ) );
        
        return factor * Math::rsqrt( sq_dist );
    }

protected:
    //! Evaluate \f$ D_x \gamma(x, y) \f$ on with \f$x\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dx  ( const idx_t              tri_idx,
                             const T3Point &          y,
                             const tri_quad_rule_t &  quad_rule,
                             std::vector< real > &    values ) const
    {
        const TGrid *            grid = this->ansatz_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            x[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
    
        _eval_dx_impl( quad_rule, x, y, values );
    }
    
    //! Evaluate \f$ D_y \gamma(x, y) \f$ on with \f$y\f$ defined by quadrature
    //! points on triangle \a tri_idx. The computed values for each quadrature
    //! point i are stored on \a values[i].
    virtual void  eval_dy  ( const T3Point &          x,
                             const idx_t              tri_idx,
                             const tri_quad_rule_t &  quad_rule,
                             std::vector< real > &    values ) const
    {
        const TGrid *            grid = this->test_space()->grid();
        const TGrid::triangle_t  tri  = grid->triangle( tri_idx );
        const T3Point            y[3] = { grid->vertex( tri.vtx[0] ),
                                          grid->vertex( tri.vtx[1] ),
                                          grid->vertex( tri.vtx[2] ) };
        const T3Point            normal( grid->tri_normal( tri_idx ) );
    
        _eval_dy_impl( quad_rule, x, y, normal, values );
    }
};

}// namespace HLIB

#endif  // __HLIB_TLAPLACEBF_HH
