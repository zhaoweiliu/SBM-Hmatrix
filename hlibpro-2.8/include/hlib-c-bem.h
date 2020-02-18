#ifndef __HLIB_C_BEM_H
#define __HLIB_C_BEM_H
/*
 * Project     : HLib
 * File        : hlib-c-bem.hh
 * Description : C interface to HLib (BEM functions)
 * Author      : Ronald Kriemann
 * Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
 */

#include <math.h>

#include <hlib-c.h>

#ifdef __cplusplus
extern "C" {
#endif

/***************************************************************
 ***************************************************************
 **
 ** types
 **
 ***************************************************************
 ***************************************************************/

/*
 * represents a grid, function spaces and bilinear forms
 */
typedef struct HLIBPF(grid_s) *     HLIBPF(grid_t);
typedef struct HLIBPF(fnspace_s) *  HLIBPF(fnspace_t);
typedef struct HLIBPF(bemrbf_s) *   HLIBPF(bemrbf_t);
typedef struct HLIBPF(bemcbf_s) *   HLIBPF(bemcbf_t);

/*
 * grid based function for BEM applications on a 3d surface,
 * evaluated at points <x> and (optional) normal direction <n>
 * result is returned in <res>
 */
typedef void (* HLIBPF(gridfn_t))   ( const double *       x,
                                      const double *       n,
                                      void *               arg,
                                      double *             res );
typedef void (* HLIBPF(gridcfn_t))  ( const double *       x,
                                      const double *       n,
                                      void *               arg,
                                      HLIBPF(complex_t) *  res );

/***************************************************************
 ***************************************************************
 **
 ** grids
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 **
 ** management
 **
 ************************************************/

/* free grid and all related resources */
HLIB_FNDECL void
HLIBPF(grid_free)                  ( HLIBPF(grid_t)            grid,
                                     int *                     info );

/* return amount of memory in bytes used by grid */
HLIB_FNDECL size_t
HLIBPF(grid_bytesize)              ( const HLIBPF(grid_t)      grid,
                                     int *                     info );

/************************************************
 **
 ** grid transformation
 **
 ************************************************/

/* translate <grid> by given vector <t> in R^3 */
HLIB_FNDECL void
HLIBPF(grid_translate)             ( HLIBPF(grid_t)            grid,
                                     const double              t[3],
                                     int *                     info );

/* scale <grid> coordinates by x, y and z components
 * defined by given vector <s> */
HLIB_FNDECL void
HLIBPF(grid_scale)                 ( HLIBPF(grid_t)            grid,
                                     const double              s[3],
                                     int *                     info );

/* rotate <grid> coordinates around <v> by by angle
 * <alpha> in radians */
HLIB_FNDECL void
HLIBPF(grid_rotate)                ( HLIBPF(grid_t)            grid,
                                     const double              v[3],
                                     const double              alpha,
                                     int *                     info );

/************************************************
 **
 ** input
 **
 ************************************************/

/* read grid in HLib format from file <filename> */
HLIB_FNDECL HLIBPF(grid_t)
HLIBPF(hformat_load_grid)          ( const char *              filename,
                                     int *                     info );

/* read grid in Ply format from file <filename> */
HLIB_FNDECL HLIBPF(grid_t)
HLIBPF(ply_load_grid)              ( const char *              filename,
                                     int *                     info );

/*
 * autodetect format
 */
HLIB_FNDECL HLIBPF(grid_t)
HLIBPF(load_grid)                  ( const char *              filename,
                                     int *                     info );

/************************************************
 **
 ** visualisation
 **
 ************************************************/

/* print grid in PostScript format to file <filename>,
 * <view> defines viewing direction                   */
HLIB_FNDECL void
HLIBPF(grid_print_ps)              ( const HLIBPF(grid_t)      grid,
                                     const double              view[3],
                                     const char *              filename,
                                     int *                     info );

/* print grid in VRML format to file <filename> */
HLIB_FNDECL void
HLIBPF(grid_print_vrml)            ( const HLIBPF(grid_t)      grid,
                                     const char *              filename,
                                     int *                     info );

/* print grid in VTK format to file <filename> */
HLIB_FNDECL void
HLIBPF(grid_print_vtk)             ( const HLIBPF(grid_t)      grid,
                                     const char *              filename,
                                     int *                     info );

/* print grid in PostScript format to file <filename>
 * - values in <v> for function in \a fnspace are
 *   used to colour triangles         
 * - <view> defines viewing direction */
HLIB_FNDECL void
HLIBPF(grid_val_print_ps)          ( const HLIBPF(grid_t)      grid,
                                     const HLIBPF(fnspace_t)   fnspace,
                                     const HLIBPF(vector_t)    v,
                                     const double              view[3],
                                     const char *              filename,
                                     int *                     info );

/* print grid in VRML format to file <filename>
 * values in <v> for function in \a fnspace are
 * used to colour triangles */
HLIB_FNDECL void
HLIBPF(grid_val_print_vrml)        ( const HLIBPF(grid_t)      grid,
                                     const HLIBPF(fnspace_t)   fnspace,
                                     const HLIBPF(vector_t)    v,
                                     const char *              filename,
                                     int *                     info );

/* print grid in VTK format to file <filename>
 * values in <v> for function in \a fnspace are
 * used to colour triangles */
HLIB_FNDECL void
HLIBPF(grid_val_print_vtk)         ( const HLIBPF(grid_t)      grid,
                                     const HLIBPF(fnspace_t)   fnspace,
                                     const HLIBPF(vector_t)    v,
                                     const char *              filename,
                                     int *                     info );

/***************************************************************
 ***************************************************************
 **
 ** function spaces
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 **
 ** construction
 **
 ************************************************/

/*
 * build a functions space with piecewise constant basis functions over <grid >
 */
HLIB_FNDECL HLIBPF(fnspace_t)
HLIBPF(fnspace_build_const)        ( const HLIBPF(grid_t)      grid,
                                     int *                     info );

/*
 * build a functions space with piecewise linear basis functions over <grid >
 */
HLIB_FNDECL HLIBPF(fnspace_t)
HLIBPF(fnspace_build_linear)       ( const HLIBPF(grid_t)      grid,
                                     int *                     info );

/*
 * build a functions space for constant normal linear tangential (CN/LT) 
 * edge elements
 */
HLIB_FNDECL HLIBPF(fnspace_t)
HLIBPF(fnspace_build_constedge)    ( const HLIBPF(grid_t)      grid,
                                     int *                     info );

/*
 * return the dimension of the given function space <fnspace>,
 *  e.g. the number of unknowns
 */
HLIB_FNDECL unsigned int
HLIBPF(fnspace_dim)                ( const HLIBPF(fnspace_t)   fnspace,
                                     int *                     info );

/*
 * return coordinates of all indices/unknowns in functions space
 */
HLIB_FNDECL HLIBPF(coord_t)
HLIBPF(fnspace_coord)              ( const HLIBPF(fnspace_t)   fnspace,
                                     int *                     info );

/************************************************
 **
 ** management
 **
 ************************************************/

/* free function space and all related resources */
HLIB_FNDECL void
HLIBPF(fnspace_free)               ( HLIBPF(fnspace_t)         fnspace,
                                     int *                     info );

/* return amount of memory in bytes used by function space */
HLIB_FNDECL size_t
HLIBPF(fnspace_bytesize)           ( const HLIBPF(fnspace_t)   fnspace,
                                     int *                     info );

/************************************************
 **
 ** evaluation
 **
 ************************************************/

/*
 * evaluate given BEM function <fn> at indices of <fnspace> and
 * build corresponding vector; <arg> is the optional argument to <fn>
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(fnspace_eval)               ( const HLIBPF(fnspace_t)   fnspace,
                                     const HLIBPF(gridfn_t)    fn,
                                     void *                    arg,
                                     int *                     info );

HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(fnspace_ceval)              ( const HLIBPF(fnspace_t)   fnspace,
                                     const HLIBPF(gridcfn_t)   fn,
                                     void *                    arg,
                                     int *                     info );

/*
 * build vector for right-hand-side function <rhs> and indices
 * defined by <fnspace>, e.g. evaluate ∫ <φ_i,rhs>, ∀ i
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(fnspace_build_rhs)          ( const HLIBPF(fnspace_t)      fnspace,
                                     const HLIBPF(gridfn_t)       rhs,
                                     void *                       arg,
                                     int *                        info );

HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(fnspace_build_crhs)         ( const HLIBPF(fnspace_t)      fnspace,
                                     const HLIBPF(gridcfn_t)      rhs,
                                     void *                       arg,
                                     int *                        info );

/***************************************************************
 ***************************************************************
 **
 ** bilinear forms
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 **
 ** construction
 **
 ************************************************/

/* return bilinear form for mass matrix over <ansatz> and <test> space */
HLIB_FNDECL HLIBPF(bemrbf_t)
HLIBPF(bembf_mass)                 ( const HLIBPF(fnspace_t)   ansatz,
                                     const HLIBPF(fnspace_t)   test,
                                     int *                     info );
    
/* return bilinear form for Laplace SLP over <ansatz> and <test> space */
HLIB_FNDECL HLIBPF(bemrbf_t)
HLIBPF(bembf_laplace_slp)          ( const HLIBPF(fnspace_t)   ansatz,
                                     const HLIBPF(fnspace_t)   test,
                                     int *                     info );
    
/* return bilinear form for Laplace DLP over <ansatz> and <test> space */
HLIB_FNDECL HLIBPF(bemrbf_t)
HLIBPF(bembf_laplace_dlp)          ( const HLIBPF(fnspace_t)   ansatz,
                                     const HLIBPF(fnspace_t)   test,
                                     int *                     info );
    
/* return bilinear form for Helmholtz SLP with wavenumber <kappa>
 * over <ansatz> and <test> space */
HLIB_FNDECL HLIBPF(bemcbf_t)
HLIBPF(bembf_helmholtz_slp)        ( const HLIBPF(fnspace_t)   ansatz,
                                     const HLIBPF(fnspace_t)   test,
                                     const HLIBPF(complex_t)   kappa,
                                     int *                     info );
    
/* return bilinear form for Helmholtz DLP with wavenumber <kappa>
 * over <ansatz> and <test> space */
HLIB_FNDECL HLIBPF(bemcbf_t)
HLIBPF(bembf_helmholtz_dlp)        ( const HLIBPF(fnspace_t)   ansatz,
                                     const HLIBPF(fnspace_t)   test,
                                     const HLIBPF(complex_t)   kappa,
                                     int *                     info );
    
/* return bilinear form for acoustic scattering with wavenumber <kappa>
 * over <ansatz> and <test> space */
HLIB_FNDECL HLIBPF(bemcbf_t)
HLIBPF(bembf_acousticscatter)      ( const HLIBPF(fnspace_t)   ansatz,
                                     const HLIBPF(fnspace_t)   test,
                                     const HLIBPF(complex_t)   kappa,
                                     int *                     info );
    
/* return bilinear form for Electric Field Integral Equation for Maxwell
 *  with wavenumber <kappa> over <ansatz> and <test> space */
HLIB_FNDECL HLIBPF(bemcbf_t)
HLIBPF(bembf_maxwell_efie)         ( const HLIBPF(fnspace_t)   ansatz,
                                     const HLIBPF(fnspace_t)   test,
                                     const HLIBPF(complex_t)   kappa,
                                     const HLIBPF(real_t)      eta,
                                     int *                     info );
    
/* return bilinear form for Magnetic Field Integral Equation for Maxwell
 * with wavenumber <kappa> over <ansatz> and <test> space */
HLIB_FNDECL HLIBPF(bemcbf_t)
HLIBPF(bembf_maxwell_mfie)         ( const HLIBPF(fnspace_t)   ansatz,
                                     const HLIBPF(fnspace_t)   test,
                                     const HLIBPF(complex_t)   kappa,
                                     int *                     info );
    
/************************************************
 **
 ** matrix construction
 **
 ************************************************/

/*
 * call back functions for "hlib_matrix_build_coeff";
 * the <arg> parameter must be the bilinear form
 */
HLIB_FNDECL void
HLIBPF(eval_bemrbf_cb)  ( const size_t         n,
                          const int *          rowidx,
                          const size_t         m,
                          const int *          colidx,
                          HLIBPF(real_t) *     matrix,
                          void *               arg );

HLIB_FNDECL void
HLIBPF(eval_bemcbf_cb)  ( const size_t         n,
                          const int *          rowidx,
                          const size_t         m,
                          const int *          colidx,
                          HLIBPF(complex_t) *  matrix,
                          void *               arg );

/************************************************
 **
 ** management
 **
 ************************************************/

/* free bilinear form */
HLIB_FNDECL void
HLIBPF(bemrbf_free)                ( HLIBPF(bemrbf_t)          bembf,
                                     int *                     info );
HLIB_FNDECL void
HLIBPF(bemcbf_free)                ( HLIBPF(bemcbf_t)          bembf,
                                     int *                     info );

/***************************************************************
 ***************************************************************
 **
 ** Misc.
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 **
 ** Quadrature Rules
 **
 ************************************************/

/*
 * return Gaussian quadrature points (<points>) and
 * weights (<weights>) of order <order> for the interval [0,1]
 */
HLIB_FNDECL void
HLIBPF(gauss_quad_1d)              ( const unsigned int        order,
                                     double *                  points,  /* array of dim <order>     */
                                     double *                  weights, /* array of dim <order>     */
                                     int *                     info );

/*
 * 3d triangle Gauss quadrature rules
 *   - <tri> contains the three coordinates for the triangle
 *   - <qtri> will contain the 3d quadrature points
 *   - <qwghts> will contain the quadrature weights
 *   - <qtri> and <qwghts> must have dimension <order>^2
 */
HLIB_FNDECL void
HLIBPF(gauss_quad_tri)             ( const unsigned int        order,
                                     const double *            tri[3],  /* coord. of triangle        */
                                     double *                  qtri[3], /* quad. points for triangle */
                                     double *                  qwghts,  /* quadrature weights        */
                                     int *                     info );

/*
 * 3d triangle pair quadrature rules by Stefan Sauter
 *   - constructs quadrature points and weights for both triangles
 *   - <tri1> and <tri2> contain the 3d coord. of the triangles
 *   - <qtri1> and <qtri2> will contain the 3d quadrature points
 *   - <qwghts> will contain the quadrature weights
 *   - <qtri1>, <qtri2> and <qwghts> must be of dimension 6*<order>^4 !!!
 */
HLIB_FNDECL void
HLIBPF(sauter_quad)                ( const unsigned int        order,
                                     const double *            tri1[3], /* coord. of triangle 1         */
                                     const double *            tri2[3], /* coord. of triangle 2         */
                                     double *                  qtri1[3],/* quadrature points for tri. 1 */
                                     double *                  qtri2[3],/* quadrature points for tri. 2 */
                                     double *                  qwghts,  /* quadrature weights           */
                                     int *                     info );

/*
 * 2d triangle pair quadrature rules by Stefan Sauter
 *   - ncommon defines number of common vertices in triangles
 *   - order defines quadrature order
 *   - constructs 2d points for both triangles (<tri1_pts>,
 *     <tri2_pts>) and corresponding weights (<weights>)
 *   - coordinates are for triangle (0,0) -- (1,0) -- (0,1)
 */

/* for equal triangles */
HLIB_FNDECL void
HLIBPF(sauter_quad_eq)             ( const unsigned int        order,
                                     double *                  tri1_pts[2], /* array of dim 6*<order>^4 */
                                     double *                  tri2_pts[2], /* array of dim 6*<order>^4 */
                                     double *                  weights,     /* array of dim 6*<order>^4 */
                                     int *                     info );

/* for triangles with common edge */
HLIB_FNDECL void
HLIBPF(sauter_quad_edge)           ( const unsigned int        order,
                                     double *                  tri1_pts[2], /* array of dim 5*<order>^4 */
                                     double *                  tri2_pts[2], /* array of dim 5*<order>^4 */
                                     double *                  weights,     /* array of dim 5*<order>^4 */
                                     int *                     info );

/* for triangles with common vertex */
HLIB_FNDECL void
HLIBPF(sauter_quad_vtx)            ( const unsigned int        order,
                                     double *                  tri1_pts[2], /* array of dim 2*<order>^4 */
                                     double *                  tri2_pts[2], /* array of dim 2*<order>^4 */
                                     double *                  weights,     /* array of dim 2*<order>^4 */
                                     int *                     info );

/* for separated triangles */
HLIB_FNDECL void
HLIBPF(sauter_quad_sep)            ( const unsigned int        order,
                                     double *                  tri1_pts[2], /* array of dim <order>^4   */
                                     double *                  tri2_pts[2], /* array of dim <order>^4   */
                                     double *                  weights,     /* array of dim <order>^4   */
                                     int *                     info );

#ifdef __cplusplus
}/* extern */
#endif

#endif  /* __HLIB_C_BEM_H */
