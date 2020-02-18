#ifndef __HLIB_C_H
#define __HLIB_C_H
/*
 * Project     : HLib
 * File        : hlib-c.hh
 * Description : C interface to HLib
 * Author      : Ronald Kriemann
 * Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
 */

/** \defgroup CBinding_Module C Bindings
 *   @{
 *
 *   This modules defines functions for the usage of \HLIBpro in
 *   C programs instead of C++.
 *
 *   To include all standard \mcH-matrix related functions add
 *   \code
 *   #include <hlib-c.h>
 *   \endcode
 *   to your source files. For using BEM functions in C, you also
 *   have to add
 *   \code
 *   #include <hlib-c-bem.h>
 *   \endcode
 */

#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "hlib-config.h"

/**************************************************************
 **************************************************************
 **
 ** defines
 **
 **************************************************************
 **************************************************************/

/** @cond */

/* function call declaration */
#ifdef WINDOWS
#define HLIB_FNDECL __declspec(dllexport)
#else
#define HLIB_FNDECL
#endif

/** @endcond */

/**************************************************************
 **************************************************************
 **
 ** types
 **
 **************************************************************
 **************************************************************/

/***********************************************
 *
 * basic types
 *
 ***********************************************/

/* real data type */
#if HLIB_SINGLE_PREC == 1
typedef float  HLIBPF(real_t);
#else
typedef double HLIBPF(real_t);
#endif

/* complex data type */
typedef struct { HLIBPF(real_t) re, im; } HLIBPF(complex_t);

/*
 * (anonymous) representation of various types
 */
typedef struct HLIBPF(coord_s) *             HLIBPF(coord_t);
typedef struct HLIBPF(admcond_s) *           HLIBPF(admcond_t);

typedef struct HLIBPF(permutation_s) *       HLIBPF(permutation_t);

typedef struct HLIBPF(cluster_s) *           HLIBPF(cluster_t);
typedef struct HLIBPF(clustertree_s) *       HLIBPF(clustertree_t);
typedef struct HLIBPF(blockcluster_s) *      HLIBPF(blockcluster_t);
typedef struct HLIBPF(blockclustertree_s) *  HLIBPF(blockclustertree_t);
typedef struct HLIBPF(linearoperator_s) *    HLIBPF(linearoperator_t);
typedef struct HLIBPF(matrix_s) *            HLIBPF(matrix_t);
typedef struct HLIBPF(vector_s) *            HLIBPF(vector_t);
typedef struct HLIBPF(solver_s) *            HLIBPF(solver_t);

/***********************************************
 *
 * matrix operation type
 *
 ***********************************************/

/** defines operation on matrices */
typedef enum { HLIB_MATOP_NORM  = 'N',    /*!< normal matrix, no modification    */
               HLIB_MATOP_TRANS = 'T',    /*!< transposed matrix                 */
               HLIB_MATOP_ADJ   = 'C'     /*!< adjoint operator: transpose (real)
                                               or conjugate transpose (complex)  */
} HLIBPF(matop_t);

/***********************************************
 *
 * types for building cluster trees and
 * block cluster trees
 *
 ***********************************************/

/** partitioning algorithms for binary sparse partitioning */
typedef enum { HLIB_BSP_AUTO,        /*!< automatic choice                        */
               HLIB_BSP_GEOM_MAX,    /*!< geometrically balanced clustering       */
               HLIB_BSP_GEOM_REG,    /*!<   w.r.t. maximal dimension or regularly */
               HLIB_BSP_CARD_MAX,    /*!< cardinality balanced clustering         */
               HLIB_BSP_CARD_REG,    /*!<   w.r.t. maximal dimension or regularly */
               HLIB_BSP_PCA,         /*!< clustering based on principle component */
               HLIB_BSP_GEOM,        /*!< deprecated, now HLIB_BSP_GEOM_MAX       */
               HLIB_BSP_CARD,        /*!< deprecated, now HLIB_BSP_CARD_MAX       */
               HLIB_BSP_REGULAR      /*!< deprecated, now HLIB_BSP_GEOM_REG       */
} HLIBPF(bsp_t);

/** partitioning algorithms for algebraic clustering */
typedef enum { HLIB_ALG_AUTO,        /*!< automatic choice                        */
               HLIB_ALG_BFS,         /*!< BFS graph partitioning                  */
               HLIB_ALG_ML,          /*!< MultiLevel graph partitioning           */
               HLIB_ALG_METIS,       /*!< METIS library                           */
               HLIB_ALG_SCOTCH       /*!< SCOTCH library                          */
} HLIBPF(alg_t);

/** admissibility criteria */
typedef enum { HLIB_ADM_AUTO,        /*!< automatically decide suitable adm. condition  */
               HLIB_ADM_STD_MIN,     /*!< standard admissibility with min. of diameters */
               HLIB_ADM_STD_MAX,     /*!< standard admissibility with max. of diameters */
               HLIB_ADM_WEAK         /*!< weak admissibility                            */
} HLIBPF(adm_t);

/***********************************************
 *
 * accuracy handling
 *
 ***********************************************/

/** different types of accuracy */
typedef enum { HLIB_ACC_FIXED_EPS,   /*!< fixed accuracy for all blocks          */
               HLIB_ACC_FIXED_RANK,  /*!< fixed rank for all blocks              */
               HLIB_ACC_BLOCKED,     /*!< different accuracy for one block level */
               HLIB_ACC_STRING       /*!< arbitrary accuracy description         */
} HLIBPF(acc_type_t);

/*
 * all accuracy types MUST have \c type as first member to ensure
 * consistency in the hlib_acc_t type
 */

/* general accuracy description type (forward decl) */
union HLIBPF(acc_u);

/** accuracy type for fixed accuracy */
typedef struct {
    HLIBPF(acc_type_t)     type;      /* accuracy type */
    HLIBPF(real_t)         eps;       /* single, global accuracy */
} HLIBPF(acc_fixed_eps_t);

/** accuracy type for fixed rank */
typedef struct {
    HLIBPF(acc_type_t)     type;      /* accuracy type */
    unsigned int           rank;      /* single, global rank */
} HLIBPF(acc_fixed_rank_t);

/** accuracy type for block-wise accuracy */
typedef struct {
    HLIBPF(acc_type_t)     type;      /* accuracy type */
    const
    union HLIBPF(acc_u) *  subacc;    /* array of dimension nrows*ncolumns for
                                       * accuracy on next block level (column
                                       * wise ordering)
                                       */
} HLIBPF(acc_blocked_t);

/** general accuracy description type */
typedef union HLIBPF(acc_u) {
    HLIBPF(acc_type_t)        type;   /* accuracy type */

    HLIBPF(acc_fixed_eps_t)   fixed_eps;
    HLIBPF(acc_fixed_rank_t)  fixed_rank;
    HLIBPF(acc_blocked_t)     blocked;
} HLIBPF(acc_t);

/***********************************************
 *
 * types for building matrices
 *
 ***********************************************/

/** different types of low rank approximation */
typedef enum { HLIB_LRAPX_SVD,     /*!< use singular value decomposition             */
               HLIB_LRAPX_ACA,     /*!< use adaptive cross approximation             */
               HLIB_LRAPX_ACAPLUS, /*!< use advanced adaptive cross approximation    */
               HLIB_LRAPX_ACAFULL, /*!< use adaptive cross approximation with        
                                        full pivot search                            */
               HLIB_LRAPX_HCA,     /*!< use hybrid cross approximation               */
               HLIB_LRAPX_ZERO,    /*!< build empty low-rank blocks (nearfield only) */
               HLIB_LRAPX_RANDSVD, /*!< use randomized SVD                           */
               HLIB_LRAPX_RRQR     /*!< uses rank revealing QR                       */
} HLIBPF(lrapx_t);

/*
 * matrix coefficient functions for (complex) dense matrices:
 * evalute the subblock rowidx x colidx of size n x m the matrix
 */
typedef void (* HLIBPF(coeff_t))   ( const size_t           n,
                                     const int *            rowidx,
                                     const size_t           m,
                                     const int *            colidx,
                                     HLIBPF(real_t) *       matrix,
                                     void *                 arg );
typedef void (* HLIBPF(ccoeff_t))  ( const size_t           n,
                                     const int *            rowidx,
                                     const size_t           m,
                                     const int *            colidx,
                                     HLIBPF(complex_t) *    matrix,
                                     void *                 arg );

/***********************************************
 *
 * solving information type
 *
 ***********************************************/

typedef struct {
    unsigned int    converged;       /* 1 if converged, 0 otherwise */
    unsigned int    failed;          /* 1 if failed, 0 otherwise */
    unsigned int    steps;           /* number of iteration steps   */
    HLIBPF(real_t)  res_norm;        /* final norm of residual      */
    HLIBPF(real_t)  conv_rate;       /* average convergence rate    */
} HLIBPF(solve_info_t);

/***********************************************
 *
 * misc. types
 *
 ***********************************************/

/** matrix printing options (can be combined) */
enum { HLIB_MATIO_SVD     = 0x1,   /*!< print singular value decomposition in each block */
       HLIB_MATIO_ENTRY   = 0x2,   /*!< print each entry of matrix                       */
       HLIB_MATIO_PATTERN = 0x4    /*!< print sparsity pattern (non-zero entries)        */
};

/**
 * list of error codes
 */
enum { HLIB_NO_ERROR            = 0,         /*!< no error occured */

       HLIB_ERR_INIT            = 100,       /*!< not initialised */
       HLIB_ERR_LICENSE         = 101,       /*!< invalid license */
       HLIB_ERR_NOT_IMPL        = 102,       /*!< functionality not implemented */
       HLIB_ERR_CONSISTENCY     = 103,       /*!< general consistency error */
       HLIB_ERR_COMM            = 104,       /*!< communication error */
       HLIB_ERR_PERM            = 105,       /*!< permission denied */

       HLIB_ERR_REAL            = 200,       /*!< data is real valued */
       HLIB_ERR_NREAL           = 201,       /*!< data is not real valued */
       HLIB_ERR_COMPLEX         = 202,       /*!< data is complex valued */
       HLIB_ERR_NCOMPLEX        = 203,       /*!< data is not complex valued */
       HLIB_ERR_REAL_CMPLX      = 204,       /*!< invalid mixing of real and complex data */
       HLIB_ERR_DIV_ZERO        = 205,       /*!< division by zero */
       HLIB_ERR_NEG_SQRT        = 206,       /*!< sqrt of negative number */
       HLIB_ERR_INF             = 207,       /*!< infinity occured */
       HLIB_ERR_NAN             = 208,       /*!< not-a-number occured */
       HLIB_ERR_NCONVERGED      = 209,       /*!< iteration did not converge */

       HLIB_ERR_ARG             = 300,       /*!< error with argument */
       HLIB_ERR_MEM             = 301,       /*!< insufficient memory available */
       HLIB_ERR_NULL            = 302,       /*!< null pointer encountered */
       HLIB_ERR_SIZE            = 303,       /*!< size of data incorrect */
       HLIB_ERR_INDEXSET        = 304,       /*!< size of data incorrect */
       HLIB_ERR_DIM             = 305,       /*!< invalid or incompatible dimension */
       HLIB_ERR_ARR_BOUND       = 306,       /*!< out-of-bound error in array */
       HLIB_ERR_DIAG_ENTRY      = 307,       /*!< entry is not on diagonal */

       HLIB_ERR_COORD_INVALID   = 400,       /*!< invalid coordinates */

       HLIB_ERR_CT_INVALID      = 500,       /*!< invalid cluster tree */
       HLIB_ERR_CT_TYPE         = 501,       /*!< wrong type of cluster tree */
       HLIB_ERR_CT_STRUCT       = 502,       /*!< invalid structure of cluster tree */
       HLIB_ERR_CT_INCOMP       = 503,       /*!< given cluster trees are incompatible */
       HLIB_ERR_CT_SPARSE       = 504,       /*!< missing sparse matrix for given cluster tree */

       HLIB_ERR_BCT_INVALID     = 600,       /*!< invalid block cluster tree */
       HLIB_ERR_BCT_STRUCT      = 601,       /*!< invalid block cluster tree structure */

       HLIB_ERR_VEC_INVALID     = 700,       /*!< invalid vector */
       HLIB_ERR_VEC_TYPE        = 701,       /*!< wrong vector type */
       HLIB_ERR_VEC_STRUCT      = 702,       /*!< invalid vector structure */
       HLIB_ERR_VEC_SIZE        = 703,       /*!< invalid size of vector */
       HLIB_ERR_VEC_INCOMP      = 704,       /*!< vector with incompatible dimension */
       HLIB_ERR_VEC_NSCALAR     = 705,       /*!< vector is not a scalar vector */

       HLIB_ERR_MAT_TYPE        = 800,       /*!< invalid matrix type */
       HLIB_ERR_MAT_STRUCT      = 801,       /*!< invalid structure of matrix */
       HLIB_ERR_MAT_SIZE        = 802,       /*!< invalid size of matrix */
       HLIB_ERR_MAT_SINGULAR    = 803,       /*!< singular matrix detected */
       HLIB_ERR_MAT_NSPARSE     = 804,       /*!< matrix not a sparse matrix */
       HLIB_ERR_MAT_NDENSE      = 805,       /*!< matrix not a dense matrix */
       HLIB_ERR_MAT_NHMAT       = 806,       /*!< matrix not an H-matrix */
       HLIB_ERR_MAT_INCOMP_TYPE = 807,       /*!< matrices with incompatible type */
       HLIB_ERR_MAT_INCOMP_CT   = 808,       /*!< matrices with incompatible cluster tree */
       HLIB_ERR_MAT_INVALID     = 809,       /*!< invalid matrix */
       HLIB_ERR_MAT_NSYM        = 810,       /*!< matrix not symmetric */
       HLIB_ERR_MAT_NHERM       = 811,       /*!< matrix not hermitian */
       HLIB_ERR_MAT_NPOSDEF     = 812,       /*!< matrix not positiv definite */

       HLIB_ERR_FMT_UNKNOWN     = 900,       /*!< error while parsing HLIBpro format */
       HLIB_ERR_FMT_HFORMAT     = 901,       /*!< error while parsing HLIBpro format */
       HLIB_ERR_FMT_SAMG        = 902,       /*!< error while parsing SAMG format */
       HLIB_ERR_FMT_MATLAB      = 903,       /*!< error while parsing Matlab format */
       HLIB_ERR_FMT_PLTMG       = 904,       /*!< error while parsing PLTMG format */
       HLIB_ERR_FMT_HB          = 905,       /*!< error while parsing Harwell Boeing format */
       HLIB_ERR_FMT_MTX         = 906,       /*!< error while parsing Matrix Market format */
       HLIB_ERR_FMT_PLY         = 907,       /*!< error while parsing Ply format */

       HLIB_ERR_GRID_FORMAT     = 1000,      /*!< invalid format of grid file */
       HLIB_ERR_GRID_DATA       = 1001,      /*!< invalid data in grid file */

       HLIB_ERR_FOPEN           = 1100,      /*!< could not open file */
       HLIB_ERR_FCLOSE          = 1101,      /*!< could not close file */
       HLIB_ERR_FWRITE          = 1102,      /*!< could not write to file */
       HLIB_ERR_FREAD           = 1103,      /*!< could not read from file */
       HLIB_ERR_FSEEK           = 1104,      /*!< could not seek in file */
       HLIB_ERR_FNEXISTS        = 1105,      /*!< file does not exists */

       HLIB_ERR_BS_SIZE         = 1200,      /*!< size of bytestream too small */
       HLIB_ERR_BS_WRITE        = 1201,      /*!< error while writing to bytestream */
       HLIB_ERR_BS_READ         = 1202,      /*!< error while reading from bytestream */
       HLIB_ERR_BS_TYPE         = 1203,      /*!< type error in bytestream */
       HLIB_ERR_BS_DATA         = 1204,      /*!< general data error in bytestream */

       HLIB_ERR_NOZLIB          = 1300,      /*!< no zlib support compiled in */
       HLIB_ERR_ZLIB_UNZIP      = 1301,      /*!< error during zlib uncompression */
       HLIB_ERR_NOMETIS         = 1302,      /*!< no METIS support compiled in */
       HLIB_ERR_NOSCOTCH        = 1303,      /*!< no Scotch support compiled in */
       HLIB_ERR_SCOTCH          = 1304,      /*!< error in call to Scotch function */
       HLIB_ERR_NOCHACO         = 1305,      /*!< no Chaco support compiled in */
       HLIB_ERR_NOLIBGRAPH      = 1306,      /*!< no libGraph support compiled in */
       HLIB_ERR_NOFFTW3         = 1307,      /*!< no FFTW3 support compiled in */
       HLIB_ERR_NOCAIRO         = 1308,      /*!< no Cairo support compiled in */
       HLIB_ERR_NOHDF5          = 1309,      /*!< no HDF5 support compiled in */

       HLIB_ERR_MPI             = 1400,      /*!< error in call to MPI function */

       HLIB_ERR_SOLVER_INVALID  = 1500,      /*!< invalid solver */
       HLIB_ERR_LRAPX_INVALID   = 1501,      /*!< invalid low-rank approximation type */
       HLIB_ERR_GRID_INVALID    = 1502,      /*!< invalid grid */
       HLIB_ERR_FNSPACE_INVALID = 1503       /*!< invalid function space */
};

/**
 * error/warning handler function types
 */
typedef void (* HLIBPF(errorfn_t))      ( const int         errcode,
                                          const char *      errmsg );

/**
 * callback function for progress information
 */
typedef void (* HLIBPF(progressfn_t))   ( const double *    values,   /* min/max/current value of progress */
                                          int *             cancel,   /* signal cancelation */
                                          void *            arg );    /* user data for callback function */

/** fields in value array progress bar callback function */
enum { HLIB_PROGRESS_MIN = 0,   /*!< minimal value of progress bar */
       HLIB_PROGRESS_MAX = 1,   /*!< maximal value of progress bar */
       HLIB_PROGRESS_VAL = 2    /*!< current value of progress bar */
};

/**************************************************************
 **************************************************************
 **
 ** general HLIBpro functions
 **
 **************************************************************
 **************************************************************/

/************************************************
 **
 ** initialisation and finalisation
 **
 ************************************************/

/**
 * initialise HLIBpro
 */
HLIB_FNDECL void
HLIBPF(init)                       ( int *                  info );

/**
 * finalise HLIBpro
 */
HLIB_FNDECL void
HLIBPF(done)                       ( int *                  info );

/**
 * return 1 if HLIBpro is initialised and 0 otherwise
 */
HLIB_FNDECL int
HLIBPF(is_init)                    ();

/************************************************
 **
 ** error handling
 **
 ************************************************/

/** copy description of current error to buffer */
HLIB_FNDECL void
HLIBPF(error_desc)                 ( char *                 desc,   /* buffer to copy description to */
                                     const size_t           size ); /* size of \a desc in bytes      */



/***********************************************************//**
 ***************************************************************
 ** 
 ** @{ \name Coordinates
 **
 ***************************************************************
 ***************************************************************/

/**
 * import given coordinates into HLIBpro
 * - data is directly referenced, not copied, i.e., any later changes to
 *   \a coord will also change coordinate set,
 * - optional periodicity of the coordinates defined by \a period,
 *   e.g. vector containing stride in all spatial dimensions or
 *   NULL if not defined
 */
HLIB_FNDECL HLIBPF(coord_t)
HLIBPF(coord_import)               ( const size_t           n,      /* number of indices                */
                                     const unsigned int     dim,    /* spatial dimension of coordinates */
                                     double **              coord,  /* coordinates of indices           */
                                     const double *         period, /* periodicity of coordinates       */
                                     int *                  info ); /* to return error code             */

/****************************************************
 *
 * data access
 *
 ****************************************************/

/** return number of coordinates in given coordinate set \a coord */
HLIB_FNDECL size_t
HLIBPF(coord_ncoord)               ( const HLIBPF(coord_t)  coord,
                                     int *                  info );

/** return dimension of coordinates in given coordinate set \a coord */
HLIB_FNDECL unsigned int
HLIBPF(coord_dim)                  ( const HLIBPF(coord_t)  coord,
                                     int *                  info );

/**
 * return pointer too coordinate \a i in given coordinate set \a coord
 * - the returned value points to an array of size dim(coord)
 * - to pointer is no copy, so all changes to the data are direct
 *   changes to the coordinates
 */
HLIB_FNDECL double *
HLIBPF(coord_get)                  ( const HLIBPF(coord_t)  coord,
                                     const size_t           i,
                                     int *                  info );

/**
 * return pointer too coordinate of \a i'th minimal bounding box
 */
HLIB_FNDECL double *
HLIBPF(coord_bbmin)                ( const HLIBPF(coord_t)  coord,
                                     const size_t           i,
                                     int *                  info );

/**
 * return pointer too coordinate of \a i'th maximal bounding box
 */
HLIB_FNDECL double *
HLIBPF(coord_bbmax)                ( const HLIBPF(coord_t)  coord,
                                     const size_t           i,
                                     int *                  info );

/**
 * return 1, if coordinate set has bounding box info and 0 otherwise
 */
HLIB_FNDECL int
HLIBPF(coord_has_bbox)             ( const HLIBPF(coord_t)  coord,
                                     int *                  info );

/** set perdiodicity of coordinates to \a period
 * - \a period must have the same dimensions as the coordinates
 */
HLIB_FNDECL void
HLIBPF(coord_set_period)           ( HLIBPF(coord_t)        coord,
                                     const double           period[],
                                     int *                  info );

/** store periodicity vector in \a period
 * - \a period must have the same dimensions as the coordinates
 */
HLIB_FNDECL void
HLIBPF(coord_get_period)           ( HLIBPF(coord_t)        coord,
                                     double                 period[],
                                     int *                  info );

/****************************************************
 *
 * management
 *
 ****************************************************/

/** free resources coupled with coordinates
 * - if "HLIBPF(coord_import)" was used to create \a coord,
 *   the memory occupied by coordinate array will NOT be freed
 */
HLIB_FNDECL void
HLIBPF(coord_free)                 ( HLIBPF(coord_t)        coord,
                                     int *                  info );

/** return size of memory in bytes used by coordinates */
HLIB_FNDECL size_t
HLIBPF(coord_bytesize)             ( const HLIBPF(coord_t)  coord,
                                     int *                  info );

/****************************************************
 *
 * input/output
 *
 ****************************************************/

/** print coordinates in VRML format to file \a filename */
HLIB_FNDECL void
HLIBPF(coord_print_vrml)           ( const HLIBPF(coord_t)  coord,
                                     const char *           filename,
                                     int *                  info );

/** print coordinates in VTK format to file \a filename */
HLIB_FNDECL void
HLIBPF(coord_print_vtk)            ( const HLIBPF(coord_t)  coord,
                                     const char *           filename,
                                     int *                  info );

/** @} */

/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Admissibility condition
 **
 ***************************************************************
 ***************************************************************/

/**
 * create admissibility condition based on geometrical data
 */
HLIB_FNDECL HLIBPF(admcond_t)
HLIBPF(admcond_geom)               ( const HLIBPF(adm_t)      crit,
                                     const double             eta,
                                     int *                    info );
    
/**
 * create admissibility condition based on geometrical data
 * with additional periodicity in geometry defined by vector \a period
 * - \a period must have same dimension \a dim as coordinates used for 
 *   constructing cluster trees
 */
HLIB_FNDECL HLIBPF(admcond_t)
HLIBPF(admcond_geom_period)        ( const HLIBPF(adm_t)      crit,
                                     const double             eta,
                                     const double *           period,
                                     const unsigned int       dim,
                                     int *                    info );
    
/**
 * create admissibility condition for high/low frequency case by limiting
 * the number of wavelengths (nwaves) per cluster
 */
HLIB_FNDECL HLIBPF(admcond_t)
HLIBPF(admcond_geom_hilo)          ( const HLIBPF(complex_t)  kappa,
                                     const unsigned int       nwaves,
                                     const double             eta,
                                     int *                    info );
    
/**
 * create admissibility condition based on algebraic connectivity in \a S
 * - \a row_perm and \a col_perm define the (optional) mapping between 
 *   internal to external ordering for row and column index sets, respectively
 *   (see hlib_ct_perm_i2e); set to NULL if not present
 */
HLIB_FNDECL HLIBPF(admcond_t)
HLIBPF(admcond_alg)                ( const HLIBPF(adm_t)          crit,
                                     const double                 eta,
                                     const HLIBPF(matrix_t)       S,
                                     const HLIBPF(permutation_t)  row_perm,
                                     const HLIBPF(permutation_t)  col_perm,
                                     int *                        info );
    
/**
 * free admissibility condition
 */
HLIB_FNDECL void
HLIBPF(admcond_free)               ( const HLIBPF(admcond_t)  ac,
                                     int *                    info );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Clusters
 **
 ** Nodes and sub trees in a cluster tree.
 **
 ***************************************************************
 ***************************************************************/

/****************************************************
 *
 * index set functions
 *
 ****************************************************/

/** return first index in index set */
HLIB_FNDECL int
HLIBPF(cl_first)                   ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/** return last index in index set */
HLIB_FNDECL int
HLIBPF(cl_last)                    ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/** return size of index set */
HLIB_FNDECL size_t
HLIBPF(cl_size)                    ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/****************************************************
 *
 * tree functions
 *
 ****************************************************/

/** return number of sons of cluster */
HLIB_FNDECL size_t
HLIBPF(cl_nsons)                   ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/** return \a i'th son of cluster \a cl */
HLIB_FNDECL HLIBPF(cluster_t)
HLIBPF(cl_son)                     ( const HLIBPF(cluster_t)  cl,
                                     const unsigned int       i,
                                     int *                    info );

/** return 1 if cluster is leaf and 0 otherwise */
HLIB_FNDECL int
HLIBPF(cl_is_leaf)                 ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/** return number of nodes in cluster tree */
HLIB_FNDECL size_t
HLIBPF(cl_nnodes)                  ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/** return depth of cluster tree */
HLIB_FNDECL size_t
HLIBPF(cl_depth)                   ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/****************************************************
 *
 * management
 *
 ****************************************************/

/** return copy of sub tree defined by \a cl */
HLIB_FNDECL HLIBPF(cluster_t)
HLIBPF(cl_copy)                    ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/** free resources coupled with cluster cl */
HLIB_FNDECL void
HLIBPF(cl_free)                    ( HLIBPF(cluster_t)        cl,
                                     int *                    info );

/** return size of memory in bytes used by cluster */
HLIB_FNDECL size_t
HLIBPF(cl_bytesize)                ( const HLIBPF(cluster_t)  cl,
                                     int *                    info );

/****************************************************
 *
 * input/output
 *
 ****************************************************/

/** print cluster tree in PostScript format to file \a filename */
HLIB_FNDECL void
HLIBPF(cl_print_ps)                ( const HLIBPF(cluster_t)  cl,
                                     const char *             filename,
                                     int *                    info );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Cluster Trees
 **
 ** Root of a cluster tree with stored index permutations.
 **
 ***************************************************************
 ***************************************************************/

/**
 * build cluster tree using binary space partitioning
 * - optional periodicity of the coordinates defined by \a period,
 *   e.g. vector containing stride in all spatial dimensions or
 *   NULL if not defined
 */
HLIB_FNDECL HLIBPF(clustertree_t)
HLIBPF(clt_build_bsp)              ( const HLIBPF(coord_t)    coord,     /*!< coordinates of indices  */
                                     const HLIBPF(bsp_t)      bsptype,   /*!< partitioning type       */
                                     const unsigned int       nmin,      /*!< minimal cluster size    */
                                     int *                    info );

/**
 * build cluster tree using binary space partitioning and
 * nested dissection based on given sparse matrix S
 * - for \a period see HLIBPF(clt_build_bsp)
 */
HLIB_FNDECL HLIBPF(clustertree_t)
HLIBPF(clt_build_bsp_nd)           ( const HLIBPF(coord_t)    coord,     /*!< coordinates of indices       */
                                     const HLIBPF(matrix_t)   S,         /*!< sparse mat. for connectivity */
                                     const HLIBPF(bsp_t)      bsptype,   /*!< partitioning type            */
                                     const unsigned int       nmin,      /*!< minimal cluster size         */
                                     int *                    info );

/**
 * build cluster tree using binary space partitioning but
 * use a predefined partition for first step, e.g. first
 * divide indices into groups defined by \a partition and
 * apply BSP on these sets
 * - \a partition must have size "size(coord)" and contain ids ∈ [0,k-1],
 *   where k is the total number of partitions
 * - index "i" will go into group \a partition[i]
 * - if \a eq_depth ≠ 0, the depth of subsequent sub trees
 *   (for each partition) is adjusted according to depth of
 *   previous sub trees
 */
HLIB_FNDECL HLIBPF(clustertree_t)
HLIBPF(clt_build_bsp_part)         ( const HLIBPF(coord_t)    coord,     /*!< coordinates of indices      */
                                     const unsigned int *     partition, /*!< predifined partition        */
                                     const int                eq_depth,  /*!< equalise depths of subtrees */
                                     const HLIBPF(bsp_t)      bsptype,   /*!< partitioning type           */
                                     const unsigned int       nmin,      /*!< minimal cluster size        */
                                     int *                    info );

/**
 * build clustertree via algebraic clustering based on
 * given sparse matrix \a S and partitioning type \a algtype.
 */
HLIB_FNDECL HLIBPF(clustertree_t)
HLIBPF(clt_build_alg)              ( const HLIBPF(matrix_t)   S,
                                     const HLIBPF(alg_t)      algtype,
                                     const unsigned int       nmin,     
                                     int *                    info );

/**
 * build clustertree via algebraic clustering with nested
 * dissection based on given sparse matrix \a S and partitioning
 * type \a algtype.
 */
HLIB_FNDECL HLIBPF(clustertree_t)
HLIBPF(clt_build_alg_nd)           ( const HLIBPF(matrix_t)   S,
                                     const HLIBPF(alg_t)      algtype,
                                     const unsigned int       nmin,     
                                     int *                    info );

/**
 * build clustertree via algebraic clustering based on matrix \a S 
 * and use a predefined partition for first step
 * - \a partition must have size "rows(S)" and contain ids ∈ [0,k-1],
 *   where k is the total number of partitions
 * - index "i" will go into group \a partition[i]
 */
HLIB_FNDECL HLIBPF(clustertree_t)
HLIBPF(clt_build_alg_part)         ( const HLIBPF(matrix_t)   S,
                                     const unsigned int *     partition,
                                     const HLIBPF(alg_t)      algtype,
                                     const unsigned int       nmin,     
                                     int *                    info );

/****************************************************
 *
 * management
 *
 ****************************************************/

/** return root node in cluster tree */
HLIB_FNDECL HLIBPF(cluster_t)
HLIBPF(clt_root)                   ( HLIBPF(clustertree_t)        ct,
                                     int *                        info );

/** return mapping from internal (in cluster tree) to external ordering
 *  - returned object is reference to internal permutation: do NOT free
 */
HLIB_FNDECL HLIBPF(permutation_t)
HLIBPF(clt_perm_i2e)               ( HLIBPF(clustertree_t)        ct,
                                     int *                        info );

/** return mapping from external to internal (in cluster tree) ordering
 *  - returned object is reference to internal permutation: do NOT free
 */
HLIB_FNDECL HLIBPF(permutation_t)
HLIBPF(clt_perm_e2i)               ( HLIBPF(clustertree_t)        ct,
                                     int *                        info );

/** free resources coupled with cluster tree ct */
HLIB_FNDECL void
HLIBPF(clt_free)                   ( HLIBPF(clustertree_t)        ct,
                                     int *                        info );

/** return size of memory in bytes used by cluster */
HLIB_FNDECL size_t
HLIBPF(clt_bytesize)               ( const HLIBPF(clustertree_t)  ct,
                                     int *                        info );

/** return number of nodes in cluster tree */
HLIB_FNDECL size_t
HLIBPF(clt_nnodes)                 ( const HLIBPF(clustertree_t)  ct,
                                     int *                        info );

/** return depth of cluster tree */
HLIB_FNDECL size_t
HLIBPF(clt_depth)                  ( const HLIBPF(clustertree_t)  ct,
                                     int *                        info );

/****************************************************
 *
 * input/output
 *
 ****************************************************/

/** print cluster tree in PostScript format to file \a filename */
HLIB_FNDECL void
HLIBPF(clt_print_ps)                ( const HLIBPF(clustertree_t)  ct,
                                     const char *                 filename,
                                     int *                        info );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Block Clusters
 **
 ** Nodes and sub trees of block cluster trees.
 **
 ***************************************************************
 ***************************************************************/

/****************************************************
 *
 * tree functions
 *
 ****************************************************/

/** return parent block cluster of block cluster \a bc */
HLIB_FNDECL HLIBPF(blockcluster_t)
HLIBPF(bc_parent)                 ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/** return total number of sons of block cluster \a bc */
HLIB_FNDECL size_t
HLIBPF(bc_nsons)                  ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/** return son (\a i,\a j) of block cluster \a bc
 * - the indices \a i and \a j are with respect to the numbering in
 *   the row and column cluster trees
 */
HLIB_FNDECL HLIBPF(blockcluster_t)
HLIBPF(bc_son)                    ( const HLIBPF(blockcluster_t)  bc,
                                    const unsigned int            i,
                                    const unsigned int            j,
                                    int *                         info );

/** return 1 if block cluster is leaf and 0 otherwise */
HLIB_FNDECL int
HLIBPF(bc_is_leaf)                ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/** return 1 if block cluster is admissible and 0 otherwise */
HLIB_FNDECL int
HLIBPF(bc_is_adm)                 ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/** return number of nodes in block cluster */
HLIB_FNDECL size_t
HLIBPF(bc_nnodes)                 ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/** return depth of block cluster */
HLIB_FNDECL size_t
HLIBPF(bc_depth)                  ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/****************************************************
 *
 * management
 *
 ****************************************************/

/** return row cluster of given block cluster */
HLIB_FNDECL HLIBPF(cluster_t)
HLIBPF(bc_rowcl)                  ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/** return column cluster of given block cluster */
HLIB_FNDECL HLIBPF(cluster_t)
HLIBPF(bc_colcl)                  ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/** return copy of sub tree defined by \a bc */
HLIB_FNDECL HLIBPF(blockcluster_t)
HLIBPF(bc_copy)                   ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/** free resources coupled with block cluster tree bct */
HLIB_FNDECL void
HLIBPF(bc_free)                   ( HLIBPF(blockcluster_t)        bc,
                                    int *                         info );

/** return size of memory in bytes used by block cluster */
HLIB_FNDECL size_t
HLIBPF(bc_bytesize)               ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/****************************************************
 *
 * misc. functions
 *
 ****************************************************/

/** compute sparsity constant of block cluster tree */
HLIB_FNDECL int
HLIBPF(bc_csp)                    ( const HLIBPF(blockcluster_t)  bc,
                                    int *                         info );

/****************************************************
 *
 * input/output
 *
 ****************************************************/

/** print block cluster tree in PostScript format to file \a filename */
HLIB_FNDECL void
HLIBPF(bc_print_ps)               ( const HLIBPF(blockcluster_t)  bc,
                                    const char *                  filename,
                                    int *                         info );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Block Cluster Trees
 **
 ** Root of block cluster trees with access to index permutations.
 **
 ***************************************************************
 ***************************************************************/

/****************************************************
 *
 * building
 *
 ****************************************************/

/**
 * build block cluster tree over given row and column clusters
 * using admissibility condition \a ac
 */
HLIB_FNDECL HLIBPF(blockclustertree_t)
HLIBPF(bct_build)                  ( const HLIBPF(clustertree_t)       rowct,
                                     const HLIBPF(clustertree_t)       colct,
                                     const HLIBPF(admcond_t)           ac,
                                     int *                             info );

/****************************************************
 *
 * tree functions
 *
 ****************************************************/

/** return number of nodes in block cluster tree */
HLIB_FNDECL size_t
HLIBPF(bct_nnodes)                 ( const HLIBPF(blockclustertree_t)  bct,
                                     int *                             info );

/** return depth of block cluster tree */
HLIB_FNDECL size_t
HLIBPF(bct_depth)                  ( const HLIBPF(blockclustertree_t)  bct,
                                     int *                             info );

/****************************************************
 *
 * management
 *
 ****************************************************/

/** return row cluster tree of given block cluster tree */
HLIB_FNDECL HLIBPF(clustertree_t)
HLIBPF(bct_rowct)                  ( const HLIBPF(blockclustertree_t)  bct,
                                     int *                             info );

/** return column cluster tree of given block cluster tree */
HLIB_FNDECL HLIBPF(clustertree_t)
HLIBPF(bct_colct)                  ( const HLIBPF(blockclustertree_t)  bct,
                                     int *                             info );

/** free resources coupled with cluster tree ct */
HLIB_FNDECL void
HLIBPF(bct_free)                   ( HLIBPF(blockclustertree_t)        bct,
                                     int *                             info );

/** return size of memory in bytes used by block cluster */
HLIB_FNDECL size_t
HLIBPF(bct_bytesize)               ( const HLIBPF(blockclustertree_t)  bct,
                                     int *                             info );

/****************************************************
 *
 * distribution functions
 *
 ****************************************************/

/** distribute block cluster tree onto \a p processors based on
 *  block-wise partition with unique processor per block
 *  - assuming compatible structure of \a bct
 *  - \a min_lvl controls minimal level for blocks in bct to
 *    schedule
 *  - if \a symmetric ≠ 0, the scheduling is performed for the
 *    lower left part and mirrored to the upper right part
 *    (for symmetric/hermitian matrices)
 */
HLIB_FNDECL void
HLIBPF(bct_distribute_block)       ( HLIBPF(blockclustertree_t)        bct,
                                     const unsigned int                p,
                                     const unsigned int                min_lvl,
                                     const int                         symmetric,
                                     int *                             info );
 
/****************************************************
 *
 * misc. functions
 *
 ****************************************************/

/** compute sparsity constant of block cluster tree */
HLIB_FNDECL int
HLIBPF(bct_csp)                    ( const HLIBPF(blockclustertree_t)  bct,
                                     int *                             info );

/****************************************************
 *
 * input/output
 *
 ****************************************************/

/** print block cluster tree in PostScript format to file \a filename */
HLIB_FNDECL void
HLIBPF(bct_print_ps)               ( const HLIBPF(blockclustertree_t)  bct,
                                     const char *                      filename,
                                     int *                             info );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Vectors
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 **
 ** general vector functions
 **
 ************************************************/

/** return size of vectors */
HLIB_FNDECL size_t
HLIBPF(vector_size)                ( const HLIBPF(vector_t)     v,
                                     int *                      info );

/** copy vectormatrix */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(vector_copy)                ( const HLIBPF(vector_t)     v,
                                     int *                      info );

/** return size of matrix in bytes */
HLIB_FNDECL size_t
HLIBPF(vector_bytesize)            ( const HLIBPF(vector_t)     v,
                                     int *                      info );

/** free resources coupled with vector v */
HLIB_FNDECL void
HLIBPF(vector_free)                ( HLIBPF(vector_t)           v,
                                     int *                      info );

/************************************************
 **
 ** access vector data
 **
 ************************************************/

/** get single entry at position \a i in vector \a x */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(vector_entry_get)           ( const HLIBPF(vector_t)     x,
                                     const size_t               i,
                                     int *                      info );

/** get single entry at position \a i in vector \a x */
HLIB_FNDECL HLIBPF(complex_t)
HLIBPF(vector_centry_get)          ( const HLIBPF(vector_t)     x,
                                     const size_t               i,
                                     int *                      info );

/** set single entry at position \a i in vector \a x */
HLIB_FNDECL void
HLIBPF(vector_entry_set)           ( const HLIBPF(vector_t)     x,
                                     const size_t               i,
                                     const HLIBPF(real_t)       f,
                                     int *                      info );

/** set single entry at position \a i in vector \a x */
HLIB_FNDECL
void
HLIBPF(vector_centry_set)          ( const HLIBPF(vector_t)     x,
                                     const size_t               i,
                                     const HLIBPF(complex_t)    f,
                                     int *                      info );

/************************************************
 **
 ** vector import/building
 **
 ************************************************/

/**
 * create scalar vectors of size \a size initialised with 0
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(vector_build)               ( const size_t               size,
                                     int *                      info );

/**
 * create scalar vectors of size \a size initialised with 0
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(vector_cbuild)              ( const size_t               size,
                                     int *                      info );

/**
 * copy data from C array \a arr into vector \a v
 * - arr has to be of size hlib_vector_size( v )
 */
HLIB_FNDECL void
HLIBPF(vector_import)              ( HLIBPF(vector_t)           v,
                                     const HLIBPF(real_t) *     arr,
                                     int *                      info );

/**
 * copy data from C array \a arr into vector \a v
 * - arr has to be of size hlib_vector_size( v )
 */
HLIB_FNDECL void
HLIBPF(vector_cimport)             ( HLIBPF(vector_t)           v,
                                     const HLIBPF(complex_t) *  arr,
                                     int *                      info );

/**
 * copy data from vector \a v into given C arrays
 * - arr has to be of size hlib_vector_size( v )
 */
HLIB_FNDECL void
HLIBPF(vector_export)              ( const HLIBPF(vector_t)     v,
                                     HLIBPF(real_t) *           arr,
                                     int *                      info );

/**
 * copy data from vector \a v into given C arrays
 * - arr has to be of size hlib_vector_size( v )
 */
HLIB_FNDECL void
HLIBPF(vector_cexport)             ( const HLIBPF(vector_t)     v,
                                     HLIBPF(complex_t) *        arr,
                                     int *                      info );

/************************************************
 **
 ** algebraic vector functions
 **
 ************************************************/

/** fill vector \a x with constant value \a f */
HLIB_FNDECL void
HLIBPF(vector_fill)                ( HLIBPF(vector_t)           x,
                                     const HLIBPF(real_t)       f,
                                     int *                      info );

/** fill vector \a x with constant value \a f */
HLIB_FNDECL void
HLIBPF(vector_cfill)               ( HLIBPF(vector_t)           x,
                                     const HLIBPF(complex_t)    f,
                                     int *                      info );

/** fill vector \a x with random values */
HLIB_FNDECL void
HLIBPF(vector_fill_rand)           ( HLIBPF(vector_t)           x,
                                     int *                      info );

/** copy content of vector \a x to vector \a y
 *  \a x and \a y have to be of equal type */
HLIB_FNDECL void
HLIBPF(vector_assign)              ( HLIBPF(vector_t)           y,
                                     const HLIBPF(vector_t)     x,
                                     int *                      info );

/** scale vector \a x by \a f */
HLIB_FNDECL void
HLIBPF(vector_scale)               ( HLIBPF(vector_t)           x,
                                     const HLIBPF(real_t)       f,
                                     int *                      info );

/** scale vector \a x by \a f */
HLIB_FNDECL void
HLIBPF(vector_cscale)              ( HLIBPF(vector_t)           x,
                                     const HLIBPF(complex_t)    f,
                                     int *                      info );

/** compute y := y + a x  */
HLIB_FNDECL void
HLIBPF(vector_axpy)                ( const HLIBPF(real_t)       alpha,
                                     const HLIBPF(vector_t)     x,
                                     HLIBPF(vector_t)           y,
                                     int *                      info );

/** compute y := y + a x  */
HLIB_FNDECL void
HLIBPF(vector_caxpy)               ( const HLIBPF(complex_t)    alpha,
                                     const HLIBPF(vector_t)     x,
                                     HLIBPF(vector_t)           y,
                                     int *                      info );

/** compute dot product <\a x,\a y> = \a x^H * \a y */
HLIB_FNDECL HLIBPF(complex_t)
HLIBPF(vector_dot)                 ( const HLIBPF(vector_t)     x,
                                     const HLIBPF(vector_t)     y,
                                     int *                      info );

/** compute euklidean norm of \a x */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(vector_norm2)               ( const HLIBPF(vector_t)     x,
                                     int *                      info );

/** compute infinity norm of \a x */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(vector_norm_inf)            ( const HLIBPF(vector_t)     x,
                                     int *                      info );

/************************************************
 **
 ** misc vector functions
 **
 ************************************************/

/** return vector containing only real parts *
 *  of the entries of the given vector
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(vector_restrict_re)         ( const HLIBPF(vector_t)     v,
                                     int *                      info );

/** return vector containing only imaginary parts *
 *  of the entries of the given vector
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(vector_restrict_im)         ( const HLIBPF(vector_t)     v,
                                     int *                      info );

/** reorder vector entries with given permutation */
HLIB_FNDECL void
HLIBPF(vector_permute)             ( HLIBPF(vector_t)             v,
                                     const HLIBPF(permutation_t)  perm,
                                     int *                        info );

/************************************************
 **
 ** optional vector FFT
 **
 ************************************************/

/** compute FFT of vector \a v inplace */
HLIB_FNDECL void
HLIBPF(vector_fft)                 ( const HLIBPF(vector_t)     v,
                                     int *                      info );

/** compute inverse FFT of vector \a v inplace */
HLIB_FNDECL void
HLIBPF(vector_ifft)                ( const HLIBPF(vector_t)     v,
                                     int *                      info );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Matrices
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 **
 ** linear operator functions
 **
 ************************************************/

/** free resources coupled with linear operator A */
HLIB_FNDECL void
HLIBPF(linearoperator_free)          ( HLIBPF(linearoperator_t)        A,
                                       int *                           info );

/** create vectors from the range of the linear operator \a A */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(linearoperator_range_vector)  ( const HLIBPF(linearoperator_t)  A,
                                       int *                           info );

/** create vectors from the domain of the linear operator \a A */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(linearoperator_domain_vector) ( const HLIBPF(linearoperator_t)  A,
                                       int *                           info );

/** represents permuted operator \f$P \cdot A \cdot R\f$ */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(perm_linearoperator)          ( const HLIBPF(permutation_t)     P,
                                       const HLIBPF(linearoperator_t)  A,
                                       const HLIBPF(permutation_t)     R,
                                       int *                           info );

/** computes \f$ y := op(A) \cdot x \f$ */
HLIB_FNDECL void
HLIBPF(linearoperator_apply)         ( const HLIBPF(linearoperator_t)  A,
                                       const HLIBPF(vector_t)          x,
                                       HLIBPF(vector_t)                y,
                                       const HLIBPF(matop_t)           matop,
                                       int *                           info );

/** computes \f$ y := y + \alpha op(A) \cdot x \f$ */
HLIB_FNDECL void
HLIBPF(linearoperator_apply_add)     ( const HLIBPF(real_t)            alpha,
                                       const HLIBPF(linearoperator_t)  A,
                                       const HLIBPF(vector_t)          x,
                                       HLIBPF(vector_t)                y,
                                       const HLIBPF(matop_t)           matop,
                                       int *                           info );

/** computes \f$ y := y + \alpha op(A) \cdot x \f$ */
HLIB_FNDECL void
HLIBPF(linearoperator_capply_add)    ( const HLIBPF(complex_t)         alpha,
                                       const HLIBPF(linearoperator_t)  A,
                                       const HLIBPF(vector_t)          x,
                                       HLIBPF(vector_t)                y,
                                       const HLIBPF(matop_t)           matop,
                                       int *                           info );

/************************************************
 **
 ** general matrix functions
 **
 ************************************************/

/** return number of rows of matrix */
HLIB_FNDECL size_t
HLIBPF(matrix_rows)                ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/** return number of columns of matrix */
HLIB_FNDECL size_t
HLIBPF(matrix_cols)                ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/** return copy of matrix \a A */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_copy)                ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/** return copy of matrix \a A with accuracy \a acc;
 * if \a coarsen != 0, copy is coarsend */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_copy_acc)            ( const HLIBPF(matrix_t)        A,
                                     const HLIBPF(acc_t)           acc,
                                     const int                     coarsen,
                                     int *                         info );

/** copy matrix \a A to \a B */
HLIB_FNDECL void
HLIBPF(matrix_copyto)              ( const HLIBPF(matrix_t)        A,
                                     HLIBPF(matrix_t)              B,
                                     int *                         info );

/** copy matrix \a A to \a B with accuracy \a acc;
 * if \a coarsen != 0, \a B is coarsend */
HLIB_FNDECL void
HLIBPF(matrix_copyto_acc)          ( const HLIBPF(matrix_t)        A,
                                     HLIBPF(matrix_t)              B,
                                     const HLIBPF(acc_t)           acc,
                                     const int                     coarsen,
                                     int *                         info );

/** copy blockdiagonal part of first \a lvl levels or blocks of size \a blocksize
 *  of matrix */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_copy_blockdiag)      ( const HLIBPF(matrix_t)        A,
                                     const unsigned int            lvl,
                                     const size_t                  blocksize,
                                     int *                         info );

/** copy blockdiagonal part of first \a lvl levels or blocks of size \a blocksize
 *  of matrix with given accuracy \a acc */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_copy_blockdiag_acc)  ( const HLIBPF(matrix_t)        A,
                                     const HLIBPF(acc_t)           acc,
                                     const unsigned int            lvl,
                                     const size_t                  blocksize,
                                     int *                         info );

/** copy nearfield part of matrix
 * - if \a without_farfield is non-zero, farfield matrices (e.g., low-rank) will 
 *   be skipped, otherwise set to rank 0 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_copy_nearfield)      ( const HLIBPF(matrix_t)        A,
                                     const int                     without_farfield,
                                     int *                         info );

/** restrict given matrix to nearfield part
 * - if \a delete_farfield is non-zero, farfield matrices (e.g., low-rank) will 
 *   be removed, otherwise set to rank 0 */
HLIB_FNDECL void
HLIBPF(matrix_restrict_nearfield)  ( HLIBPF(matrix_t)              A,
                                     const int                     delete_farfield,
                                     int *                         info );

/** returns nearfield part of given matrix in CRS format
    - CRS is stored in \a nnz, \a rowptr, \a colind and \a coeffs
 */
HLIB_FNDECL void
HLIBPF(matrix_nearfield_to_crs)    ( const HLIBPF(matrix_t)        A,
                                     size_t *                      nnz,    /*!< number of non-zero coeff.    */
                                     int **                        rowptr, /*!< array of row pointers        */
                                     int **                        colind, /*!< array of column indices      */
                                     HLIBPF(real_t) **             coeffs, /*!< array of matrix coefficients */
                                     int *                         info );

/** returns nearfield part of given matrix in CRS format
    - complex valued version */
HLIB_FNDECL void
HLIBPF(matrix_nearfield_to_ccrs)   ( const HLIBPF(matrix_t)        A,
                                     size_t *                      nnz,    /*!< number of non-zero coeff.    */
                                     int **                        rowptr, /*!< array of row pointers        */
                                     int **                        colind, /*!< array of column indices      */
                                     HLIBPF(complex_t) **          coeffs, /*!< array of matrix coefficients */
                                     int *                         info );

/** return size of matrix in bytes */
HLIB_FNDECL size_t
HLIBPF(matrix_bytesize)            ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/** free resources coupled with matrix A */
HLIB_FNDECL void
HLIBPF(matrix_free)                ( HLIBPF(matrix_t)              A,
                                     int *                         info );

/** return 1 if matrix \a A is symmetric, e.g. A = A^T */
HLIB_FNDECL int
HLIBPF(matrix_is_sym)              ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/** return 1 if matrix \a A is hermitian, e.g. A = A^H */
HLIB_FNDECL int
HLIBPF(matrix_is_herm)             ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/** return 1 if matrix \a A is complex valued and 0 otherwise */
HLIB_FNDECL int
HLIBPF(matrix_is_complex)          ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/** switch matrix to real value (if possible); if \a force is 1, the
 *  conversion is performed independent of the current state
 *  (applies to block matrices) */
HLIB_FNDECL void
HLIBPF(matrix_to_real)             ( const HLIBPF(matrix_t)        A,
                                     const int                     force,
                                     int *                         info );

/** switch matrix to complex value; if \a force is 1, the
 *  conversion is performed independent of the current state
 *  (applies to block matrices) */
HLIB_FNDECL void
HLIBPF(matrix_to_complex)          ( const HLIBPF(matrix_t)        A,
                                     const int                     force,
                                     int *                         info );

/** get single entry at position (\a i,\a j) in matrix \a A */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(matrix_entry_get)           ( const HLIBPF(matrix_t)        A,
                                     const size_t                  i,
                                     const size_t                  j,
                                     int *                         info );

/** get single entry at position (\a i,\a j) in matrix \a A */
HLIB_FNDECL HLIBPF(complex_t)
HLIBPF(matrix_centry_get)          ( const HLIBPF(matrix_t)        A,
                                     const size_t                  i,
                                     const size_t                  j,
                                     int *                         info );

/************************************************
 **
 ** matrix dependent vector building
 **
 ************************************************/

/** create vectors matching the row indexsets of matrix \a A */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(matrix_row_vector)          ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/** create vectors matching the column indexsets of matrix \a A */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(matrix_col_vector)          ( const HLIBPF(matrix_t)        A,
                                     int *                         info );

/************************************************
 **
 ** matrix import/building
 **
 ************************************************/

/**
 * import given real sparse matrix in CRS format into
 * internal sparse matrix format
 * - the data is copied into internal representation
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_import_crs)          ( const size_t                  rows,   /*!< number of rows and           */
                                     const size_t                  cols,   /*!< columns of the matrix        */
                                     const size_t                  nnz,    /*!< number of non-zero coeff.    */
                                     const int *                   rowptr, /*!< array of row pointers        */
                                     const int *                   colind, /*!< array of column indices      */
                                     const HLIBPF(real_t) *        coeffs, /*!< array of matrix coefficients */
                                     const int                     sym,    /*!< if non-zero, matrix is sym.  */
                                     int *                         info );

/**
 * import given complex sparse matrix in CRS format into
 * internal sparse matrix format
 * - the data is copied into internal representation
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_import_ccrs)         ( const size_t                  rows,   /*!< number of rows and           */
                                     const size_t                  cols,   /*!< columns of the matrix        */
                                     const size_t                  nnz,    /*!< number of non-zero coeff.    */
                                     const int *                   rowptr, /*!< array of row pointers        */
                                     const int *                   colind, /*!< array of column indices      */
                                     const HLIBPF(complex_t) *     coeffs, /*!< array of matrix coeff.       */
                                     const int                     sym,    /*!< if non-zero, matrix is sym.  */
                                     int *                         info );

/**
 * import given real sparse matrix in CCS format into
 * internal sparse matrix format
 * - the data is copied into internal representation
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_import_ccs)          ( const size_t                  rows,   /*!< number of rows and           */
                                     const size_t                  cols,   /*!< columns of the matrix        */
                                     const size_t                  nnz,    /*!< number of non-zero coeff.    */
                                     const int *                   colptr, /*!< array of column pointers     */
                                     const int *                   rowind, /*!< array of row indices         */
                                     const HLIBPF(real_t) *        coeffs, /*!< array of matrix coefficients */
                                     const int                     sym,    /*!< if non-zero, matrix is sym.  */
                                     int *                         info );

/**
 * import given complex sparse matrix in CCS format into
 * internal sparse matrix format
 * - the data is copied into internal representation
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_import_cccs)         ( const size_t                  rows,   /*!< number of rows and           */
                                     const size_t                  cols,   /*!< columns of the matrix        */
                                     const size_t                  nnz,    /*!< number of non-zero coeff.    */
                                     const int *                   colptr, /*!< array of column pointers     */
                                     const int *                   rowind, /*!< array of row indices         */
                                     const HLIBPF(complex_t) *     coeffs, /*!< array of matrix coeff.       */
                                     const int                     sym,    /*!< if non-zero, matrix is sym.  */
                                     int *                         info );

/**
 * import given real dense matrix in column major format
 * into internal dense matrix format
 * - the data is copied into internal representation
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_import_dense)        ( const size_t                  rows,   /*!< number of rows and            */
                                     const size_t                  cols,   /*!< columns of the matrix         */
                                     const HLIBPF(real_t) *        D,      /*!< accuracy of the approximation */
                                     const int                     sym,    /*!< if non-zero, matrix is sym.   */
                                     int *                         info );

/**
 * import given complex dense matrix in column major format
 * into internal dense matrix format
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_import_cdense)       ( const size_t                  rows,   /*!< number of rows and            */
                                     const size_t                  cols,   /*!< columns of the matrix         */
                                     const HLIBPF(complex_t) *     D,      /*!< accuracy of the approximation */
                                     const int                     sym,    /*!< if non-zero, matrix is sym.   */
                                     int *                         info );

/**
 * build H-matrix over block cluster tree \a bct with contents
 * defined by given sparse matrix \a S
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_build_sparse)        ( const HLIBPF(blockclustertree_t)  bct,    /*!< block ct. defining indexsets  */
                                     const HLIBPF(matrix_t)            S,      /*!< sparse mat. to be converted   */
                                     const HLIBPF(acc_t)               acc,    /*!< accuracy of the approximation */
                                     int *                             info );

/**
 * build H-matrix over block cluster tree \a bct with contents
 * defined by given dense matrix \a D with low-rank approximation
 * method \a lrapx and block-wise accuracy \a acc
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_build_dense)         ( const HLIBPF(blockclustertree_t)  bct,    /*!< block ct. defining indexsets  */
                                     const HLIBPF(matrix_t)            D,      /*!< dense mat. to be converted    */
                                     const HLIBPF(lrapx_t)             lrapx,  /*!< low-rank approximation method */
                                     const HLIBPF(acc_t)               acc,    /*!< accuracy of the approximation */
                                     int *                             info );

/**
 * build H-matrix over block cluster tree \a bct of precision \a acc out
 * of dense matrix given by a matrix coefficient function \a f by using a
 * low-rank approximation for admissible blocks
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_build_coeff)         ( const HLIBPF(blockclustertree_t)  bct,    /*!< block ct. defining indexsets  */
                                     const HLIBPF(coeff_t)             f,      /*!< coeff. function defining mat. */
                                     void *                            arg,    /*!< add. argument to coeff. fn.   */
                                     const HLIBPF(lrapx_t)             lrapx,  /*!< low-rank approximation method */
                                     const HLIBPF(acc_t)               acc,    /*!< accuracy of the approximation */
                                     const int                         sym,    /*!< if non-zero, matrix is sym.   */
                                     int *                             info );

/**
 * build H-matrix over block cluster tree \a bct of precision \a acc out
 * of dense matrix given by a matrix coefficient function \a f by using a
 * low-rank approximation for admissible blocks
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_build_ccoeff)        ( const HLIBPF(blockclustertree_t)  bct,    /*!< block ct. defining indexsets  */
                                     const HLIBPF(ccoeff_t)            f,      /*!< coeff. function defining mat. */
                                     void *                            arg,    /*!< add. argument to coeff. fn.   */
                                     const HLIBPF(lrapx_t)             lrapx,  /*!< low-rank approximation method */
                                     const HLIBPF(acc_t)               acc,    /*!< accuracy of the approximation */
                                     const int                         sym,    /*!< if non-zero, matrix is sym.   */
                                     int *                             info );

/**
 * build H-matrix over block cluster tree \a bct representing identity
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_build_identity)      ( const HLIBPF(blockclustertree_t)  bct,   
                                     int *                             info );

/**
 * assemble \a brows × \a bcols block matrix out of sub blocks \a submat
 * - sub matrices have to be stored column wise in \a submat, e.g.,
 *       \f{eqnarray*}{ A_{00} & A_{01} \\
 *                      A_{10} & A_{11} \f}
 *   is stored { A_00, A_10, A_01, A_11 },
 * - if \a copy is non-zero, sub matrices are copied, otherwise used
 *   directly in new block matrix,
 * - index sets of sub matrices are adjusted to be consistent in new
 *   block matrix,
 * - if sub matrices are H-matrices, permutations from internal to
 *   external ordering (and vice versa) are also adjusted.
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_assemble_block)      ( const size_t                      brows,  /*!< number of block rows         */
                                     const size_t                      bcols,  /*!< number of block columns      */
                                     HLIBPF(matrix_t) *                submat, /*!< array of sub matrices        */
                                     const int                         copy,   /*!< copy matrices if true        */
                                     int *                             info );

/************************************************
 **
 ** manage index ordering
 **
 ************************************************/

/**
 * return mapping from internal (in H-matrix) to external (user) ordering
 * for rows and columns.
 * - returned object is reference to internal permutation: do NOT free
 * - if matrix has no permutation, the return value is NULL
 */
HLIB_FNDECL HLIBPF(permutation_t)
HLIBPF(matrix_row_perm_i2e)        ( const HLIBPF(matrix_t)  A,
                                     int *                   info );

HLIB_FNDECL HLIBPF(permutation_t)
HLIBPF(matrix_col_perm_i2e)        ( const HLIBPF(matrix_t)  A,
                                     int *                   info );

/**
 * return mapping from external (user) to internal (in H-matrix) ordering
 * for rows and columns.
 * - returned object is reference to internal permutation: do NOT free
 * - if matrix has no permutation, the return value is NULL
 */
HLIB_FNDECL HLIBPF(permutation_t)
HLIBPF(matrix_row_perm_e2i)        ( const HLIBPF(matrix_t)  A,
                                     int *                   info );

HLIB_FNDECL HLIBPF(permutation_t)
HLIBPF(matrix_col_perm_e2i)        ( const HLIBPF(matrix_t)  A,
                                     int *                   info );

/** reorder matrix \a A entries with given row and column permutations
 *  - only supported for dense matrices
 *  - NULL permutation will be handled as identity
 */
HLIB_FNDECL void
HLIBPF(matrix_permute)             ( HLIBPF(matrix_t)             A,
                                     const HLIBPF(permutation_t)  row_perm,
                                     const HLIBPF(permutation_t)  col_perm,
                                     int *                        info );

/**
 * represents permuted matrix \f$P \cdot A \cdot R\f$
 * - all objects (\a P, \a A and \a R) are only referenced, not copied
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(perm_matrix)                ( const HLIBPF(permutation_t)  P,
                                     const HLIBPF(matrix_t)       A,
                                     const HLIBPF(permutation_t)  R,
                                     int *                        info );

/**
 * represents product of matrices \f$A_1 \cdot A_2\f$
 * - all matrix objects are only referenced, not copied
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_product2)            ( const HLIBPF(matop_t)           matop1,
                                     const HLIBPF(linearoperator_t)  A1,
                                     const HLIBPF(matop_t)           matop2,
                                     const HLIBPF(linearoperator_t)  A2,
                                     int *                           info );

/**
 * represents product of matrices \f$A_1 \cdot A_2 \cdot A3\f$
 * - all matrix objects are only referenced, not copied
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_product3)            ( const HLIBPF(matop_t)           matop1,
                                     const HLIBPF(linearoperator_t)  A1,
                                     const HLIBPF(matop_t)           matop2,
                                     const HLIBPF(linearoperator_t)  A2,
                                     const HLIBPF(matop_t)           matop3,
                                     const HLIBPF(linearoperator_t)  A3,
                                     int *                           info );

/**
 * represents product of matrices \f$A_1 + A_2\f$
 * - all matrix objects are only referenced, not copied
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_sum2)                ( const HLIBPF(real_t)            alpha1,
                                     const HLIBPF(matop_t)           matop1,
                                     const HLIBPF(linearoperator_t)  A1,
                                     const HLIBPF(real_t)            alpha2,
                                     const HLIBPF(matop_t)           matop2,
                                     const HLIBPF(linearoperator_t)  A2,
                                     int *                           info );

/**
 * represents product of matrices \f$A_1 + A_2 + A3\f$
 * - all matrix objects are only referenced, not copied
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_sum3)                ( const HLIBPF(real_t)            alpha1,
                                     const HLIBPF(matop_t)           matop1,
                                     const HLIBPF(linearoperator_t)  A1,
                                     const HLIBPF(real_t)            alpha2,
                                     const HLIBPF(matop_t)           matop2,
                                     const HLIBPF(linearoperator_t)  A2,
                                     const HLIBPF(real_t)            alpha3,
                                     const HLIBPF(matop_t)           matop3,
                                     const HLIBPF(linearoperator_t)  A3,
                                     int *                           info );

/************************************************
 **
 ** matrix conversion and approximation 
 **
 ************************************************/

/** convert matrix to dense format */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_to_dense)            ( const HLIBPF(matrix_t)          A,
                                     int *                           info );

/** convert matrix to lowrank format using best approximation (SVD)
 *  with truncation defined by acc */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_to_rank)             ( const HLIBPF(matrix_t)          A,
                                     const HLIBPF(acc_t)             acc,
                                     int *                           info );

/** convert matrix to lowrank format using approximation by ACA
 *  with accuracy defined by acc
 *  - currently only dense matrices supported
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_approx_rank_aca)     ( const HLIBPF(matrix_t)          A,
                                     const HLIBPF(acc_t)             acc,
                                     int *                           info );

/************************************************
 **
 ** input/output
 **
 ************************************************/

/**
 * print matrix \a A in postscript format to file \a filename
 * options (may be ORed) :
 * - HLIB_MATIO_SVD     : print singular value decomposition in each block
 * - HLIB_MATIO_ENTRY   : print each entry of matrix           
 * - HLIB_MATIO_PATTERN : print sparsity pattern (non-zero entries)
 */
HLIB_FNDECL void
HLIBPF(matrix_print_ps)            ( const HLIBPF(matrix_t)        A,
                                     const char *                  filename,
                                     const int                     options,  /*!< options defining output    */
                                     int *                         info );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Algebra
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 **
 ** matrix vector multiplication
 **
 ************************************************/

/**
 * compute y := alpha op(A) x + beta y, where
 * op(A) is either A or A^T
 */
HLIB_FNDECL void
HLIBPF(matrix_mulvec)              ( const HLIBPF(real_t)     alpha,
                                     const HLIBPF(matrix_t)   A,
                                     const HLIBPF(vector_t)   x,      /*!< source vector of dimension columns(op(A)) */
                                     const HLIBPF(real_t)     beta,
                                     HLIBPF(vector_t)         y,      /*!< dest. vector of dimension rows(op(A))    */
                                     const HLIBPF(matop_t)    matop,  /*!< matrix operation, e.g. transpose   */
                                     int *                    info );

/**
 * compute y := alpha op(A) x + beta y, where
 * op(A) is either A, A^T or A^H
 */
HLIB_FNDECL void
HLIBPF(matrix_cmulvec)             ( const HLIBPF(complex_t)  alpha,
                                     const HLIBPF(matrix_t)   A,
                                     const HLIBPF(vector_t)   x,      /*!< source vector of dimension columns(op(A)) */
                                     const HLIBPF(complex_t)  beta,
                                     HLIBPF(vector_t)         y,      /*!< dest. vector of dimension rows(op(A))    */
                                     const HLIBPF(matop_t)    matop,  /*!< matrix operation, e.g. transpose   */
                                     int *                    info );

/************************************************
 **
 ** matrix operations
 **
 ************************************************/

/**
 * transpose matrix A
 */
HLIB_FNDECL void
HLIBPF(matrix_transpose)           ( const HLIBPF(matrix_t)   A,
                                     int *                    info );

/**
 * conjugate matrix A
 */
HLIB_FNDECL void
HLIBPF(matrix_conjugate)           ( const HLIBPF(matrix_t)   A,
                                     int *                    info );

/**
 * scale matrix \a A by factor \a f
 */
HLIB_FNDECL void
HLIBPF(matrix_scale)               ( const HLIBPF(real_t)     f,
                                     const HLIBPF(matrix_t)   A,
                                     int *                    info );

/**
 * scale matrix \a A by factor \a f
 */
HLIB_FNDECL void
HLIBPF(matrix_cscale)              ( const HLIBPF(complex_t)  f,
                                     const HLIBPF(matrix_t)   A,
                                     int *                    info );

/**
 * compute B := alpha A + beta B with block-wise precision \a acc
 */
HLIB_FNDECL void
HLIBPF(matrix_add)                 ( const HLIBPF(real_t)     alpha,
                                     const HLIBPF(matrix_t)   A,
                                     const HLIBPF(real_t)     beta,
                                     HLIBPF(matrix_t)         B,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * compute B := alpha A + beta B with block-wise precision \a acc
 */
HLIB_FNDECL void
HLIBPF(matrix_cadd)                ( const HLIBPF(complex_t)  alpha,
                                     const HLIBPF(matrix_t)   A,
                                     const HLIBPF(complex_t)  beta,
                                     HLIBPF(matrix_t)         B,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * compute C := alpha op(A) op(B) + beta C with block-wise precision \a acc
 */
HLIB_FNDECL void
HLIBPF(matrix_mul)                 ( const HLIBPF(real_t)     alpha,
                                     const HLIBPF(matop_t)    matop_A,
                                     const HLIBPF(matrix_t)   A,
                                     const HLIBPF(matop_t)    matop_B,
                                     const HLIBPF(matrix_t)   B,
                                     const HLIBPF(real_t)     beta,
                                     HLIBPF(matrix_t)         C,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * compute C := alpha op(A) op(B) + beta C with block-wise precision \a acc
 */
HLIB_FNDECL void
HLIBPF(matrix_cmul)                ( const HLIBPF(complex_t)  alpha,
                                     const HLIBPF(matop_t)    matop_A,
                                     const HLIBPF(matrix_t)   A,
                                     const HLIBPF(matop_t)    matop_B,
                                     const HLIBPF(matrix_t)   B,
                                     const HLIBPF(complex_t)  beta,
                                     HLIBPF(matrix_t)         C,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * compute inverse of A with block-wise precision \a acc;
 * A will be overwritten with A^-1 upon exit
 * - the return value is either A, or a new representation
 *   for the inverse of A (A is still needed then)
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matrix_inv)                 ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * compute diagonal of inverse of \a A with local accuracy \a acc
 * - \a A will be overwritten during computation
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(matrix_inv_diag)            ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * factorise matrix \a A up to block-wise precision \a acc
 * - depending on form of \a A, e.g. if unsymmetric, symmetric
 *   or hermitian, an appropriate factorisation method is chosen
 * - \a A will be overwritten by the factors
 * - the return value is a matrix which can be used to
 *   evaluate the factorisation, e.g. for matrix-vector mult.
 *   (this is not possible with A after factorisation
 *    as it holds just the data, albeit, it is still needed)
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_factorise)           ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * same as hlib_matrix_factorise but return matrix object for
 * evaluation of inverse of \a A.
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_factorise_inv)       ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * LU factorise matrix A up to block-wise precision \a acc
 * - A will be overwritten by L*U, which only allows
 *   matrix vector mult.
 * - the return value is a matrix which can be used to
 *   evaluate LU, e.g. for matrix-vector multiplication
 *   (this is not possible with A after factorisation
 *    as it holds just the data, albeit, it is still needed)
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_lu)                  ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * LDL factorisation of matrix A up to block-wise precision \a acc
 * - A will be overwritten by L*D*L^T, which only allows
 *   matrix vector mult.
 * - in case of hermitian matrices, LDL^H is computed
 * - for return value see "hlib_matrix_lu"
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_ldl)                 ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * Cholesky factorisation of matrix A up to block-wise precision \a acc
 * - A will be overwritten by L*L^T, which only allows
 *   matrix vector mult.
 * - in case of hermitian matrices, L*L^H is computed
 * - for return value see "hlib_matrix_lu"
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_ll)                  ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_chol)                ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * Factorise matrix A into \f$WAZ=I\f$, e.g., factorization of inverse of A
 * - A will be overwritten by \f$W^{-1} Z^{-1}\f$.
 * - the return value is a linear operator which can be used to
 *   evaluate \f$A^{-1}\f$, e.g. for matrix-vector multiplication
 *   (this is not possible with A after factorisation
 *    as it holds just the data, albeit, it is still needed)
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_waz)                 ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * LU factorise matrix A up to block-wise precision \a acc
 * - A will be overwritten by L*U, which only allows
 *   matrix vector mult.
 * - the return value is a matrix which can be used to
 *   evaluate (LU)^-1, e.g. for matrix-vector multiplication
 *   (this is not possible with A after factorisation
 *    as it holds just the data, albeit, it is still needed)
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_lu_inv)              ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * LDL factorisation of matrix A up to block-wise precision \a acc
 * - A will be overwritten by L*D*L^T, which only allows
 *   matrix vector mult.
 * - in case of hermitian matrices, LDL^H is computed
 * - for return value see "hlib_matrix_lu_inv"
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_ldl_inv)             ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * Cholesky factorisation of matrix A up to block-wise precision \a acc
 * - A will be overwritten by L*L^T, which only allows
 *   matrix vector mult.
 * - in case of hermitian matrices, L*L^H is computed
 * - for return value see "hlib_matrix_lu_inv"
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_ll_inv)              ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(matrix_chol_inv)            ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(acc_t)      acc,
                                     int *                    info );

/**
 * solve L B = A with lower triangular L and known A
 * - B will overwrite A
 * - if \a pointwise is true, evaluation is performed pointwise, otherwise
 *   blockwise, e.g. diagonal leaf matrices are considered full and not
 *   triangular
 * - if \a unit_diag is true, the pointwise/blockwise diagonal is considered
 *   to be the identity
 */
HLIB_FNDECL void
HLIBPF(matrix_solve_lower_left)    ( const HLIBPF(matrix_t)     L,
                                     const HLIBPF(matop_t)      op_L,
                                     HLIBPF(matrix_t)           A,
                                     const HLIBPF(acc_t)        acc,
                                     const int                  pointwise,
                                     const int                  unit_diag,
                                     int *                      info );

/**
 * solve B op(L) = A with lower triangular L and known A
 * - B will overwrite A
 * - if \a pointwise is true, evaluation is performed pointwise, otherwise
 *   blockwise, e.g. diagonal leaf matrices are considered full and not
 *   triangular
 * - if \a unit_diag is true, the pointwise/blockwise diagonal is considered
 *   to be the identity
 */
HLIB_FNDECL void
HLIBPF(matrix_solve_lower_right)   ( HLIBPF(matrix_t)           A,
                                     const HLIBPF(matrix_t)     L,
                                     const HLIBPF(matop_t)      op_L,
                                     const HLIBPF(acc_t)        acc,
                                     const int                  pointwise,
                                     const int                  unit_diag,
                                     int *                      info );

/**
 * solve B U = A with upper triangular U and known A
 * - B will overwrite A
 * - if \a pointwise is true, evaluation is performed pointwise, otherwise
 *   blockwise, e.g. diagonal leaf matrices are considered full and not
 *   triangular
 * - if \a unit_diag is true, the pointwise/blockwise diagonal is considered
 *   to be the identity
 */
HLIB_FNDECL void
HLIBPF(matrix_solve_upper_right)   ( HLIBPF(matrix_t)           A,
                                     const HLIBPF(matrix_t)     U,
                                     const HLIBPF(acc_t)        acc,
                                     const int                  pointwise,
                                     const int                  unit_diag,
                                     int *                      info );

/**
 * solve op(D) B = A with diagonal D and known A
 * - B will overwrite A
 * - if \a pointwise is true, evaluation is performed pointwise, otherwise
 *   blockwise, e.g. diagonal leaf matrices are considered full and not
 *   diagonal
 */
HLIB_FNDECL void
HLIBPF(matrix_solve_diag_left)     ( const HLIBPF(matrix_t)     D,
                                     const HLIBPF(matop_t)      op_D,
                                     HLIBPF(matrix_t)           A,
                                     const HLIBPF(acc_t)        acc,
                                     const int                  pointwise,
                                     int *                      info );

/**
 * solve B op(D) = A with diagonal D and known A
 * - B will overwrite A
 * - if \a pointwise is true, evaluation is performed pointwise, otherwise
 *   blockwise, e.g. diagonal leaf matrices are considered full and not
 *   diagonal
 */
HLIB_FNDECL void
HLIBPF(matrix_solve_diag_right)    ( HLIBPF(matrix_t)           A,
                                     const HLIBPF(matrix_t)     D,
                                     const HLIBPF(matop_t)      op_D,
                                     const HLIBPF(acc_t)        acc,
                                     const int                  pointwise,
                                     int *                      info );

/**
 * compute A := A + lambda I
 */
HLIB_FNDECL void
HLIBPF(matrix_add_identity)        ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(real_t)     lambda,
                                     int *                    info );

/**
 * compute A := A + lambda I
 */
HLIB_FNDECL void
HLIBPF(matrix_cadd_identity)       ( HLIBPF(matrix_t)         A,
                                     const HLIBPF(complex_t)  lambda,
                                     int *                    info );

/************************************************
 **
 ** norm computations
 **
 ************************************************/

/** compute Frobenius norm of matrix */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(matrix_norm_frobenius)      ( const HLIBPF(matrix_t)   A,
                                     int *                    info );

/** compute Frobenius norm of A-B */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(matrix_norm_frobenius_diff) ( const HLIBPF(matrix_t)   A,
                                     const HLIBPF(matrix_t)   B,
                                     int *                    info );

/** compute spectral norm of matrix */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(matrix_norm_spectral)       ( const HLIBPF(matrix_t)   A,
                                     int * info );

/** compute spectral norm of inverse matrix */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(matrix_norm_spectral_inv)   ( const HLIBPF(matrix_t)   A,
                                     int *                    info );

/** compute relative spectral norm of A-B, e.g. |A-B|_2/|A|_2 */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(matrix_norm_spectral_diff)  ( const HLIBPF(matrix_t)   A,
                                     const HLIBPF(matrix_t)   B,
                                     int *                    info );

/** compute inversion error |I-BA|_2, where B is supposed to be A^-1 */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(matrix_norm_inv_approx)     ( const HLIBPF(matrix_t)   A,
                                     const HLIBPF(matrix_t)   B,
                                     int *                    info );

/************************************************
 **
 ** norm computations with linear operators
 **
 ************************************************/

/** compute spectral norm of matrix */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(linearoperator_norm_spectral)       ( const HLIBPF(linearoperator_t)   A,
                                             int *                            info );

/** compute spectral norm of inverse matrix */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(linearoperator_norm_spectral_inv)   ( const HLIBPF(linearoperator_t)   A,
                                             int *                            info );

/** compute relative spectral norm of A-B, e.g. |A-B|_2/|A|_2 */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(linearoperator_norm_spectral_diff)  ( const HLIBPF(linearoperator_t)   A,
                                             const HLIBPF(linearoperator_t)   B,
                                             int *                            info );

/** compute inversion error |I-BA|_2, where B is supposed to be A^-1 */
HLIB_FNDECL HLIBPF(real_t)
HLIBPF(linearoperator_norm_inv_approx)     ( const HLIBPF(linearoperator_t)   A,
                                             const HLIBPF(linearoperator_t)   B,
                                             int *                            info );

/** @} */


/***********************************************************//**
 ***************************************************************
 *
 * @{ \name Solver
 *
 ***************************************************************
 ***************************************************************/

/************************************************
 *
 * solver types
 *
 ************************************************/

/** create automatic solver */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_auto)                ( int *                  info );

/** create linear iteration */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_linear_iteration)    ( int *                  info );

/** create Richardson solver */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_richardson)          ( int *                  info );

/** create CG solver */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_cg)                  ( int *                  info );

/** create CGS solver */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_cgs)                 ( int *                  info );

/** create BiCG-Stab solver */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_bicgstab)            ( int *                  info );

/** create TFQMR solver */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_tfqmr)               ( int *                  info );

/** create MINRES solver */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_minres)              ( int *                  info );

/** create GMRES solver with restart after \a restart steps */
HLIB_FNDECL HLIBPF(solver_t)
HLIBPF(solver_gmres)               ( const int              restart,
                                     int *                  info );

/************************************************
 *
 * solving
 *
 ************************************************/

/**
 * solve A x = b using \a solver with optional preconditioning,
 * e.g. W A x = W b if W != NULL; information about the solving
 * process is stored in \a solve_info
 */
HLIB_FNDECL void
HLIBPF(solver_solve)               ( const HLIBPF(solver_t)           solver,
                                     const HLIBPF(linearoperator_t)   A,
                                     HLIBPF(vector_t)                 x,
                                     const HLIBPF(vector_t)           b,
                                     const HLIBPF(linearoperator_t)   W,
                                     HLIBPF(solve_info_t) *           solve_info,
                                     int *                            info );

/**
 * set stopping criterion for \a solver
 * - \a maxit   : maximal number of iterations
 * - \a abs_red : iterate until ||residual||_2 reaches \a abs_red
 * - \a rel_red : iterate until ||residual||_2 reaches \a rel_red * ||start-res.||_2
 * - in case of a preconditioner r = W(Ax - b), e.g. the preconditioned res.
 */
HLIB_FNDECL void
HLIBPF(solver_stopcrit)            ( HLIBPF(solver_t)                 solver,
                                     const int                        maxit,
                                     const HLIBPF(real_t)             abs_red,
                                     const HLIBPF(real_t)             rel_red,
                                     int *                            info );

/**
 * turn initialisation of start value during iteration on (\a flag != 0)
 * or off (\a flag != 0)
 */
HLIB_FNDECL void
HLIBPF(solver_initialise_start_value) ( HLIBPF(solver_t)              solver,
                                        const int                     flag,
                                        int *                         info );

/**
 * turn computation of exact residual during iteration on (\a flag != 0)
 * or off (\a flag != 0)
 */
HLIB_FNDECL void
HLIBPF(solver_use_exact_residual)  ( HLIBPF(solver_t)                 solver,
                                     const int                        flag,
                                     int *                            info );

/************************************************
 *
 * management
 *
 ************************************************/

/** free \a solver */
HLIB_FNDECL void
HLIBPF(solver_free)                ( HLIBPF(solver_t)         solver,
                                     int *                    info );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Input/Output
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 *
 * HLib format
 *
 ************************************************/

/**
 * read matrix from file \a filename
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(hformat_load_matrix)         ( const char *             filename,
                                      int *                    info );

/**
 * read matrix from file \a filename
 */
HLIB_FNDECL HLIBPF(linearoperator_t)
HLIBPF(hformat_load_linearoperator) ( const char *             filename,
                                      int *                    info );

/**
 * read vector from file \a filename
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(hformat_load_vector)         ( const char *             filename,
                                      int *                    info );

/**
 * save linear operator \a A to file \a filename
 */
HLIB_FNDECL void
HLIBPF(hformat_save_linearoperator) ( const HLIBPF(linearoperator_t)  A,
                                      const char *                    filename,
                                      int *                           info );

/**
 * save matrix \a A to file \a filename
 */
HLIB_FNDECL void
HLIBPF(hformat_save_matrix)         ( const HLIBPF(matrix_t)   A,
                                      const char *             filename,
                                      int *                    info );

/**
 * save vector \a v to file \a filename
 */
HLIB_FNDECL void
HLIBPF(hformat_save_vector)         ( const HLIBPF(vector_t)   v,
                                      const char *             filename,
                                      int *                    info );

/**
 * read coordinates from file \a filename
 */
HLIB_FNDECL HLIBPF(coord_t)
HLIBPF(hformat_load_coord)          ( const char *             filename,
                                      int *                    info );

/**
 * save coordinates to file \a filename
 */
HLIB_FNDECL void
HLIBPF(hformat_save_coord)          ( const HLIBPF(coord_t)    coord,
                                      const char *             filename,
                                      int *                    info );


/************************************************
 *
 * SAMG format
 *
 ************************************************/

/**
 * read matrix \a A from file "\a filename" in SAMG format
 * - expecting format definition in basename(\a filename)>.frm
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(samg_load_matrix)           ( const char *             filename,
                                     int        *             info );

/**
 * save matrix \a A to file "\a filename" in SAMG format
 * - also writing format definition to basename(\a filename)>.frm
 */
HLIB_FNDECL void
HLIBPF(samg_save_matrix)           ( const HLIBPF(matrix_t)   A,
                                     const char *             filename,
                                     int *                    info );

/**
 * read vector from file \a filename in SAMG format
 * - expecting format definition in basename(\a filename)>.frm
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(samg_load_vector)           ( const char *             filename,
                                     int        *             info );

/**
 * save vector to file \a filename in SAMG format
 */
HLIB_FNDECL void
HLIBPF(samg_save_vector)            ( const HLIBPF(vector_t)  x,
                                      const char *            filename,
                                      int *                   info );

/**
 * write coordinates \a coord to file \a filename
 */
HLIB_FNDECL void
HLIBPF(samg_save_coord)            ( const HLIBPF(coord_t)    coord,
                                     const char *             filename,
                                     int *                    info );

/**
 * read coordinates from file \a filename
 */
HLIB_FNDECL HLIBPF(coord_t)
HLIBPF(samg_load_coord)            ( const char *             filename,
                                     int *                    info );

/************************************************
 *
 * Matlab format V7 (with libz)
 *
 ************************************************/

/**
 * read matrix named \a matname in Matlab V7 format
 * from file \a filename; if \a matname == NULL or \a matname == "",
 * the first matrix in \a filename will be read
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(matlab_load_matrix)         ( const char *             filename,
                                     const char *             matname,
                                     int *                    info );

/**
 * save matrix \a M named \a matname in Matlab V7 format to
 * file \a filename
 */
HLIB_FNDECL void
HLIBPF(matlab_save_matrix)         ( const HLIBPF(matrix_t)   M,
                                     const char *             filename,
                                     const char *             matname,
                                     int *                    info );

/**
 * read vector named \a vname from file \a filename in Matlab V7 format
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(matlab_load_vector)         ( const char *             filename,
                                     const char *             vecname,
                                     int *                    info );

/**
 * save vector named \a vname to file \a filename in Matlab V7 format
 */
HLIB_FNDECL void
HLIBPF(matlab_save_vector)         ( const HLIBPF(vector_t)   v,
                                     const char *             filename,
                                     const char *             vecname,
                                     int *                    info );

/************************************************
 *
 * PLTMG format
 *
 ************************************************/

/**
 * read matrix from file \a filename
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(pltmg_load_matrix)          ( const char *             filename,
                                     int *                    info );

/************************************************
 *
 * Harwell Boeing format
 *
 ************************************************/

/**
 * read matrix from file \a filename
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(hb_load_matrix)             ( const char *             filename,
                                     int *                    info );

/**
 * save matrix \a M to file \a filename
 */
HLIB_FNDECL void
HLIBPF(hb_save_matrix)             ( const HLIBPF(matrix_t)   M,
                                     const char *             filename,
                                     int *                    info );

/************************************************
 *
 * Matrix Market format
 *
 ************************************************/

/**
 * read matrix from file \a filename
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(mm_load_matrix)            ( const char *              filename,
                                    int *                     info );


/************************************************
 *
 * general form (autodetect format)
 *
 ************************************************/

/**
 * read matrix from file \a filename
 */
HLIB_FNDECL HLIBPF(matrix_t)
HLIBPF(load_matrix)                ( const char *             filename,
                                     int *                    info );

/**
 * read vector from file \a filename
 */
HLIB_FNDECL HLIBPF(vector_t)
HLIBPF(load_vector)                ( const char *             filename,
                                     int *                    info );

/**
 * read coordinates from file \a filename
 */
HLIB_FNDECL HLIBPF(coord_t)
HLIBPF(load_coord)                 ( const char *             filename,
                                     int *                    info );

/** @} */


/***********************************************************//**
 ***************************************************************
 *
 * @{ \name Accuracy Management
 *
 ***************************************************************
 ***************************************************************/

/**
 * return accuracy object with fixed accuracy \a eps
 */
HLIB_FNDECL HLIBPF(acc_t)
HLIBPF(acc_fixed_eps)              ( const HLIBPF(real_t)   eps );

/**
 * return accuracy object with fixed rank \a k
 */
HLIB_FNDECL HLIBPF(acc_t)
HLIBPF(acc_fixed_rank)             ( const unsigned int     k );

/**
 * return accuracy object with blockwise accuracy defined
 * by array \a blockacc of dimension blockrows × blockcolumns
 * corresponding to first partitioning level; the array
 * must be stored column wise
 */
HLIB_FNDECL HLIBPF(acc_t)
HLIBPF(acc_blocked)                ( const HLIBPF(acc_t) *  blockacc );

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Misc.
 **
 ***************************************************************
 ***************************************************************/

/************************************************
 **
 ** timing
 **
 ************************************************/

/** return current walltime in seconds */
HLIB_FNDECL double
HLIBPF(walltime)                   ();

/** return execution time of program in seconds */
HLIB_FNDECL double
HLIBPF(cputime)                    ();

/************************************************
 **
 ** change internal parameters
 **
 ************************************************/

/** set maximal leaf size in cluster trees to \a n */
HLIB_FNDECL void
HLIBPF(set_n_min)                  ( const unsigned int       n );

/** define verbosity of HLib */
HLIB_FNDECL void
HLIBPF(set_verbosity)              ( const unsigned int       verb );

/** define minimal boundary for singular values, e.g. not smaller */
HLIB_FNDECL void
HLIBPF(set_abs_eps)                ( const HLIBPF(real_t)     eps );

/** enable (1) or disable (0) coarsening during H-matrix building
 *  and arithmetic (default: build=1, arith=0) */
HLIB_FNDECL void
HLIBPF(set_coarsening)             ( const int                build,
                                     const int                arith );

/** enable (1) or disable (0) recompression during H-matrix building
 *  (default: on) */
HLIB_FNDECL void
HLIBPF(set_recompress)             ( const int                recompress );

/** enable (1) or disable (0) diagonal scaling during LU, etc.
 * (default: off) */
HLIB_FNDECL void
HLIBPF(set_diag_scale)             ( const int                scale );

/** set maximal number of threads to use in H arithmetic */
HLIB_FNDECL void
HLIBPF(set_nthreads)               ( const unsigned int       p );

/** set callback function for progress information     
 *  if fn == NULL, HLIBpro uses default progress output */
HLIB_FNDECL void
HLIBPF(set_progress_cb)            ( HLIBPF(progressfn_t)     fn,
                                     void *                   arg );

/** general function to change configuration variables:
 *  set \a option to \a value */
HLIB_FNDECL void
HLIBPF(set_config)                 ( const char *             option,
                                     const char *             value,
                                     int *                    info );

/** return current value of \a option, stored as string in \a value */
HLIB_FNDECL void
HLIBPF(get_config)                 ( const char *             option,
                                     char *                   value,
                                     const size_t             len,
                                     int *                    info );

/** print all internal parameters of HLIBpro with their current value */
HLIB_FNDECL void
HLIBPF(print_parameters)           ();

/************************************************
 **
 ** error/warning handling
 **
 ************************************************/

/** set call back function in case of an error */
HLIB_FNDECL void
HLIBPF(set_error_fn)               ( const HLIBPF(errorfn_t)  errorfn );

/** set call back function in case of a warning */
HLIB_FNDECL void
HLIBPF(set_warning_fn)             ( const HLIBPF(errorfn_t)  warnfn );

/************************************************
 **
 ** version information
 **
 ************************************************/

/** return major version number of \HLIBpro */
HLIB_FNDECL unsigned int
HLIBPF(major_version)              ();

/** return minor version number of \HLIBpro */
HLIB_FNDECL unsigned int
HLIBPF(minor_version)              ();

/** @} */


/***********************************************************//**
 ***************************************************************
 **
 ** @{ \name Functions for hlib_complex_t
 **
 ***************************************************************
 ***************************************************************/

#if !defined(__cplusplus) && defined(__STDC_VERSION__) &&  __STDC_VERSION__ >= 199901L
#include <complex.h>
#endif

/** @cond */

#if defined(__GNUC__) || defined(__ICC) || defined(__ECC) || defined(__clang__)
#define HLIB_INLINE  inline
#else
#define HLIB_INLINE  static
#endif

/* disable annoying warnings */
#if defined(__cplusplus) && (defined(__GNUC__) || defined(__clang__))
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#endif

/** @endcond */

/**
 * construct hlib_complex_t variable
 */
HLIB_INLINE HLIBPF(complex_t)
HLIBPF(complex)                    ( const HLIBPF(real_t)  re,
                                     const HLIBPF(real_t)  im )
{
    const HLIBPF(complex_t)  z = { ((HLIBPF(real_t)) re),
                                   ((HLIBPF(real_t)) im) };

    return z;
}

/**
 * conversion to/from C99 complex if available
 */
#if !defined(__cplusplus) && defined(__STDC_VERSION__) &&  __STDC_VERSION__ >= 199901L
HLIB_INLINE HLIBPF(complex_t)
HLIBPF(from_c99cmplx)              ( const double _Complex    f )
{
    const HLIBPF(complex_t)  z = { creal( f ), cimag( f ) };

    return z;
}

HLIB_INLINE double _Complex
HLIBPF(to_c99cmplx)                ( const HLIBPF(complex_t)  f )
{
    return f.re + f.im * _Complex_I;
}
#endif

/**
 * standard mathematical operators
 */
HLIB_INLINE HLIBPF(complex_t)
HLIBPF(cmplx_add)                  ( const HLIBPF(complex_t)  a,
                                     const HLIBPF(complex_t)  b )
{
    const HLIBPF(complex_t)  c = { a.re + b.re,
                                a.im + b.im };
    return c;
}

HLIB_INLINE HLIBPF(complex_t)
HLIBPF(cmplx_sub)                  ( const HLIBPF(complex_t)  a,
                                     const HLIBPF(complex_t)  b )
{
    const HLIBPF(complex_t)  c = { a.re - b.re,
                                a.im - b.im };
    return c;
}

HLIB_INLINE HLIBPF(complex_t)
HLIBPF(cmplx_mul)                  ( const HLIBPF(complex_t)  a,
                                     const HLIBPF(complex_t)  b )
{
    const HLIBPF(complex_t)  c = { a.re * b.re - a.im * b.im,
                                a.re * b.im + a.im * b.re };
    return c;
}

HLIB_INLINE HLIBPF(complex_t)
HLIBPF(cmplx_div)                  ( const HLIBPF(complex_t)  a,
                                     const HLIBPF(complex_t)  b )
{
    const HLIBPF(real_t)     n = b.re * b.re + b.im * b.im;
    const HLIBPF(complex_t)  c = { (a.re * b.re + a.im * b.im) / n,
                                (a.im * b.re - a.re * b.im) / n };
    return c;
}

/**
 * return conjugate value
 */
HLIB_INLINE HLIBPF(complex_t)
HLIBPF(cmplx_conj)                 ( const HLIBPF(complex_t)  a )
{
    const HLIBPF(complex_t)  c = { a.re, - a.im };
    return c;
}

/**
 * return absolute value
 */
HLIB_INLINE HLIBPF(real_t)
HLIBPF(cmplx_abs)                  ( const HLIBPF(complex_t)  a )
{
    const HLIBPF(real_t) s = ((HLIBPF(real_t)) ( fabs(a.re) > fabs(a.im) ? fabs(a.re) : fabs(a.im) ));

    if ( s == 0.0 )
        return 0.0;
    else
    {
        const HLIBPF(real_t) x = a.re / s;
        const HLIBPF(real_t) y = a.im / s;

        return s * ((HLIBPF(real_t)) sqrt( x * x + y * y));
    }
}

/**
 * return square root
 */
HLIB_INLINE HLIBPF(complex_t)
HLIBPF(cmplx_sqrt)                 ( const HLIBPF(complex_t)  a )
{
    const HLIBPF(real_t) x = a.re;
    const HLIBPF(real_t) y = a.im;

    if ( x == 0 )
    {
        const HLIBPF(real_t)     t = ((HLIBPF(real_t)) sqrt( fabs( y ) / 2.0 ));
        const HLIBPF(complex_t)  c = { t, ( y < 0 ? -t : t ) };

        return c;
    }
    else
    {
        const HLIBPF(real_t) t = ((HLIBPF(real_t)) sqrt( 2.0 * ( HLIBPF(cmplx_abs)( a ) + fabs(x) ) ));
        const HLIBPF(real_t) u = t / ((HLIBPF(real_t)) 2);

        if ( x > 0 )
        {
            const HLIBPF(complex_t)  c = { u, y / t };
            return c;
        }
        else
        {
            const HLIBPF(complex_t)  c = { ((HLIBPF(real_t)) fabs(y) / t),
                                           ( y < 0 ? -u : u ) };
            return c;
        }
    }
}

/**
 * exponential function e^a
 */
HLIB_INLINE HLIBPF(complex_t)
HLIBPF(cmplx_exp)                  ( const HLIBPF(complex_t)  a )
{
    HLIBPF(complex_t)  res;

    if ( a.re == 0.0 )
    {
        res.re = ((HLIBPF(real_t)) cos(a.im));
        res.im = ((HLIBPF(real_t)) sin(a.im));
    }
    else
    {
        const HLIBPF(real_t) t = ((HLIBPF(real_t)) exp( a.re ));

        res.re = t * ((HLIBPF(real_t)) cos(a.im));
        res.im = t * ((HLIBPF(real_t)) sin(a.im));
    }

    return res;
}

/* disable annoying warnings */
#if defined(__cplusplus) && (defined(__GNUC__) || defined(__clang__))
#pragma clang diagnostic pop
#endif

#ifdef __cplusplus
}/* extern */
#endif

/** @} */

/** @} */    /* finish C-Binding doxygen module */

#endif  /* __HLIB_C_H */
