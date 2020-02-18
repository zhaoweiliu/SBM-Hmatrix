#ifndef __HLIB_BLAS_LAPACK_HH
#define __HLIB_BLAS_LAPACK_HH
//
// Project     : HLib
// File        : lapack.hh
// Description : definition of LAPACK functions in c-format
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hlib-config.h"

// prevents issues with Windows build environment and dot-wrappers
#if USE_MKL == 1
#include <mkl_cblas.h>
#endif

#include "hpro/base/traits.hh"
#include "hpro/base/types.hh"
#include "hpro/base/complex.hh"

#include "hpro/blas/flops.hh"

namespace HLIB
{

namespace BLAS
{

//
// FLOP counting
//

extern double  FLOPS;

#if HLIB_COUNT_FLOPS == 1
#  define ADD_FLOPS( n )  HLIB::BLAS::FLOPS += (n)
#else
#  define ADD_FLOPS( n )
#endif

}// namespace blas

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// BLAS integer type
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

// define 32-bit integers vs. 64-bit integers
#if HAS_ILP64 == 1
using  blas_int_t = long;   // ILP64
#else
using  blas_int_t = int;    // LP64
#endif

/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
//
// declaration of external BLAS functions
//
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////

//
// Win32 DLL export declaration
///
#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
// #  define  CFUNCDECL __declspec(dllimport)
#  define  CFUNCDECL
#else
#  define  CFUNCDECL
#endif

extern "C" {

///////////////////////////////////////////////////////////////////
//
// real-valued functions (single precision)
//
///////////////////////////////////////////////////////////////////

// set (part of) A to alpha/beta
CFUNCDECL
void
slaset_  ( const char *         UPLO,
           const blas_int_t *   M,
           const blas_int_t *   N,
           const float *        ALPHA,
           const float *        BETA,
           float *              A,
           const blas_int_t *   LDA   );

// copy (part of) A to B
CFUNCDECL
void
slacpy_  ( const char *         UPLO,
           const blas_int_t *   M,
           const blas_int_t *   N,
           const float *        A,
           const blas_int_t *   LDA,
           float *              B,
           const blas_int_t *   LDB  );

// compute dotproduct between x and y
CFUNCDECL
float
sdot_ ( const blas_int_t *      n,
        const float *           dx,
        const blas_int_t *      incx,
        const float *           dy,
        const blas_int_t *      incy);

// copy vector x into vector y
CFUNCDECL
void
scopy_ ( const blas_int_t *     n,
         const float *          dx,
         const blas_int_t *     incx,
         float       *          dy,
         const blas_int_t *     incy);

// compute y = y + a * x
CFUNCDECL
void
saxpy_ ( const blas_int_t *     n,
         const float *          da,
         const float *          dx,
         const blas_int_t *     incx,
         float *                dy,
         const blas_int_t *     incy );

// compute sum of absolut values of x
CFUNCDECL
float
sasum_  ( const blas_int_t *    n,
          const float *         dx,
          const blas_int_t *    incx );

// return euclidean norm of x
CFUNCDECL
float
snrm2_  ( const blas_int_t *    n,
          const float *         dx,
          const blas_int_t *    incx );

// scale x by a
CFUNCDECL
void
sscal_ ( const blas_int_t *     n,
         const float *          da,
         float *                dx,
         const blas_int_t *     incx );

// interchange x and y
CFUNCDECL
void
sswap_ ( const blas_int_t *     n,
         float *                dx,
         const blas_int_t *     incx,
         float *                dy,
         const blas_int_t *     incy );

// finds index with element having maximal absolut value
CFUNCDECL
blas_int_t
isamax_ ( const blas_int_t *    n,
          const float   *,
          const blas_int_t * );

// y = alpha A x + beta y
CFUNCDECL
void
sgemv_ ( const char *           trans,
         const blas_int_t *     M,
         const blas_int_t *     N,
         const float *          alpha,
         const float *          A,
         const blas_int_t *     lda,
         const float *          dx,
         const blas_int_t *     incx,
         const float *          beta,
         float *                dy,
         const blas_int_t *     incy );

// compute c = alpha * a * b + beta * c
CFUNCDECL
void
sgemm_ ( const char *           transa,
         const char *           transb,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const blas_int_t *     k,
         const float *          alpha,
         const float *          a,
         const blas_int_t *     lda,
         const float *          b,
         const blas_int_t *     ldb,
         const float *          beta,
         float *                c,
         const blas_int_t *     ldc );

// computes y = op(A)*x for upper/lower triangular A (y overwrites x)
blas_int_t
strmv_ ( const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const float *          A,
         const blas_int_t *     ldA,
         float *                x,
         const blas_int_t *     incx );

// computes B = alpha * op(A) * B or B = alpha * B * op(A)
CFUNCDECL
void
strmm_ ( const char *           side,
         const char *           uplo,
         const char *           transa,
         const char *           diag,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const float *          alpha,
         const float *          A,
         const blas_int_t *     lda,
         float *                B,
         const blas_int_t *     ldb );

// solves A x = b or A^T x = b with triangular A
CFUNCDECL
void
strsv_ ( const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const float *          A,
         const blas_int_t *     lda,
         float *                x,
         const blas_int_t *     incx );

// solves A X = B or A^T X = B with triangular A
CFUNCDECL
void
strsm_ ( const char *           side,
         const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const float *          alpha,
         const float *          A,
         const blas_int_t *     lda,
         float *                B,
         const blas_int_t *     ldb );

// A = alpha * x * y^T + A, A \in \R^{m x n}
CFUNCDECL
void
sger_ ( const blas_int_t *      m,
        const blas_int_t *      n,
        const float *           alpha,
        const float *           x,
        const blas_int_t *      incx,
        const float *           y,
        const blas_int_t *      incy,
        float *                 A,
        const blas_int_t *      lda );

///////////////////////////////////////////////////////////////////
//
// real-valued functions (double precision)
//
///////////////////////////////////////////////////////////////////

// set (part of) A to alpha/beta
CFUNCDECL
void
dlaset_  ( const char *         UPLO,
           const blas_int_t *   M,
           const blas_int_t *   N,
           const double *       ALPHA,
           const double *       BETA,
           double *             A,
           const blas_int_t *   LDA   );

// copy (part of) A to B
CFUNCDECL
void
dlacpy_  ( const char *         UPLO,
           const blas_int_t *   M,
           const blas_int_t *   N,
           const double *       A,
           const blas_int_t *   LDA,
           double *             B,
           const blas_int_t *   LDB  );

// compute dotproduct between x and y
CFUNCDECL
double
ddot_ ( const blas_int_t *      n,
        const double *          dx,
        const blas_int_t *      incx,
        const double *          dy,
        const blas_int_t *      incy);
    
// copy vector x into vector y
CFUNCDECL
void
dcopy_ ( const blas_int_t *     n,
         const double *         dx,
         const blas_int_t *     incx,
         double       *         dy,
         const blas_int_t *     incy);

// compute y = y + a * x
CFUNCDECL
void
daxpy_ ( const blas_int_t *     n,
         const double *         da,
         const double *         dx,
         const blas_int_t *     incx,
         double *               dy,
         const blas_int_t *     incy );

// compute sum of absolut values of x
CFUNCDECL
double
dasum_  ( const blas_int_t *    n,
          const double *        dx,
          const blas_int_t *    incx );

// return euclidean norm of x
CFUNCDECL
double
dnrm2_  ( const blas_int_t *    n,
          const double *        dx,
          const blas_int_t *    incx );

// scale x by a
CFUNCDECL
void
dscal_ ( const blas_int_t *     n,
         const double *         da,
         double *               dx,
         const blas_int_t *     incx );

// interchange x and y
CFUNCDECL
void
dswap_ ( const blas_int_t *     n,
         double *               dx,
         const blas_int_t *     incx,
         double *               dy,
         const blas_int_t *     incy );

// finds index with element having maximal absolut value
CFUNCDECL
blas_int_t
idamax_ ( const blas_int_t *    n,
          const double   *,
          const blas_int_t * );
    
// y = alpha A x + beta y
CFUNCDECL
void
dgemv_ ( const char *           trans,
         const blas_int_t *     M,
         const blas_int_t *     N,
         const double *         alpha,
         const double *         A,
         const blas_int_t *     lda,
         const double *         dx,
         const blas_int_t *     incx,
         const double *         beta,
         double *               dy,
         const blas_int_t *     incy );

// compute c = alpha * a * b + beta * c
CFUNCDECL
void
dgemm_ ( const char *           transa,
         const char *           transb,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const blas_int_t *     k,
         const double *         alpha,
         const double *         a,
         const blas_int_t *     lda,
         const double *         b,
         const blas_int_t *     ldb,
         const double *         beta,
         double *               c,
         const blas_int_t *     ldc );

// computes y = op(A)*x for upper/lower triangular A (y overwrites x)
void
dtrmv_ ( const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const double *         A,
         const blas_int_t *     ldA,
         double *               x,
         const blas_int_t *     incx );

// computes B = alpha * op(A) * B or B = alpha * B * op(A)
CFUNCDECL
void
dtrmm_ ( const char *           side,
         const char *           uplo,
         const char *           transa,
         const char *           diag,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const double *         alpha,
         const double *         A,
         const blas_int_t *     lda,
         double *               B,
         const blas_int_t *     ldb );

// solves A x = b or A^T x = b with triangular A
CFUNCDECL
void
dtrsv_ ( const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const double *         A,
         const blas_int_t *     lda,
         double *               x,
         const blas_int_t *     incx );
    
// solves A X = B or A^T X = B with triangular A
CFUNCDECL
void
dtrsm_ ( const char *           side,
         const char *           uplo,
         const char *           trans,
         const char *           diag,
         const blas_int_t *     n,
         const blas_int_t *     m,
         const double *         alpha,
         const double *         A,
         const blas_int_t *     lda,
         double *               B,
         const blas_int_t *     ldb );

// A = alpha * x * y^T + A, A \in \R^{m x n}
CFUNCDECL
void
dger_ ( const blas_int_t *      m,
        const blas_int_t *      n,
        const double *          alpha,
        const double *          x,
        const blas_int_t *      incx,
        const double *          y,
        const blas_int_t *      incy,
        double *                A,
        const blas_int_t *      lda );

///////////////////////////////////////////////////////////////////
//
// complex-valued functions (single precision)
//
///////////////////////////////////////////////////////////////////

// set (part of) A to alpha/beta
CFUNCDECL
void
claset_  ( const char *                 UPLO,
           const blas_int_t *           M,
           const blas_int_t *           N,
           const Complex<float> *       ALPHA,
           const Complex<float> *       BETA,
           Complex<float> *             A,
           const blas_int_t *           LDA   );

// copy (part of) A to B
CFUNCDECL
void
clacpy_  ( const char *                 UPLO,
           const blas_int_t *           M,
           const blas_int_t *           N,
           const Complex<float> *       A,
           const blas_int_t *           LDA,
           Complex<float> *             B,
           const blas_int_t *           LDB  );

// compute dotproduct between x and y
CFUNCDECL
void
xcdotu_ ( const blas_int_t *            n,
          const Complex<float> *        dx,
          const blas_int_t *            incx,
          const Complex<float> *        dy,
          const blas_int_t *            incy,
          Complex<float> *              retval );
CFUNCDECL
void
xcdotc_ ( const blas_int_t *            n,
          const Complex<float> *        dx,
          const blas_int_t *            incx,
          const Complex<float> *        dy,
          const blas_int_t *            incy,
          Complex<float> *              retval );

// copy vector x into vector y
CFUNCDECL
void
ccopy_ ( const blas_int_t *             n,
         const Complex<float> *         dx,
         const blas_int_t *             incx,
         Complex<float> *               dy,
         const blas_int_t *             incy);

// compute y = y + a * x
CFUNCDECL
void
caxpy_ ( const blas_int_t *             n,
         const Complex<float> *         da,
         const Complex<float> *         dx,
         const blas_int_t *             incx,
         Complex<float> *               dy,
         const blas_int_t *             incy );

// compute euclidean norm
CFUNCDECL
float
scnrm2_ ( const blas_int_t *            n,
          const Complex<float> *        x,
          const blas_int_t *            incx );

// scale x by a
CFUNCDECL
void
cscal_  ( const blas_int_t *            n,
          const Complex<float> *        da,
          Complex<float> *              dx,
          const blas_int_t *            incx );
CFUNCDECL
void
csscal_ ( const blas_int_t *            n,
          const float *                 da,
          Complex<float> *              dx,
          const blas_int_t *            incx );

// interchange x and y
CFUNCDECL
void
cswap_ ( const blas_int_t *     n,
         Complex<float> *       dx,
         const blas_int_t *     incx,
         Complex<float> *       dy,
         const blas_int_t *     incy );

// finds index with element having maximal absolut value
CFUNCDECL
blas_int_t
icamax_ ( const blas_int_t *    n,
          const Complex<float> *,
          const blas_int_t * );

// y = alpha A x + beta y
CFUNCDECL
void
cgemv_ ( const char *                   trans,
         const blas_int_t *             M,
         const blas_int_t *             N,
         const Complex<float> *         alpha,
         const Complex<float> *         A,
         const blas_int_t *             lda,
         const Complex<float> *         dx,
         const blas_int_t *             incx,
         const Complex<float> *         beta,
         Complex<float> *               dy,
         const blas_int_t *             incy );

// compute c = alpha * a * b + beta * c
CFUNCDECL
void
cgemm_ ( const char *                   transa,
         const char *                   transb,
         const blas_int_t *             n,
         const blas_int_t *             m,
         const blas_int_t *             k,
         const Complex<float> *         alpha,
         const Complex<float> *         a,
         const blas_int_t *             lda,
         const Complex<float> *         b,
         const blas_int_t *             ldb,
         const Complex<float> *         beta,
         Complex<float> *               c,
         const blas_int_t *             ldc );

// computes y = op(A)*x for upper/lower triangular A (y overwrites x)
void
ctrmv_ ( const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const Complex<float> *         A,
         const blas_int_t *             ldA,
         Complex<float> *               x,
         const blas_int_t *             incx );

// computes B = alpha * op(A) * B or B = alpha * B * op(A)
CFUNCDECL
void
ctrmm_ ( const char *                   side,
         const char *                   uplo,
         const char *                   transa,
         const char *                   diag,
         const blas_int_t *             m,
         const blas_int_t *             n,
         const Complex<float> *         alpha,
         const Complex<float> *         a,
         const blas_int_t *             lda,
         Complex<float> *               b,
         const blas_int_t *             ldb );

// solves A x = b or A^T x = b with triangular A
CFUNCDECL
void
ctrsv_ ( const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const Complex<float> *         a,
         const blas_int_t *             lda,
         Complex<float> *               x,
         const blas_int_t *             incx );

// solves A X = B or A^T X = B with triangular A
CFUNCDECL
void
ctrsm_ ( const char *                   side,
         const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const blas_int_t *             m,
         const Complex< float > *       alpha,
         const Complex< float > *       A,
         const blas_int_t *             lda,
         Complex< float > *             B,
         const blas_int_t *             ldb );

// rank-1 update: A = alpha * x * y^T + A
CFUNCDECL
void
cgeru_ ( const blas_int_t *             m,
         const blas_int_t *             n,
         const Complex<float> *         alpha,
         const Complex<float> *         x,
         const blas_int_t *             incx,
         const Complex<float> *         y,
         const blas_int_t *             incy,
         Complex<float> *               A,
         const blas_int_t *             lda );

// rank-1 update: A = alpha * x * conj(y^T) + A
CFUNCDECL
void
cgerc_ ( const blas_int_t *             m,
         const blas_int_t *             n,
         const Complex<float> *         alpha,
         const Complex<float> *         x,
         const blas_int_t *             incx,
         const Complex<float> *         y,
         const blas_int_t *             incy,
         Complex<float> *               A,
         const blas_int_t *             lda );

///////////////////////////////////////////////////////////////////
//
// complex-valued functions (double precision)
//
///////////////////////////////////////////////////////////////////

// set (part of) A to alpha/beta
CFUNCDECL
void
zlaset_  ( const char *                 UPLO,
           const blas_int_t *           M,
           const blas_int_t *           N,
           const Complex<double> *      ALPHA,
           const Complex<double> *      BETA,
           Complex<double> *            A,
           const blas_int_t *           LDA   );

// copy (part of) A to B
CFUNCDECL
void
zlacpy_  ( const char *                 UPLO,
           const blas_int_t *           M,
           const blas_int_t *           N,
           const Complex<double> *      A,
           const blas_int_t *           LDA,
           Complex<double> *            B,
           const blas_int_t *           LDB  );

// compute dotproduct between x and y
CFUNCDECL
void
xzdotu_ ( const blas_int_t *            n,
          const Complex<double> *       dx,
          const blas_int_t *            incx,
          const Complex<double> *       dy,
          const blas_int_t *            incy,
          Complex<double> *             retval );
CFUNCDECL
void
xzdotc_ ( const blas_int_t *            n,
          const Complex<double> *       dx,
          const blas_int_t *            incx,
          const Complex<double> *       dy,
          const blas_int_t *            incy,
          Complex<double> *             retval );

// copy vector x into vector y
CFUNCDECL
void
zcopy_ ( const blas_int_t *             n,
         const Complex<double> *        dx,
         const blas_int_t *             incx,
         Complex<double>       *        dy,
         const blas_int_t *             incy);

// compute y = y + a * x
CFUNCDECL
void
zaxpy_ ( const blas_int_t *             n,
         const Complex<double> *        da,
         const Complex<double> *        dx,
         const blas_int_t *             incx,
         Complex<double> *              dy,
         const blas_int_t *             incy );

// compute euclidean norm
CFUNCDECL
double
dznrm2_ ( const blas_int_t *            n,
          const Complex<double> *       x,
          const blas_int_t *            incx );

// scale x by a
CFUNCDECL
void
zscal_  ( const blas_int_t *            n,
          const Complex<double> *       da,
          Complex<double> *             dx,
          const blas_int_t *            incx );
CFUNCDECL
void
zdscal_ ( const blas_int_t *            n,
          const double *                da,
          Complex<double> *             dx,
          const blas_int_t *            incx );

// interchange x and y
CFUNCDECL
void
zswap_ ( const blas_int_t *             n,
         Complex<double> *              dx,
         const blas_int_t *             incx,
         Complex<double> *              dy,
         const blas_int_t *             incy );

// finds index with element having maximal absolut value
CFUNCDECL
blas_int_t
izamax_ ( const blas_int_t *            n,
          const Complex<double> *,
          const blas_int_t * );
    
// y = alpha A x + beta y
CFUNCDECL
void
zgemv_ ( const char *                   trans,
         const blas_int_t *             M,
         const blas_int_t *             N,
         const Complex<double> *        alpha,
         const Complex<double> *        A,
         const blas_int_t *             lda,
         const Complex<double> *        dx,
         const blas_int_t *             incx,
         const Complex<double> *        beta,
         Complex<double> *              dy,
         const blas_int_t *             incy );

// compute c = alpha * a * b + beta * c
CFUNCDECL
void
zgemm_ ( const char *                   transa,
         const char *                   transb,
         const blas_int_t *             n,
         const blas_int_t *             m,
         const blas_int_t *             k,
         const Complex<double> *        alpha,
         const Complex<double> *        a,
         const blas_int_t *             lda,
         const Complex<double> *        b,
         const blas_int_t *             ldb,
         const Complex<double> *        beta,
         Complex<double> *              c,
         const blas_int_t *             ldc );

// computes y = op(A)*x for upper/lower triangular A (y overwrites x)
void
ztrmv_ ( const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const Complex<double> *        A,
         const blas_int_t *             ldA,
         Complex<double> *              x,
         const blas_int_t *             incx );

// computes B = alpha * op(A) * B or B = alpha * B * op(A)
CFUNCDECL
void
ztrmm_ ( const char *                   side,
         const char *                   uplo,
         const char *                   transa,
         const char *                   diag,
         const blas_int_t *             m,
         const blas_int_t *             n,
         const Complex<double> *        alpha,
         const Complex<double> *        a,
         const blas_int_t *             lda,
         Complex<double> *              b,
         const blas_int_t *             ldb );

// solves A x = b or A^T x = b with triangular A
CFUNCDECL
void
ztrsv_ ( const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const Complex<double> *        a,
         const blas_int_t *             lda,
         Complex<double> *              x,
         const blas_int_t *             incx );

// solves A X = B or A^T X = B with triangular A
CFUNCDECL
void
ztrsm_ ( const char *                   side,
         const char *                   uplo,
         const char *                   trans,
         const char *                   diag,
         const blas_int_t *             n,
         const blas_int_t *             m,
         const Complex< double > *      alpha,
         const Complex< double > *      A,
         const blas_int_t *             lda,
         Complex< double > *            B,
         const blas_int_t *             ldb );

// rank-1 update: A = alpha * x * y^T + A
CFUNCDECL
void
zgeru_ ( const blas_int_t *             m,
         const blas_int_t *             n,
         const Complex<double> *        alpha,
         const Complex<double> *        x,
         const blas_int_t *             incx,
         const Complex<double> *        y,
         const blas_int_t *             incy,
         Complex<double> *              A,
         const blas_int_t *             lda );

// rank-1 update: A = alpha * x * conj(y^T) + A
CFUNCDECL
void
zgerc_ ( const blas_int_t *             m,
         const blas_int_t *             n,
         const Complex<double> *        alpha,
         const Complex<double> *        x,
         const blas_int_t *             incx,
         const Complex<double> *        y,
         const blas_int_t *             incy,
         Complex<double> *              A,
         const blas_int_t *             lda );
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// definition of external Lapack functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

extern "C" {
    
//////////////////////////////////////////////////////////////
//
// real-valued functions (single precision)
//
//////////////////////////////////////////////////////////////

// solve linear system of equations
CFUNCDECL
void
sgesv_   ( const blas_int_t *   n,
           const blas_int_t *   nrhs,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           float *              B,
           const blas_int_t *   ldb,
           blas_int_t *         info );

// compute inverse of triangular system
CFUNCDECL
void
strtri_  ( const char *         uplo,
           const char *         diag,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         info );

// compute eigenvalues and eigenvectors of a tridiagonal, symmetric matrix
CFUNCDECL
void
sstev_   ( const char *         jobz,
           const blas_int_t *   n,
           float *              D,
           float *              E,
           float *              Z,
           const blas_int_t *   ldz,
           float *              work,
           blas_int_t *         info );

// compute eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
ssyev_   ( const char *         jobz,
           const char *         uplo,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           float *              w,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );
    
// compute selected eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void ssyevx_ ( const char * jobz, const char * range, const char * uplo,
               const blas_int_t * n, float * A, const blas_int_t * ldA,
               const float * vl, const float * vu, const blas_int_t * il, const blas_int_t * iu,
               const float * abstol, blas_int_t * m, float * W, float * Z, const blas_int_t * ldZ,
               float * work, const blas_int_t * lwork, blas_int_t * iwork, blas_int_t * ifail,
               blas_int_t * info );
    
// compute singular-value-decomposition
CFUNCDECL
void
sgesvd_  ( const char *         jobu,
           const char *         jobv,
           const blas_int_t *   n,
           const blas_int_t *   m,
           float *              A,
           const blas_int_t *   lda,
           float *              S,
           float *              U,
           const blas_int_t *   ldu,
           float *              V,
           const blas_int_t *   ldv,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sgesdd_  ( const char *         job,
           const blas_int_t *   n,
           const blas_int_t *   m,
           float *              A,
           const blas_int_t *   lda,
           float *              S,
           float *              U,
           const blas_int_t *   ldu,
           float *              VT,
           const blas_int_t *   ldvt,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         info );

CFUNCDECL
void
sgesvj_  ( const char *         joba,
           const char *         jobu,
           const char *         jobv,
           const blas_int_t *   m,
           const blas_int_t *   n,
           float *              a,
           const blas_int_t *   lda,
           float *              sva,
           const blas_int_t *   mv,
           float *              v,
           const blas_int_t *   ldv,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sgejsv_  ( const char *         joba,
           const char *         jobu,
           const char *         jobv,
           const char *         jobr,
           const char *         jobt,
           const char *         jobp,
           const blas_int_t *   m,
           const blas_int_t *   n,
           float *              a,
           const blas_int_t *   lda,
           float *              sva,
           float *              u,
           const blas_int_t *   ldu,
           float *              v,
           const blas_int_t *   ldv,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         info );

// compute QR-factorisation
CFUNCDECL
void
sgeqrf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           float *              tau,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sorgqr_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           const blas_int_t *   k,
           float *              A,
           const blas_int_t *   lda,
           float *              tau,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
sgeqp3_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         jpvt,
           float *              tau,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

#if HAS_GEQP3_TRUNC == 1

CFUNCDECL
void
sgeqp3trunc_  ( const blas_int_t *   m,
                const blas_int_t *   n,
                float *              A,
                const blas_int_t *   lda,
                blas_int_t *         jpvt,
                float *              tau,
                blas_int_t *         ntrunc,
                float *              atrunc,
                float *              rtrunc,
                float *              work,
                const blas_int_t *   lwork,
                blas_int_t *         info );

#endif

// compute LU factorisation of given matrix A
CFUNCDECL
void
sgetrf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           blas_int_t *         info );

// compute inverse of A (using result from getrf)
CFUNCDECL
void
sgetri_  ( const blas_int_t *   n,
           float *              A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           float *              work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

// determine machine parameters
CFUNCDECL
float
slamch_  ( char *               cmach );

// generate plane-rotation
CFUNCDECL
blas_int_t
slartg_  ( float *              f,
           float *              g,
           float *              cs,
           float *              sn,
           float *              r );

// compute householder reflection
CFUNCDECL
void
slarfg_ ( const blas_int_t *      n,
          const float *           alpha,
          const float *           x,
          const blas_int_t *      incx,
          const float *           tau );

// apply householder reflection
CFUNCDECL
void
slarf_  ( const char *            side,
          const blas_int_t *      n,
          const blas_int_t *      m,
          const float *           V,
          const blas_int_t *      incv,
          const float *           tau,
          float *                 C,
          const blas_int_t *      ldc,
          const float *           work );

//////////////////////////////////////////////////////////////
//
// real-valued functions (double precision)
//
//////////////////////////////////////////////////////////////

// solve linear system of equations
CFUNCDECL
void
dgesv_   ( const blas_int_t *   n,
           const blas_int_t *   nrhs,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           double *             B,
           const blas_int_t *   ldb,
           blas_int_t *         info );

// compute inverse of triangular system
CFUNCDECL
void
dtrtri_  ( const char *         uplo,
           const char *         diag,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         info );
    
// compute eigenvalues and eigenvectors of a tridiagonal, symmetric matrix
CFUNCDECL
void
dstev_   ( const char *         jobz,
           const blas_int_t *   n,
           double *             D,
           double *             E,
           double *             Z,
           const blas_int_t *   ldz,
           double *             work,
           blas_int_t *         info );
    
// compute eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
dsyev_   ( const char *         jobz,
           const char *         uplo,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           double *             w,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );
    
// compute selected eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
dsyevx_  ( const char *         jobz,
           const char *         range,
           const char *         uplo,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   ldA,
           const double *       vl,
           const double *       vu,
           const blas_int_t *   il,
           const blas_int_t *   iu,
           const double *       abstol,
           blas_int_t *         m,
           double *             W,
           double *             Z,
           const blas_int_t *   ldZ,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         ifail,
           blas_int_t *         info );
    
// compute singular-value-decomposition
CFUNCDECL
void
dgesvd_  ( const char *         jobu,
           const char *         jobv,
           const blas_int_t *   n,
           const blas_int_t *   m,
           double *             A,
           const blas_int_t *   lda,
           double *             S,
           double *             U,
           const blas_int_t *   ldu,
           double *             V,
           const blas_int_t *   ldv,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
dgesdd_  ( const char *         jobz,
           const blas_int_t *   n,
           const blas_int_t *   m,
           double *             A,
           const blas_int_t *   lda,
           double *             S,
           double *             U,
           const blas_int_t *   ldu,
           double *             VT,
           const blas_int_t *   ldvt,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         info );

CFUNCDECL
void
dgejsv_  ( const char *         joba,
           const char *         jobu,
           const char *         jobv,
           const char *         jobr,
           const char *         jobt,
           const char *         jobp,
           const blas_int_t *   m,
           const blas_int_t *   n,
           double *             a,
           const blas_int_t *   lda,
           double *             sva,
           double *             u,
           const blas_int_t *   ldu,
           double *             v,
           const blas_int_t *   ldv,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         iwork,
           blas_int_t *         info );

CFUNCDECL
void
dgesvj_  ( const char *         joba,
           const char *         jobu,
           const char *         jobv,
           const blas_int_t *   m,
           const blas_int_t *   n,
           double *             a,
           const blas_int_t *   lda,
           double *             sva,
           const blas_int_t *   mv,
           double *             v,
           const blas_int_t *   ldv,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

// compute QR-factorisation
CFUNCDECL
void
dgeqrf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           double *             tau,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
dorgqr_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           const blas_int_t *   k,
           double *             A,
           const blas_int_t *   lda,
           double *             tau,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

CFUNCDECL
void
dgeqp3_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         jpvt,
           double *             tau,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );

#if HAS_GEQP3_TRUNC == 1

CFUNCDECL
void
dgeqp3trunc_  ( const blas_int_t *   m,
                const blas_int_t *   n,
                double *             A,
                const blas_int_t *   lda,
                blas_int_t *         jpvt,
                double *             tau,
                blas_int_t *         ntrunc,
                double *             atrunc,
                double *             rtrunc,
                double *             work,
                const blas_int_t *   lwork,
                blas_int_t *         info );

#endif

// compute LU factorisation of given matrix A
CFUNCDECL
void
dgetrf_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           blas_int_t *         info );
    
// compute inverse of A (using result from getrf)
CFUNCDECL
void
dgetri_  ( const blas_int_t *   n,
           double *             A,
           const blas_int_t *   lda,
           blas_int_t *         ipiv,
           double *             work,
           const blas_int_t *   lwork,
           blas_int_t *         info );
    
// determine machine parameters
CFUNCDECL
double
dlamch_  ( char *               cmach );

// generate plane-rotation
CFUNCDECL
blas_int_t
dlartg_  ( double *             f,
           double *             g,
           double *             cs,
           double *             sn,
           double *             r );

// compute householder reflection
CFUNCDECL
void
dlarfg_ ( const blas_int_t *      n,
          const double *          alpha,
          const double *          x,
          const blas_int_t *      incx,
          const double *          tau );

// apply householder reflection
CFUNCDECL
void
dlarf_  ( const char *            side,
          const blas_int_t *      n,
          const blas_int_t *      m,
          const double *          V,
          const blas_int_t *      incv,
          const double *          tau,
          double *                C,
          const blas_int_t *      ldc,
          const double *          work );

//////////////////////////////////////////////////////////////
//
// complex-valued functions (single precision)
//
//////////////////////////////////////////////////////////////

// solve linear system of equations
CFUNCDECL
void
cgesv_   ( const blas_int_t *           n,
           const blas_int_t *           nrhs,
           Complex<float> *             a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           Complex<float> *             b,
           const blas_int_t *           ldb,
           blas_int_t *                 info );

// compute inverse of triangular system
CFUNCDECL
void
ctrtri_  ( const char *                 uplo,
           const char *                 diag,
           const blas_int_t *           n,
           Complex< float > *           A,
           const blas_int_t *           lda,
           blas_int_t *                 info );
    
// compute eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
cheev_   ( const char *                 jobz,
           const char *                 uplo,
           const blas_int_t *           n,
           Complex<float> *             A,
           const blas_int_t *           lda,
           float *                      w,
           Complex<float> *             work,
           const blas_int_t *           lwork,
           float *                      rwork,
           blas_int_t *                 info );
    
// compute singular-value-decomposition
CFUNCDECL
void
cgesvd_  ( const char *                 jobu,
           const char *                 jobv,
           const blas_int_t *           n,
           const blas_int_t *           m,
           Complex<float> *             A,
           const blas_int_t *           lda,
           float *                      S,
           Complex<float> *             U,
           const blas_int_t *           ldu,
           Complex<float> *             V,
           const blas_int_t *           ldv,
           Complex<float> *             work,
           const blas_int_t *           lwork,
           float *                      rwork,
           blas_int_t *                 info );

CFUNCDECL
void
cgesdd_  ( const char *                 job,
           const blas_int_t *           n,
           const blas_int_t *           m,
           Complex<float> *             A,
           const blas_int_t *           lda,
           float *                      S,
           Complex<float> *             U,
           const blas_int_t *           ldu,
           Complex<float> *             VT,
           const blas_int_t *           ldvt,
           Complex<float> *             work,
           const blas_int_t *           lwork,
           float *                      rwork,
           const blas_int_t *           iwork,
           blas_int_t *                 info );

CFUNCDECL
void
cgesvj_  ( const char *                 joba,
           const char *                 jobu,
           const char *                 jobv,
           const blas_int_t *           n,
           const blas_int_t *           m,
           Complex<float> *             A,
           const blas_int_t *           lda,
           float *                      S,
           const blas_int_t *           mv,
           Complex<float> *             V,
           const blas_int_t *           ldv,
           Complex<float> *             cwork,
           const blas_int_t *           lwork,
           float *                      rwork,
           const blas_int_t *           lrwork,
           blas_int_t *                 info );

// compute QR-factorisation
CFUNCDECL
void
cgeqrf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           Complex<float> *             A,
           const blas_int_t *           lda,
           Complex<float> *             tau,
           Complex<float> *             work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
cungqr_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           const blas_int_t *           k,
           Complex<float> *             a,
           const blas_int_t *           lda,
           Complex<float> *             tau,
           Complex<float> *             work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
cgeqp3_  ( const blas_int_t *   m,
           const blas_int_t *   n,
           Complex<float> *     A,
           const blas_int_t *   lda,
           blas_int_t *         jpvt,
           Complex<float> *     tau,
           Complex<float> *     work,
           const blas_int_t *   lwork,
           float *              rwork,
           blas_int_t *         info );

#if HAS_GEQP3_TRUNC == 1

CFUNCDECL
void
cgeqp3trunc_  ( const blas_int_t *   m,
                const blas_int_t *   n,
                Complex<float> *     A,
                const blas_int_t *   lda,
                blas_int_t *         jpvt,
                Complex<float> *     tau,
                blas_int_t *         ntrunc,
                float *              atrunc,
                float *              rtrunc,
                Complex<float> *     work,
                const blas_int_t *   lwork,
                float *              rwork,
                blas_int_t *         info );

#endif

// compute LU factorisation of given matrix A
CFUNCDECL
void
cgetrf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           Complex<float> *             a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           blas_int_t *                 info);

// compute inverse of A (using result from getrf)
CFUNCDECL
void
cgetri_  ( const blas_int_t *           n,
           Complex<float> *             a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           Complex<float> *             work,
           const blas_int_t *           lwork,
           blas_int_t *                 info);

// compute householder reflection
CFUNCDECL
void
clarfg_ ( const blas_int_t *            n,
          const Complex< float > *      alpha,
          const Complex< float > *      x,
          const blas_int_t *            incx,
          const Complex< float > *      tau );

// apply householder reflection
CFUNCDECL
void
clarf_  ( const char *                  side,
          const blas_int_t *            n,
          const blas_int_t *            m,
          const Complex< float > *      V,
          const blas_int_t *            incv,
          const Complex< float > *      tau,
          Complex< float > *            C,
          const blas_int_t *            ldc,
          const Complex< float > *      work );

//////////////////////////////////////////////////////////////
//
// complex-valued functions (double precision)
//
//////////////////////////////////////////////////////////////

// solve linear system of equations
CFUNCDECL
void
zgesv_   ( const blas_int_t *           n,
           const blas_int_t *           nrhs,
           Complex<double> *            a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           Complex<double> *            b,
           const blas_int_t *           ldb,
           blas_int_t *                 info );

// compute inverse of triangular system
CFUNCDECL
void
ztrtri_  ( const char *                 uplo,
           const char *                 diag,
           const blas_int_t *           n,
           Complex< double > *          A,
           const blas_int_t *           lda,
           blas_int_t *                 info );
    
// compute eigenvalues and eigenvectors of a symmetric matrix
CFUNCDECL
void
zheev_   ( const char *                 jobz,
           const char *                 uplo,
           const blas_int_t *           n,
           Complex<double> *            A,
           const blas_int_t *           lda,
           double *                     w,
           Complex<double> *            work,
           const blas_int_t *           lwork,
           double *                     rwork,
           blas_int_t *                 info );
    
// compute singular-value-decomposition
CFUNCDECL
void
zgesvd_  ( const char *                 jobu,
           const char *                 jobv,
           const blas_int_t *           n,
           const blas_int_t *           m,
           Complex<double> *            A,
           const blas_int_t *           lda,
           double *                     S,
           Complex<double> *            U,
           const blas_int_t *           ldu,
           Complex<double> *            V,
           const blas_int_t *           ldv,
           Complex<double> *            work,
           const blas_int_t *           lwork,
           double *                     rwork,
           blas_int_t *                 info );

CFUNCDECL
void
zgesdd_  ( const char *                 job,
           const blas_int_t *           n,
           const blas_int_t *           m,
           Complex<double> *            A,
           const blas_int_t *           lda,
           double *                     S,
           Complex<double> *            U,
           const blas_int_t *           ldu,
           Complex<double> *            VT,
           const blas_int_t *           ldvt,
           Complex<double> *            work,
           const blas_int_t *           lwork,
           double *                     rwork,
           const blas_int_t *           iwork,
           blas_int_t *                 info );

CFUNCDECL
void
zgesvj_  ( const char *                 joba,
           const char *                 jobu,
           const char *                 jobv,
           const blas_int_t *           n,
           const blas_int_t *           m,
           Complex<double> *            A,
           const blas_int_t *           lda,
           double *                     S,
           const blas_int_t *           mv,
           Complex<double> *            V,
           const blas_int_t *           ldv,
           Complex<double> *            cwork,
           const blas_int_t *           lwork,
           double *                     rwork,
           const blas_int_t *           lrwork,
           blas_int_t *                 info );

// compute QR-factorisation
CFUNCDECL
void
zgeqrf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           Complex<double> *            A,
           const blas_int_t *           lda,
           Complex<double> *            tau,
           Complex<double> *            work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
zungqr_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           const blas_int_t *           k,
           Complex<double> *            a,
           const blas_int_t *           lda,
           Complex<double> *            tau,
           Complex<double> *            work,
           const blas_int_t *           lwork,
           blas_int_t *                 info );

CFUNCDECL
void
zgeqp3_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           Complex<double> *            A,
           const blas_int_t *           lda,
           blas_int_t *                 jpvt,
           Complex<double> *            tau,
           Complex<double> *            work,
           const blas_int_t *           lwork,
           double *                     rwork,
           blas_int_t *                 info );

#if HAS_GEQP3_TRUNC == 1

CFUNCDECL
void
zgeqp3trunc_  ( const blas_int_t *      m,
                const blas_int_t *      n,
                Complex<double> *       A,
                const blas_int_t *      lda,
                blas_int_t *            jpvt,
                Complex<double> *       tau,
                blas_int_t *            ntrunc,
                double *                atrunc,
                double *                rtrunc,
                Complex<double> *       work,
                const blas_int_t *      lwork,
                double *                rwork,
                blas_int_t *            info );

#endif

// compute LU factorisation of given matrix A
CFUNCDECL
void
zgetrf_  ( const blas_int_t *           m,
           const blas_int_t *           n,
           Complex<double> *            a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           blas_int_t *                 info);

// compute inverse of A (using result from getrf)
CFUNCDECL
void
zgetri_  ( const blas_int_t *           n,
           Complex<double> *            a,
           const blas_int_t *           lda,
           blas_int_t *                 ipiv,
           Complex<double> *            work,
           const blas_int_t *           lwork,
           blas_int_t *                 info);

// compute householder reflection
CFUNCDECL
void
zlarfg_ ( const blas_int_t *            n,
          const Complex< double > *     alpha,
          const Complex< double > *     x,
          const blas_int_t *            incx,
          const Complex< double > *     tau );

// apply householder reflection
CFUNCDECL
void
zlarf_  ( const char *                  side,
          const blas_int_t *            n,
          const blas_int_t *            m,
          const Complex< double > *     V,
          const blas_int_t *            incv,
          const Complex< double > *     tau,
          Complex< double > *           C,
          const blas_int_t *            ldc,
          const Complex< double > *     work );

//////////////////////////////////////////////////////////////
//
// misc. helpers
//
//////////////////////////////////////////////////////////////

// return problem dependent parameters
blas_int_t
ilaenv_ ( const blas_int_t *    ispec,
          const char *          name,
          const char *          opts,
          const blas_int_t *    n1,
          const blas_int_t *    n2,
          const blas_int_t *    n3,
          const blas_int_t *    n4 );

}// extern "C"



namespace BLAS
{

//! @cond

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// template wrappers for BLAS functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// *laset
//
template < typename T >
void laset ( const char        uplo,
             const blas_int_t  m,
             const blas_int_t  n,
             const T           alpha,
             const T           beta,
             T *               a,
             const blas_int_t  lda );

#define HLIB_LASET_FUNC( type, func )                                  \
    template <> inline void laset<type> ( const char       uplo,       \
                                          const blas_int_t m,          \
                                          const blas_int_t n,          \
                                          const type       alpha,      \
                                          const type       beta,       \
                                          type *           A,          \
                                          const blas_int_t lda ) {     \
        func( & uplo, & m, & n, & alpha, & beta, A, & lda ); }

HLIB_LASET_FUNC( float,             slaset_ )
HLIB_LASET_FUNC( double,            dlaset_ )
HLIB_LASET_FUNC( Complex< float >,  claset_ )
HLIB_LASET_FUNC( Complex< double >, zlaset_ )

#undef HLIB_LASET_FUNC

//
// *lacpy
//
template < typename T >
void lacpy ( const char        uplo,
             const blas_int_t  m,
             const blas_int_t  n,
             const T *         A,
             const blas_int_t  lda,
             T *               B,
             const blas_int_t  ldb );

#define HLIB_LACPY_FUNC( type, func )                       \
    template <> inline void lacpy<type> ( const char       uplo,   \
                                   const blas_int_t m,      \
                                   const blas_int_t n,      \
                                   const type *     A,      \
                                   const blas_int_t lda,    \
                                   type *           B,      \
                                   const blas_int_t ldb ) { \
        func( & uplo, & m, & n, A, & lda, B, & ldb ); }

HLIB_LACPY_FUNC( float,             slacpy_ )
HLIB_LACPY_FUNC( double,            dlacpy_ )
HLIB_LACPY_FUNC( Complex< float >,  clacpy_ )
HLIB_LACPY_FUNC( Complex< double >, zlacpy_ )

#undef HLIB_LACPY_FUNC

//
// *scal
//
template < typename T >
void scal ( const blas_int_t n,
            const T          alpha,
            T *              x,
            const blas_int_t incx );

#define HLIB_SCAL_FUNC( type, func )                                    \
    template <> inline void scal ( const blas_int_t n,                  \
                                   const type       alpha,              \
                                   type *           x,                  \
                                   const blas_int_t incx )              \
        { func( & n, & alpha, x, & incx ); }

HLIB_SCAL_FUNC( float,             sscal_ )
HLIB_SCAL_FUNC( double,            dscal_ )
HLIB_SCAL_FUNC( Complex< float >,  cscal_ )
HLIB_SCAL_FUNC( Complex< double >, zscal_ )

//
// *copy
//
template <typename T>
void copy ( const blas_int_t  n,
            const T *  x, const blas_int_t  incx,
            T *        y, const blas_int_t  incy );

#define HLIB_COPY_FUNC( type, func )                                    \
    template <> inline void copy ( const blas_int_t n, const type * x,  \
                                   const blas_int_t incx,               \
                                   type * y, const blas_int_t incy )    \
        { func( & n, x, & incx, y, & incy ); }
                                                       
HLIB_COPY_FUNC( float,             scopy_ )
HLIB_COPY_FUNC( double,            dcopy_ )
HLIB_COPY_FUNC( Complex< float >,  ccopy_ )
HLIB_COPY_FUNC( Complex< double >, zcopy_ )

//
// *swap
//
template <typename T>
void swap ( const blas_int_t  n,
            T *  x, const blas_int_t  incx,
            T *  y, const blas_int_t  incy );

#define HLIB_SWAP_FUNC( type, func )                                    \
    template <>  inline  void swap ( const blas_int_t n, type * x,      \
                                     const blas_int_t incx,             \
                                     type * y, const blas_int_t incy )  \
        { func( & n, x, & incx, y, & incy ); }
                                                       
HLIB_SWAP_FUNC( float,             sswap_ )
HLIB_SWAP_FUNC( double,            dswap_ )
HLIB_SWAP_FUNC( Complex< float >,  cswap_ )
HLIB_SWAP_FUNC( Complex< double >, zswap_ )

//
// i*amax
//
template <typename T>
blas_int_t max_idx ( const blas_int_t  n, T *  x, const blas_int_t  incx );

#define HLIB_MAX_IDX_FUNC( type, func )                                 \
    template <> inline blas_int_t max_idx ( const blas_int_t n,         \
                                            type * x, const blas_int_t incx ) \
    { return func( & n, x, & incx ); }
                                                       
HLIB_MAX_IDX_FUNC( float,             isamax_ )
HLIB_MAX_IDX_FUNC( double,            idamax_ )
HLIB_MAX_IDX_FUNC( Complex< float >,  icamax_ )
HLIB_MAX_IDX_FUNC( Complex< double >, izamax_ )

//
// *axpy
//
template <typename T>
void axpy ( const blas_int_t n,
            const T alpha, const T * x, const blas_int_t incx,
            T * y, const blas_int_t incy );

#define HLIB_AXPY_FUNC( type, func, flops )                             \
    template <> inline void axpy ( const blas_int_t n,                  \
                                   const type alpha, const type * x,    \
                                   const blas_int_t incx,               \
                                   type * y, const blas_int_t incy )    \
    {                                                                   \
        ADD_FLOPS( FLOPS_##flops( n ) );                                \
        func( & n, & alpha, x, & incx, y, & incy );                     \
    }

HLIB_AXPY_FUNC( float,             saxpy_, SAXPY )
HLIB_AXPY_FUNC( double,            daxpy_, DAXPY )
HLIB_AXPY_FUNC( Complex< float >,  caxpy_, CAXPY )
HLIB_AXPY_FUNC( Complex< double >, zaxpy_, ZAXPY )

//
// *dot(c)
//
template <typename T>
T dot ( const blas_int_t n,
        const T * x, const blas_int_t incx,
        const T * y, const blas_int_t incy );

template <> inline
float dot ( const blas_int_t  n,
            const float *  x, const blas_int_t  incx,
            const float *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_SDOT( n ) );
    return sdot_( & n, x, & incx, y, & incy );
}

template <> inline
double dot ( const blas_int_t  n,
             const double *  x, const blas_int_t  incx,
             const double *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_DDOT( n ) );
    return ddot_( & n, x, & incx, y, & incy );
}

template <> inline
Complex< float > dot ( const blas_int_t  n,
                       const Complex< float > *  x, const blas_int_t  incx,
                       const Complex< float > *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_CDOT( n ) );

    Complex< float >  res;
    
    #if USE_MKL == 1

    cblas_cdotc_sub( n, x, incx, y, incy, & res );
    
    #else
    
    xcdotc_( & n, x, & incx, y, & incy, & res );
    
    #endif
    
    return res;
}

template <> inline
Complex< double > dot ( const blas_int_t  n,
                        const Complex< double > *  x, const blas_int_t  incx,
                        const Complex< double > *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_ZDOT( n ) );

    Complex< double >  res;
    
    #if USE_MKL == 1

    cblas_zdotc_sub( n, x, incx, y, incy, & res );
    
    #else
    
    xzdotc_( & n, x, & incx, y, & incy, & res );
    
    #endif
    
    return res;
}

//
// *dotu
//
template <typename T>
T dotu ( const blas_int_t n,
         const T * x, const blas_int_t incx,
         const T * y, const blas_int_t incy )
{
    ADD_FLOPS( FLOPS_DDOT( n ) );
    
    return dot<T>( n, x, incx, y, incy );
}

template <> inline
Complex< float > dotu ( const blas_int_t  n,
                        const Complex< float > *  x, const blas_int_t  incx,
                        const Complex< float > *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_CDOT( n ) );
    
    Complex< float >  res;
    
    #if USE_MKL == 1

    cblas_cdotu_sub( n, x, incx, y, incy, & res );
    
    #else
    
    xcdotu_( & n, x, & incx, y, & incy, & res );
    
    #endif
    
    return res;
}

template <> inline
Complex< double > dotu ( const blas_int_t  n,
                         const Complex< double > *  x, const blas_int_t  incx,
                         const Complex< double > *  y, const blas_int_t  incy )
{
    ADD_FLOPS( FLOPS_ZDOT( n ) );
    
    Complex< double >  res;
    
    #if USE_MKL == 1

    cblas_zdotu_sub( n, x, incx, y, incy, & res );
    
    #else
    
    xzdotu_( & n, x, & incx, y, & incy, & res );
    
    #endif
    
    return res;
}

//
// *norm2
//
template <typename T>
typename real_type< T >::type_t
norm2 ( const blas_int_t n, const T * x, const blas_int_t incx );

#define HLIB_NORM2_FUNC( type, func, flops )                            \
    template <> inline                                                  \
    typename real_type< type >::type_t                                  \
    norm2 ( const blas_int_t n, const type * x, const blas_int_t incx ) \
    {                                                                   \
        ADD_FLOPS( FLOPS_##flops( n ) );                                \
        return func( & n, x, & incx );                                  \
    }

HLIB_NORM2_FUNC( float,             snrm2_,  SDOT )
HLIB_NORM2_FUNC( double,            dnrm2_,  DDOT )
HLIB_NORM2_FUNC( Complex< float >,  scnrm2_, CDOT )
HLIB_NORM2_FUNC( Complex< double >, dznrm2_, ZDOT )

//
// *ger
//
template <typename T>
void ger ( const blas_int_t n, const blas_int_t m, const T alpha,
           const T * x, const blas_int_t incx,
           const T * y, const blas_int_t incy,
           T * A, const blas_int_t ldA );

#define HLIB_GER_FUNC( type, func, flops )                              \
    template <> inline                                                  \
    void ger ( const blas_int_t n, const blas_int_t m, const type alpha, \
               const type * x, const blas_int_t incx,                   \
               const type * y, const blas_int_t incy,                   \
               type * A, const blas_int_t ldA )                         \
    {                                                                   \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        func( & n, & m, & alpha, x, & incx, y, & incy, A, & ldA );      \
    }

HLIB_GER_FUNC( float,             sger_,  SGER )
HLIB_GER_FUNC( double,            dger_ , DGER )
HLIB_GER_FUNC( Complex< float >,  cgerc_, CGER )
HLIB_GER_FUNC( Complex< double >, zgerc_, ZGER )

//
// *geru
//
template <typename T>
void geru ( const blas_int_t n, const blas_int_t m, const T alpha,
            const T * x, const blas_int_t incx,
            const T * y, const blas_int_t incy,
            T * A, const blas_int_t ldA );

#define HLIB_GERU_FUNC( type, func, flops )                             \
    template <> inline                                                  \
    void geru ( const blas_int_t n, const blas_int_t m, const type alpha, \
               const type * x, const blas_int_t incx,                   \
               const type * y, const blas_int_t incy,                   \
               type * A, const blas_int_t ldA )                         \
    {                                                                   \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        func( & n, & m, & alpha, x, & incx, y, & incy, A, & ldA );      \
    }

HLIB_GERU_FUNC( float,             sger_,  SGER )
HLIB_GERU_FUNC( double,            dger_,  DGER )
HLIB_GERU_FUNC( Complex< float >,  cgeru_, CGER )
HLIB_GERU_FUNC( Complex< double >, zgeru_, ZGER )

//
// *gemv
//
template <typename T>
void gemv ( const char trans, const blas_int_t n, const blas_int_t m, const T alpha,
            const T * A, const blas_int_t ldA,
            const T * x, const blas_int_t incx,
            const T beta, T * y, const blas_int_t incy );

#define HLIB_GEMV_FUNC( type, func, flops )                             \
    template <> inline                                                  \
    void gemv ( const char trans, const blas_int_t n, const blas_int_t m, const type alpha, \
                const type * A, const blas_int_t ldA,                   \
                const type * x, const blas_int_t incx,                  \
                const type beta, type * y, const blas_int_t incy ) {    \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & ltrans, & n, & m, & alpha, A, & ldA, x, & incx, & beta, y, & incy ); }

HLIB_GEMV_FUNC( float,             sgemv_, SGEMV )
HLIB_GEMV_FUNC( double,            dgemv_, DGEMV )
HLIB_GEMV_FUNC( Complex< float >,  cgemv_, CGEMV )
HLIB_GEMV_FUNC( Complex< double >, zgemv_, ZGEMV )

//
// *trmv
//
template <typename T>
void trmv ( const char uplo, const char trans, const char diag,
            const blas_int_t n, const T * A, const blas_int_t ldA,
            T * x, const blas_int_t incx );

#define HLIB_TRMV_FUNC( type, func, flops )                             \
    template <> inline                                                  \
    void trmv ( const char uplo, const char trans, const char diag,     \
                const blas_int_t n, const type * A, const blas_int_t ldA, \
                type * x, const blas_int_t incx ) {                     \
        ADD_FLOPS( FLOPS_##flops( n ) );                                \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & uplo, & ltrans, & diag, & n, A, & ldA, x, & incx ); }

HLIB_TRMV_FUNC( float,             strmv_, STRMV )
HLIB_TRMV_FUNC( double,            dtrmv_, DTRMV )
HLIB_TRMV_FUNC( Complex< float >,  ctrmv_, CTRMV )
HLIB_TRMV_FUNC( Complex< double >, ztrmv_, ZTRMV )

//
// *trsv
//
template <typename T>
void trsv ( const char uplo, const char trans, const char diag,
            const blas_int_t n, const T * A, const blas_int_t ldA,
            T * b, const blas_int_t incb );

#define HLIB_TRSV_FUNC( type, func, flops )                             \
    template <> inline                                                  \
    void trsv ( const char uplo, const char trans, const char diag,     \
                const blas_int_t n, const type * A, const blas_int_t ldA, \
                type * b, const blas_int_t incb ) {                     \
        ADD_FLOPS( FLOPS_##flops( n ) );                                \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & uplo, & ltrans, & diag, & n, A, & ldA, b, & incb ); }

HLIB_TRSV_FUNC( float,             strsv_, STRSV )
HLIB_TRSV_FUNC( double,            dtrsv_, DTRSV )
HLIB_TRSV_FUNC( Complex< float >,  ctrsv_, CTRSV )
HLIB_TRSV_FUNC( Complex< double >, ztrsv_, ZTRSV )

//
// *gemm
//
template <typename T>
void gemm ( const char transA, const char transB,
            const blas_int_t n, const blas_int_t m, const blas_int_t k,
            const T alpha,
            const T * A, const blas_int_t ldA,
            const T * B, const blas_int_t ldB,
            const T beta, T * C, const blas_int_t ldC );

#define HLIB_GEMM_FUNC( type, func, flops )                             \
    template <> inline                                                  \
    void gemm ( const char transA, const char transB,                   \
                const blas_int_t n, const blas_int_t m, const blas_int_t k, \
                const type alpha,                                       \
                const type * A, const blas_int_t ldA,                   \
                const type * B, const blas_int_t ldB,                   \
                const type beta, type * C, const blas_int_t ldC ) {     \
        ADD_FLOPS( FLOPS_##flops( n, m, k ) );                          \
        char  ltransA = transA;                                         \
        char  ltransB = transB;                                         \
        if ( ! is_complex_type< type >::value && ( transA == 'C' ) )    \
            ltransA = 'T';                                              \
        if ( ! is_complex_type< type >::value && ( transB == 'C' ) )    \
            ltransB = 'T';                                              \
        func( & ltransA, & ltransB, & n, & m, & k,                      \
              & alpha, A, & ldA, B, & ldB,                              \
              & beta, C, & ldC ); }

HLIB_GEMM_FUNC( float,             sgemm_, SGEMM )
HLIB_GEMM_FUNC( double,            dgemm_, DGEMM )
HLIB_GEMM_FUNC( Complex< float >,  cgemm_, CGEMM )
HLIB_GEMM_FUNC( Complex< double >, zgemm_, ZGEMM )

//
// *trmm
//
template <typename T>
void trmm ( const char side, const char uplo, const char trans, const char diag,
            const blas_int_t n, const blas_int_t m, const T alpha, const T * A,
            const blas_int_t ldA, T * B, const blas_int_t ldB );

#define HLIB_TRMM_FUNC( type, func, flops )                             \
    template <> inline                                                  \
    void trmm ( const char side, const char uplo, const char trans, const char diag, \
                const blas_int_t n, const blas_int_t m, const type alpha, const type * A, \
                const blas_int_t ldA, type * B, const blas_int_t ldB ) { \
        ADD_FLOPS( FLOPS_##flops( side, n, m ) );                       \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & side, & uplo, & ltrans, & diag, & n, & m, & alpha, A, & ldA, B, & ldB ); }

HLIB_TRMM_FUNC( float,             strmm_, STRMM )
HLIB_TRMM_FUNC( double,            dtrmm_, CTRMM )
HLIB_TRMM_FUNC( Complex< float >,  ctrmm_, DTRMM )
HLIB_TRMM_FUNC( Complex< double >, ztrmm_, ZTRMM )

//
// *trsm
//
template <typename T>
void trsm ( const char side, const char uplo, const char trans, const char diag,
            const blas_int_t n, const blas_int_t m, const T alpha, const T * A,
            const blas_int_t ldA, T * B, const blas_int_t ldB );

#define HLIB_TRSM_FUNC( type, func, flops )                             \
    template <> inline                                                  \
    void trsm ( const char side, const char uplo, const char trans, const char diag, \
                const blas_int_t n, const blas_int_t m, const type alpha, const type * A, \
                const blas_int_t ldA, type * B, const blas_int_t ldB ) { \
        ADD_FLOPS( FLOPS_##flops( side, n, m ) );                       \
        char  ltrans = trans;                                           \
        if ( ! is_complex_type< type >::value && ( trans == 'C' ) )     \
            ltrans = 'T';                                               \
        func( & side, & uplo, & ltrans, & diag, & n, & m, & alpha, A, & ldA, B, & ldB ); }

HLIB_TRSM_FUNC( float,             strsm_, STRSM )
HLIB_TRSM_FUNC( double,            dtrsm_, CTRSM )
HLIB_TRSM_FUNC( Complex< float >,  ctrsm_, DTRSM )
HLIB_TRSM_FUNC( Complex< double >, ztrsm_, ZTRSM )



//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//
// template wrappers for LAPACK functions
//
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////

//
// *gesv
//
template < typename T >
void gesv ( const blas_int_t  n,
            const blas_int_t  nrhs,
            T *               A,
            const blas_int_t  ldA,
            blas_int_t *      ipiv,
            T *               B,
            const blas_int_t  ldB,
            blas_int_t &      info );

#define HLIB_GESV_FUNC( type, func )                                \
    template <> inline void gesv<type> ( const blas_int_t  n,       \
                                         const blas_int_t  nrhs,    \
                                         type *            A,       \
                                         const blas_int_t  ldA,     \
                                         blas_int_t *      ipiv,    \
                                         type *            B,       \
                                         const blas_int_t  ldB,     \
                                         blas_int_t &      info ) { \
        info = 0;                                                   \
        func( & n, & nrhs, A, & ldA, ipiv, B, & ldB, & info ); }

HLIB_GESV_FUNC( float,             sgesv_ )
HLIB_GESV_FUNC( double,            dgesv_ )
HLIB_GESV_FUNC( Complex< float >,  cgesv_ )
HLIB_GESV_FUNC( Complex< double >, zgesv_ )

#undef HLIB_GESV_FUNC

//
// *trtri
//
template < typename T >
void trtri ( const char        uplo,
             const char        diag,
             const blas_int_t  n,
             T *               A,
             const blas_int_t  ldA,
             blas_int_t &      info );

#define HLIB_TRTRI_FUNC( type, func )                           \
    template <> inline  void trtri<type> ( const char        uplo,  \
                                           const char        diag,  \
                                           const blas_int_t  n,     \
                                           type *            A,     \
                                           const blas_int_t  ldA,   \
                                           blas_int_t &      info ) {   \
        info = 0;                                               \
        func( & uplo, & diag, & n, A, & ldA, & info ); }

HLIB_TRTRI_FUNC( float,             strtri_ )
HLIB_TRTRI_FUNC( double,            dtrtri_ )
HLIB_TRTRI_FUNC( Complex< float >,  ctrtri_ )
HLIB_TRTRI_FUNC( Complex< double >, ztrtri_ )

#undef HLIB_TRTRI_FUNC

//
// *getrf
//
template < typename T >
void getrf ( const blas_int_t  n,
             const blas_int_t  m,
             T *               A,
             const blas_int_t  ldA,
             blas_int_t *      ipiv,
             blas_int_t &      info );

#define HLIB_GETRF_FUNC( type, func, flops )                            \
    template <> inline void getrf<type> ( const blas_int_t  n,          \
                                          const blas_int_t  m,          \
                                          type *            A,          \
                                          const blas_int_t  ldA,        \
                                          blas_int_t *      ipiv,       \
                                          blas_int_t &      info ) {    \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, ipiv, & info ); }

HLIB_GETRF_FUNC( float,             sgetrf_, SGETRF )
HLIB_GETRF_FUNC( double,            dgetrf_, CGETRF )
HLIB_GETRF_FUNC( Complex< float >,  cgetrf_, DGETRF )
HLIB_GETRF_FUNC( Complex< double >, zgetrf_, ZGETRF )

#undef HLIB_GETRF_FUNC

//
// *getri
//
template <typename T>  
void getri ( const blas_int_t  n,
             T *               A,
             const blas_int_t  ldA,
             blas_int_t *      ipiv,
             T *               work,
             blas_int_t        lwork,
             blas_int_t &      info );

#define HLIB_GETRI_FUNC( type, func, flops )                           \
    template <> inline void getri<type> ( const blas_int_t  n,         \
                                          type *            A,         \
                                          const blas_int_t  ldA,       \
                                          blas_int_t *      ipiv,      \
                                          type *            work,      \
                                          blas_int_t        lwork,     \
                                          blas_int_t &      info ) {   \
        ADD_FLOPS( FLOPS_##flops( n ) );                               \
        info = 0;                                                      \
        func( & n, A, & ldA, ipiv, work, & lwork, & info ); }

HLIB_GETRI_FUNC( float,             sgetri_, SGETRI )
HLIB_GETRI_FUNC( double,            dgetri_, CGETRI )
HLIB_GETRI_FUNC( Complex< float >,  cgetri_, DGETRI )
HLIB_GETRI_FUNC( Complex< double >, zgetri_, ZGETRI )

#undef HLIB_GETRI_FUNC

//
// *geqrf
//
template <typename T>  
void geqrf ( const blas_int_t  n,
             const blas_int_t  m,
             T *               A,
             const blas_int_t  ldA,
             T *               tau,
             T *               work,
             blas_int_t        lwork,
             blas_int_t &      info );

#define HLIB_GEQRF_FUNC( type, func, flops )                            \
    template <> inline void geqrf<type> ( const blas_int_t  n,          \
                                          const blas_int_t  m,          \
                                          type *            A,          \
                                          const blas_int_t  ldA,        \
                                          type *            tau,        \
                                          type *            work,       \
                                          blas_int_t        lwork,      \
                                          blas_int_t &      info ) {    \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, tau, work, & lwork, & info ); }

HLIB_GEQRF_FUNC( float,             sgeqrf_, SGEQRF )
HLIB_GEQRF_FUNC( double,            dgeqrf_, CGEQRF )
HLIB_GEQRF_FUNC( Complex< float >,  cgeqrf_, DGEQRF )
HLIB_GEQRF_FUNC( Complex< double >, zgeqrf_, ZGEQRF )

#undef HLIB_GEQRF_FUNC

//
// *orgqr
//
template <typename T>  
void orgqr ( const blas_int_t  n,
             const blas_int_t  m,
             const blas_int_t  k,
             T *               A,
             const blas_int_t  ldA,
             T *               tau,
             T *               work,
             blas_int_t        lwork,
             blas_int_t &      info );

#define HLIB_ORGQR_FUNC( type, func, flops )                            \
    template <> inline void orgqr<type> ( const blas_int_t  n,          \
                                          const blas_int_t  m,          \
                                          const blas_int_t  k,          \
                                          type *            A,          \
                                          const blas_int_t  ldA,        \
                                          type *            tau,        \
                                          type *            work,       \
                                          blas_int_t        lwork,      \
                                          blas_int_t &      info ) {    \
        ADD_FLOPS( FLOPS_##flops( n, m, k ) );                          \
        info = 0;                                                       \
        func( & n, & m, & k, A, & ldA, tau, work, & lwork, & info ); }

HLIB_ORGQR_FUNC( float,             sorgqr_, SORGQR )
HLIB_ORGQR_FUNC( double,            dorgqr_, DORGQR )
HLIB_ORGQR_FUNC( Complex< float >,  cungqr_, CUNGQR )
HLIB_ORGQR_FUNC( Complex< double >, zungqr_, ZUNGQR )

#undef HLIB_ORGQR_FUNC

//
// *geqp3
//
template <typename T>  
void geqp3 ( const blas_int_t  n,
             const blas_int_t  m,
             T *               A,
             const blas_int_t  ldA,
             blas_int_t *      jpvt,
             T *               tau,
             T *               work,
             blas_int_t        lwork,
             typename real_type< T >::type_t *  rwork,
             blas_int_t &      info );

#define HLIB_GEQP3_FUNC( type, func )                                   \
    template <> inline void geqp3<type> ( const blas_int_t  n,          \
                                          const blas_int_t  m,          \
                                          type *            A,          \
                                          const blas_int_t  ldA,        \
                                          blas_int_t *      jpvt,       \
                                          type *            tau,        \
                                          type *            work,       \
                                          blas_int_t        lwork,      \
                                          typename real_type< type >::type_t *, \
                                          blas_int_t &      info ) {    \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, jpvt, tau, work, & lwork, & info ); }

HLIB_GEQP3_FUNC( float,  sgeqp3_ )
HLIB_GEQP3_FUNC( double, dgeqp3_ )

#undef HLIB_GEQP3_FUNC

#define HLIB_GEQP3_FUNC( type, func )                                   \
    template <> inline void geqp3<type> ( const blas_int_t  n,          \
                                          const blas_int_t  m,          \
                                          type *            A,          \
                                          const blas_int_t  ldA,        \
                                          blas_int_t *      jpvt,       \
                                          type *            tau,        \
                                          type *            work,       \
                                          blas_int_t        lwork,      \
                                          typename real_type< type >::type_t *  rwork, \
                                          blas_int_t &      info ) {    \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, jpvt, tau, work, & lwork, rwork, & info ); }

HLIB_GEQP3_FUNC( Complex< float >,  cgeqp3_ )
HLIB_GEQP3_FUNC( Complex< double >, zgeqp3_ )

#undef HLIB_GEQP3_FUNC

#if HAS_GEQP3_TRUNC == 1

//
// *geqp3trunc
//
template <typename T>  
void geqp3trunc ( const blas_int_t                  n,
                  const blas_int_t                  m,
                  T *                               A,
                  const blas_int_t                  ldA,
                  blas_int_t *                      jpvt,
                  T *                               tau,
                  blas_int_t &                      ntrunc,
                  typename real_type< T >::type_t   atrunc,
                  typename real_type< T >::type_t   rtrunc,
                  T *                               work,
                  blas_int_t                        lwork,
                  typename real_type< T >::type_t * rwork,
                  blas_int_t &                      info );

#define HLIB_GEQP3TRUNC_FUNC( type, func )                              \
    template <> inline void geqp3trunc<type> ( const blas_int_t  n,     \
                                               const blas_int_t  m,     \
                                               type *            A,     \
                                               const blas_int_t  ldA,   \
                                               blas_int_t *      jpvt,  \
                                               type *            tau,   \
                                               blas_int_t &      ntrunc, \
                                               typename real_type< type >::type_t  atrunc, \
                                               typename real_type< type >::type_t  rtrunc, \
                                               type *            work,  \
                                               blas_int_t        lwork, \
                                               typename real_type< type >::type_t *, \
                                               blas_int_t &      info ) { \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, jpvt, tau, & ntrunc, & atrunc, & rtrunc, work, & lwork, & info ); \
    }

HLIB_GEQP3TRUNC_FUNC( float,  sgeqp3trunc_ )
HLIB_GEQP3TRUNC_FUNC( double, dgeqp3trunc_ )

#undef HLIB_GEQP3TRUNC_FUNC

#define HLIB_GEQP3TRUNC_FUNC( type, func )                              \
    template <> inline void geqp3trunc<type> ( const blas_int_t  n,     \
                                               const blas_int_t  m,     \
                                               type *            A,     \
                                               const blas_int_t  ldA,   \
                                               blas_int_t *      jpvt,  \
                                               type *            tau,   \
                                               blas_int_t &      ntrunc, \
                                               typename real_type< type >::type_t  atrunc, \
                                               typename real_type< type >::type_t  rtrunc, \
                                               type *            work,  \
                                               blas_int_t        lwork, \
                                               typename real_type< type >::type_t *  rwork, \
                                               blas_int_t &      info ) { \
        info = 0;                                                       \
        func( & n, & m, A, & ldA, jpvt, tau, & ntrunc, & atrunc, & rtrunc, work, & lwork, rwork, & info ); }

HLIB_GEQP3TRUNC_FUNC( Complex< float >,  cgeqp3trunc_ )
HLIB_GEQP3TRUNC_FUNC( Complex< double >, zgeqp3trunc_ )

#undef HLIB_GEQP3TRUNC_FUNC

#endif

//
// *syev/*heev
//
template <typename T>  
void heev ( const char                        jobz,
            const char                        uplo,
            const blas_int_t                  n,
            T *                               A,
            const blas_int_t                  ldA,
            typename real_type< T >::type_t *  W,
            T *                               work,
            const blas_int_t                  lwork,
            typename real_type< T >::type_t *  rwork,
            blas_int_t &                      info );

#define HLIB_HEEV_FUNC( type, func )                                    \
    template <> inline void heev<type> ( const char        jobz,              \
                                         const char        uplo,        \
                                         const blas_int_t  n,           \
                                         type *            A,           \
                                         const blas_int_t  ldA,         \
                                         type *            W,           \
                                         type *            work,        \
                                         const blas_int_t  lwork,       \
                                         type *,                        \
                                         blas_int_t &      info ) {     \
        info = 0;                                                       \
        func( & jobz, & uplo, & n, A, & ldA, W, work, & lwork, & info ); }

HLIB_HEEV_FUNC( float,  ssyev_ )
HLIB_HEEV_FUNC( double, dsyev_ )

#undef HLIB_HEEV_FUNC

// #define HLIB_HEEV_FUNC( type, func )                                        
//     template <>  void heev<type> ( const char jobz, const char uplo, const blas_int_t n, 
//                                    type * A, const blas_int_t ldA, real_type< type >::type_t * W,   
//                                    type * work, const blas_int_t lwork, real_type< type >::type_t * rwork, blas_int_t & info ) { 
//         info = 0;                                                       
//         func( & jobz, & uplo, & n, A, & ldA, W, work, & lwork, rwork, & info ); }

// HLIB_HEEV_FUNC( Complex< float >,  cheev_ )
// HLIB_HEEV_FUNC( Complex< double >, zheev_ )

// #undef HLIB_HEEV_FUNC

//
// *syevx/*heevx
//
template <typename T>  
void heevx ( const char                        jobz,
             const char                        uplo,
             const blas_int_t                  n,
             T *                               A,
             const blas_int_t                  ldA,
             const blas_int_t                  il,
             const blas_int_t                  iu,
             blas_int_t &                      m,
             typename real_type< T >::type_t *  W,
             T *                               Z,
             const blas_int_t                  ldZ,
             T *                               work,
             const blas_int_t                  lwork,
             blas_int_t *                      iwork,
             blas_int_t &                      info );

#define HLIB_HEEVX_FUNC( type, func )                                   \
    template <> inline void heevx<type> ( const char                  jobz,   \
                                          const char                  uplo, \
                                          const blas_int_t            n, \
                                          type *                      A, \
                                          const blas_int_t            ldA, \
                                          const blas_int_t            il, \
                                          const blas_int_t            iu, \
                                          blas_int_t &                m, \
                                          real_type< type >::type_t *  W, \
                                          type *                      Z, \
                                          const blas_int_t            ldZ, \
                                          type *                      work, \
                                          const blas_int_t            lwork, \
                                          blas_int_t *                iwork, \
                                          blas_int_t &                info ) { \
        char                      range = 'I';                          \
        real_type< type >::type_t  vl, vu;                               \
        real_type< type >::type_t  abstol = 0;                           \
        blas_int_t                       ifail  = 0;                    \
        info = 0;                                                       \
        func( & jobz, & range, & uplo, & n, A, & ldA, & vl, & vu, & il, & iu, \
              & abstol, & m, W, Z, & ldZ, work, & lwork, iwork, & ifail, & info ); }

HLIB_HEEVX_FUNC( float,  ssyevx_ )
HLIB_HEEVX_FUNC( double, dsyevx_ )

#undef HLIB_HEEVX_FUNC

//
// sstev/dstev
//
template <typename T>  
void stev ( const char        jobz,
            const blas_int_t  n,
            T *               D,
            T *               E,
            T *               Z,
            const blas_int_t  ldZ,
            T *               work,
            blas_int_t &      info );

#define HLIB_STEV_FUNC( type, func )                            \
    template <> inline void stev<type> ( const char        jobz,      \
                                         const blas_int_t  n,         \
                                         type *            D,         \
                                         type *            E,         \
                                         type *            Z,         \
                                         const blas_int_t  ldZ,       \
                                         type *            work,      \
                                         blas_int_t &      info ) {   \
        info = 0;                                                     \
        func( & jobz, & n, D, E, Z, & ldZ, work, & info ); }

HLIB_STEV_FUNC( float,  sstev_ )
HLIB_STEV_FUNC( double, dstev_ )

#undef HLIB_STEV_FUNC

//
// *gesvd
//
template <typename T>  
void gesvd ( const char                        jobu,
             const char                        jobv,
             const blas_int_t                  n,
             const blas_int_t                  m,
             T *                               A,
             const blas_int_t                  ldA,
             typename real_type< T >::type_t *  S,
             T *                               U,
             const blas_int_t                  ldU,
             T *                               VT,
             const blas_int_t                  ldVT,
             T *                               work,
             const blas_int_t                  lwork,
             typename real_type< T >::type_t *  rwork,
             blas_int_t &                      info );

#define HLIB_GESVD_FUNC( type, func, flops )                            \
    template <> inline void gesvd<type> ( const char                  jobu, \
                                          const char                  jobv, \
                                          const blas_int_t            n, \
                                          const blas_int_t            m, \
                                          type *                      A, \
                                          const blas_int_t            ldA, \
                                          real_type< type >::type_t *  S, \
                                          type *                      U, \
                                          const blas_int_t            ldU, \
                                          type *                      VT, \
                                          const blas_int_t            ldVT, \
                                          type *                      work, \
                                          const blas_int_t            lwork, \
                                          real_type< type >::type_t *,  \
                                          blas_int_t &                info ) { \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & jobu, & jobv, & n, & m, A, & ldA, S, U, & ldU, VT, & ldVT, work, & lwork, & info ); }

HLIB_GESVD_FUNC( float,  sgesvd_, SGESVD )
HLIB_GESVD_FUNC( double, dgesvd_, DGESVD )

#undef  HLIB_GESVD_FUNC

#define HLIB_GESVD_FUNC( type, func, flops )                            \
    template <> inline void gesvd<type> ( const char                  jobu, \
                                          const char                  jobv, \
                                          const blas_int_t            n, \
                                          const blas_int_t            m, \
                                          type *                      A, \
                                          const blas_int_t            ldA, \
                                          real_type< type >::type_t *  S, \
                                          type *                      U, \
                                          const blas_int_t            ldU, \
                                          type *                      VT, \
                                          const blas_int_t            ldVT, \
                                          type *                      work, \
                                          const blas_int_t            lwork, \
                                          real_type< type >::type_t *  rwork, \
                                          blas_int_t &                info ) { \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & jobu, & jobv, & n, & m, A, & ldA, S, U, & ldU, VT, & ldVT, work, & lwork, rwork, & info ); }

HLIB_GESVD_FUNC( Complex< float >,  cgesvd_, CGESVD )
HLIB_GESVD_FUNC( Complex< double >, zgesvd_, CGESVD )

#undef  HLIB_GESVD_FUNC

//
// *gesdd
//
template <typename T>  
void gesdd ( const char                        jobz,
             const blas_int_t                  n,
             const blas_int_t                  m,
             T *                               A,
             const blas_int_t                  ldA,
             typename real_type< T >::type_t *  S,
             T *                               U,
             const blas_int_t                  ldU,
             T *                               VT,
             const blas_int_t                  ldVT,
             T *                               work,
             const blas_int_t                  lwork,
             typename real_type< T >::type_t *  rwork,
             blas_int_t *                      iwork,
             blas_int_t &                      info );

#define HLIB_GESDD_FUNC( type, func, flops )                            \
    template <> inline void gesdd<type> ( const char                  jobz, \
                                          const blas_int_t            n, \
                                          const blas_int_t            m, \
                                          type *                      A, \
                                          const blas_int_t            ldA, \
                                          real_type< type >::type_t *  S, \
                                          type *                      U, \
                                          const blas_int_t            ldU, \
                                          type *                      VT, \
                                          const blas_int_t            ldVT, \
                                          type *                      work, \
                                          const blas_int_t            lwork, \
                                          real_type< type >::type_t *,  \
                                          blas_int_t *                iwork, \
                                          blas_int_t &                info ) { \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & jobz, & n, & m, A, & ldA, S, U, & ldU, VT, & ldVT, work, & lwork, iwork, & info ); }

HLIB_GESDD_FUNC( float,  sgesdd_, SGESDD )
HLIB_GESDD_FUNC( double, dgesdd_, DGESDD )

#undef  HLIB_GESDD_FUNC

#define HLIB_GESDD_FUNC( type, func, flops )                            \
    template <> inline void gesdd<type> ( const char                  jobz, \
                                          const blas_int_t            n, \
                                          const blas_int_t            m, \
                                          type *                      A, \
                                          const blas_int_t            ldA, \
                                          real_type< type >::type_t *  S, \
                                          type *                      U, \
                                          const blas_int_t            ldU, \
                                          type *                      VT, \
                                          const blas_int_t            ldVT, \
                                          type *                      work, \
                                          const blas_int_t            lwork, \
                                          real_type< type >::type_t *  rwork, \
                                          blas_int_t *                iwork, \
                                          blas_int_t &                info ) { \
        ADD_FLOPS( FLOPS_##flops( n, m ) );                             \
        info = 0;                                                       \
        func( & jobz, & n, & m, A, & ldA, S, U, & ldU, VT, & ldVT, work, & lwork, rwork, iwork, & info ); }

HLIB_GESDD_FUNC( Complex< float >,  cgesdd_, CGESDD )
HLIB_GESDD_FUNC( Complex< double >, zgesdd_, ZGESDD )

#undef  HLIB_GESDD_FUNC

template <typename T>  
void gesvj ( const char                        joba,
             const char                        jobu,
             const char                        jobv,
             const blas_int_t                  n,
             const blas_int_t                  m,
             T *                               A,
             const blas_int_t                  ldA,
             typename real_type< T >::type_t *  S,
             const blas_int_t                  mv,
             T *                               V,
             const blas_int_t                  ldV,
             T *                               cwork,
             const blas_int_t                  lwork,
             typename real_type< T >::type_t * rwork,
             const blas_int_t                  lrwork,
             blas_int_t &                      info );

#define HLIB_GESVJ_FUNC( type, func )                                   \
    template <> inline void gesvj<type> ( const char                  joba, \
                                          const char                  jobu, \
                                          const char                  jobv, \
                                          const blas_int_t            n, \
                                          const blas_int_t            m, \
                                          type *                      A, \
                                          const blas_int_t            ldA, \
                                          real_type< type >::type_t * S, \
                                          const blas_int_t            mv, \
                                          type *                      V, \
                                          const blas_int_t            ldV, \
                                          type *                      cwork, \
                                          const blas_int_t            lwork, \
                                          real_type< type >::type_t *,  \
                                          const blas_int_t,             \
                                          blas_int_t &                info ) { \
        info = 0;                                                       \
        func( & joba, & jobu, & jobv, & n, & m, A, & ldA, S, & mv, V, & ldV, cwork, & lwork, & info ); }

HLIB_GESVJ_FUNC( float,  sgesvj_ )
HLIB_GESVJ_FUNC( double, dgesvj_ )

#undef  HLIB_GESVJ_FUNC

#define HLIB_GESVJ_FUNC( type, func )                                   \
    template <> inline void gesvj<type> ( const char                  joba, \
                                          const char                  jobu, \
                                          const char                  jobv, \
                                          const blas_int_t            n, \
                                          const blas_int_t            m, \
                                          type *                      A, \
                                          const blas_int_t            ldA, \
                                          real_type< type >::type_t * S, \
                                          const blas_int_t            mv, \
                                          type *                      V, \
                                          const blas_int_t            ldV, \
                                          type *                      cwork, \
                                          const blas_int_t            lwork, \
                                          real_type< type >::type_t * rwork, \
                                          const blas_int_t            lrwork, \
                                          blas_int_t &                info ) { \
        info = 0;                                                       \
        func( & joba, & jobu, & jobv, & n, & m, A, & ldA, S, & mv, V, & ldV, cwork, & lwork, rwork, & lrwork, & info ); }

HLIB_GESVJ_FUNC( Complex< float >,  cgesvj_ )
HLIB_GESVJ_FUNC( Complex< double >, zgesvj_ )

#undef  HLIB_GESVJ_FUNC

//
// *larfg
//
template <typename T>  
void
larfg ( const blas_int_t  n,
        T &               alpha,
        T *               x,
        const blas_int_t  incx,
        T &               tau );

#define HLIB_LARFG_FUNC( type, func )                                   \
    template <> inline void larfg<type> ( const blas_int_t  n,          \
                                          type &            alpha,      \
                                          type *            x,          \
                                          const blas_int_t  incx,       \
                                          type &            tau ) {     \
        func( & n, & alpha, x, & incx, & tau ); }

HLIB_LARFG_FUNC( float,             slarfg_ )
HLIB_LARFG_FUNC( double,            dlarfg_ )
HLIB_LARFG_FUNC( Complex< float >,  clarfg_ )
HLIB_LARFG_FUNC( Complex< double >, zlarfg_ )

#undef HLIB_LARFG_FUNC

//
// *larf
//
// apply householder reflection
template <typename T>  
void
larf  ( const char        side,
        const blas_int_t  n,
        const blas_int_t  m,
        const T *         V,
        const blas_int_t  incv,
        const T           tau,
        T *               C,
        const blas_int_t  ldc,
        T *               work );

#define HLIB_LARF_FUNC( type, func )                                    \
    template <> inline void larf<type> ( const char        side,        \
                                         const blas_int_t  n,           \
                                         const blas_int_t  m,           \
                                         const type *      V,           \
                                         const blas_int_t  incv,        \
                                         const type        tau,         \
                                         type *            C,           \
                                         const blas_int_t  ldc,         \
                                         type *            work ) {     \
        func( & side, & n, & m, V, & incv, & tau, C, & ldc, work ); }

HLIB_LARF_FUNC( float,             slarf_ )
HLIB_LARF_FUNC( double,            dlarf_ )
HLIB_LARF_FUNC( Complex< float >,  clarf_ )
HLIB_LARF_FUNC( Complex< double >, zlarf_ )

#undef HLIB_LARF_FUNC

//! @endcond

}// namespace BLAS

}// namespace HLIB

#endif // __HLIB_BLAS_LAPACK_HH
