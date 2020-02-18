#ifndef __HLIB_CONFIG_HH
#define __HLIB_CONFIG_HH
/*
 * Project     : HLib
 * File        : hlib-config.h
 * Description : configuration settings
 * Author      : Ronald Kriemann
 * Copyright   : Max-Planck-Institute MIS 2004-2018. All Rights Reserved.
 */

/* HLIBpro namespace override */
#define  HLIB                    HLIB

/* prefix all HLIBpro types/functions (except constants) */
#define  HLIBPF( name )          hlib_##name

/*
 * set to 1 if the following libraries are available and
 * set to 0 otherwise 
 */

#define USE_TBB                  1
#define USE_ZLIB                 1
#define USE_METIS                0
#define USE_SCOTCH               0
#define USE_CHACO                0
#define USE_FFTW3                1
#define USE_GSL                  1
#define USE_CAIRO                0
#define USE_HDF5                 0
#define USE_AMDLIBM              0
#define USE_ACML                 0
#define USE_SVML                 0
#define USE_LIBMVEC              0
#define USE_MKL                  0
#define USE_MKL_SEQ              0
#define USE_VAMPIRTRACE          0
#define USE_LIC_CHECK            0

/*
 * set to 1 if the following functions are available and
 * set to 0 otherwise 
 */

#define HAS_CLOCKGETTIME         1
#define HAS_GETTIMEOFDAY         1
#define HAS_GETTICKCOUNT         0
#define HAS_GETPROCESSTIMES      0
#define HAS_LOCALTIME_R          1
#define HAS_GETRUSAGE            1
#define HAS_LWPINFO              0

#define HAS_GETPAGESIZE          0 /* 1 */
#define HAS_MMAP                 0 /* yes */

#define HAS_STRERROR_R           1
#define HAS_BACKTRACE            1
#define HAS_CXXDEMANGLE          1

#define HAS_UNORDERED_MAP        0
#define HAS_BOOST_UNORDERED_MAP  1
#define HAS_BOOST_IOSTREAMS      1

#define HAS_GEJSV                0
#define HAS_GESVJ                0
#define HAS_SINCOS               0
#define HAS_GEQP3_TRUNC          0

#define HAS_ILP64                0

#define HAS_CPUID                1
#define HAS_SSE3                 1
#define HAS_AVX                  1
#define HAS_AVX2                 1
#define HAS_MIC                  0
#define HAS_AVX512F              1
#define HAS_VSX                  0

/*
 * set to 1 to activate various features or change
 * behaviour of algorithms
 */

/* activate to enable additional debugging features */
#define HLIB_DEBUG               0

/* activate argument tests in inline functions of BLAS module
 * (also activated by HLIB_DEBUG) */
#define HLIB_BLAS_TESTS          0

/* activate to enable single precision arithmetics */
#define HLIB_SINGLE_PREC         0

/* activate to enable flop counting */
#define HLIB_COUNT_FLOPS         0

/*
 * define type of network environment
 *   1 : sequential
 *   2 : MPI library
 *   3 : shared memory
 */

#define HLIB_NET_TYPE            1

#endif
