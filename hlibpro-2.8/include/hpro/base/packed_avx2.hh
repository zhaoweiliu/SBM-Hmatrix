#ifndef __HLIB_PACKED_AVX2_HH
#define __HLIB_PACKED_AVX2_HH
//
// Project     : HLib
// File        : packed_avx2.hh
// Description : datatype for packed (vector) operations using AVX2
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

/////////////////////////////////////////////////////////////
//
// AVX2 version of packed type
//

#if defined(__AVX2__)

#include <immintrin.h>

#if USE_AMDLIBM == 1
extern "C" void  amd_vrda_exp     ( int  n, double *  f, double *  res );
extern "C" void  amd_vrda_sincos  ( int  n, double *  f, double *  res_s, double *  res_c );
#endif

#if USE_ACML == 1
extern "C" void  vrda_exp         ( int  n, double *  f, double *  res );
extern "C" void  vrda_sincos      ( int  n, double *  f, double *  res_s, double *  res_c );
#endif

#if USE_SVML == 1
extern "C"  __m256d  _mm256_exp_pd     ( __m256d );
extern "C"  __m256d  _mm256_sincos_pd  ( __m256d *, __m256d );
#endif

#if USE_LIBMVEC == 1
extern "C" __m256d  _ZGVdN4v_exp  ( __m256d  x );
extern "C" __m256d  _ZGVdN4v_sin  ( __m256d  x );
extern "C" __m256d  _ZGVdN4v_cos  ( __m256d  x );
#endif

namespace HLIB
{

template <>
struct simd_traits< double, ISA_AVX2 >
{
    // SIMD base type
    using  value_t  = double;

    // SIMD vector type
    using  packed_t = __m256d;

    // SIMD instruction set
    enum { isa = ISA_AVX2 };

    // SIMD vector size
    enum { vector_size = 4 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm256_set1_pd( 0.0 ); }
    static packed_t  fill ( const value_t  f ) { return _mm256_set1_pd( f ); }
    
    static packed_t  load  ( const value_t *  f ) { return _mm256_loadu_pd( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm256_storeu_pd( f, a ); }
    
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_add_pd( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_sub_pd( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_mul_pd( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return _mm256_div_pd( x, y ); }

    static packed_t  muladd     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm256_fmadd_pd( x, y, z ); }
    static packed_t  negmuladd  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm256_fnmadd_pd( x, y, z ); }
    static packed_t  mulsub     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm256_fmsub_pd( x, y, z ); }
    static packed_t  negmulsub  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return _mm256_fnmsub_pd( x, y, z ); }

    static packed_t  sqrt   ( const packed_t  x ) { return _mm256_sqrt_pd( x ); }
    
    static packed_t  rsqrt  ( const packed_t  x )
    {
        static const packed_t  vthree = fill( 3.0 );
        static const packed_t  vhalf  = fill( 0.5 );
        packed_t               res;

        // r = 1/√(x) in single precision
        res = _mm256_cvtps_pd( _mm_rsqrt_ps( _mm256_cvtpd_ps( x ) ) );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res, mul( negmuladd( mul( res, res ), x, vthree ), vhalf ) );
    
        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res, mul( negmuladd( mul( res, res ), x, vthree ), vhalf ) );

        return res;
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if USE_SVML == 1
    
        return _mm256_exp_pd( x );
    
        #elif USE_LIBMVEC == 1
    
        return _ZGVdN4v_exp( x );
    
        #elif USE_AMDLIBM == 1

        packed_t  res;
    
        amd_vrda_exp( 4,
                      reinterpret_cast< double * >( const_cast< packed_t * >( & x ) ),
                      reinterpret_cast< double * >( & res ) );

        return res;
    
        #elif USE_ACML == 1

        packed_t  res;
    
        vrda_exp( 4,
                  reinterpret_cast< double * >( const_cast< packed_t * >( & x ) ),
                  reinterpret_cast< double * >( & res ) );

        return res;
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        double  sx[4];

        store( x, sx );

        sx[0] = std::exp( sx[0] );
        sx[1] = std::exp( sx[1] );
        sx[2] = std::exp( sx[2] );
        sx[3] = std::exp( sx[3] );

        return load( sx );
    
        #endif
    }

    static void      sincos ( const packed_t   a,
                              packed_t &       s,
                              packed_t &       c )
    {
        #if USE_SVML == 1

        s = _mm256_sincos_pd( & c, a );
    
        #elif USE_LIBMVEC == 1

        s = _ZGVdN4v_sin( a );
        c = _ZGVdN4v_cos( a );
    
        #elif USE_AMDLIBM == 1

        amd_vrda_sincos( 4,
                         reinterpret_cast< double * >( const_cast< packed_t * >( & a ) ),
                         reinterpret_cast< double * >( & s ),
                         reinterpret_cast< double * >( & c ) );
        
        #elif USE_ACML == 1

        vrda_sincos( 4,
                     reinterpret_cast< double * >( const_cast< packed_t * >( & a ) ),
                     reinterpret_cast< double * >( & s ),
                     reinterpret_cast< double * >( & c ) );
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        double  sa[4], ss[4], sc[4];

        store( a, sa );

        #  if HAS_SINCOS == 1    

        ::sincos( sa[0], & ss[0], & sc[0] );
        ::sincos( sa[1], & ss[1], & sc[1] );
        ::sincos( sa[2], & ss[2], & sc[2] );
        ::sincos( sa[3], & ss[3], & sc[3] );

        #  else

        ss[0] = std::sin( sa[0] );
        ss[1] = std::sin( sa[1] );
        ss[2] = std::sin( sa[2] );
        ss[3] = std::sin( sa[3] );

        sc[0] = std::cos( sa[0] );
        sc[1] = std::cos( sa[1] );
        sc[2] = std::cos( sa[2] );
        sc[3] = std::cos( sa[3] );

        #  endif
    
        s = load( ss );
        c = load( sc );

        #endif
    }
};

template <>
struct packed< double, ISA_AVX2 >
{
    // SIMD base type
    using  value_t  = typename simd_traits< double, ISA_AVX2 >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< double, ISA_AVX2 >::packed_t;

    // SIMD vector size
    enum { isa = simd_traits< double, ISA_AVX2 >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< double, ISA_AVX2 >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                    : x( simd_traits< double, ISA_AVX2 >::zero() )    {}
    packed ( const packed_t  y ) : x( y )                                          {}
    packed ( const value_t   f ) : x( simd_traits< double, ISA_AVX2 >::fill( f ) ) {}
    packed ( const value_t   a,
             const value_t   b,
             const value_t   c,
             const value_t   d ) : x( _mm256_setr_pd( a, b, c, d ) ) {}
};

}// namespace HLIB

#endif // __AVX2__

#endif // __HLIB_PACKED_AVX2_HH
