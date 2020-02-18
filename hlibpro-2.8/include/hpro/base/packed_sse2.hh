#ifndef __HLIB_PACKED_SSE2_HH
#define __HLIB_PACKED_SSE2_HH
//
// Project     : HLib
// File        : packed_sse2.hh
// Description : datatype for packed (vector) operations using SSE2
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

/////////////////////////////////////////////////////////////
//
// SSE2 version of packed type
//

#if defined(__SSE2__)

#include <emmintrin.h>

#if USE_AMDLIBM == 1
extern "C" __m128d  amd_vrd2_exp     ( __m128d );
extern "C" void     amd_vrd2_sincos  ( __m128d, __m128d *, __m128d * );
#endif

#if USE_ACML == 1
extern "C" __m128d  __vrd2_exp       ( __m128d );
extern "C" void     __vrd2_sincos    ( __m128d, __m128d *, __m128d * );
#endif

#if USE_SVML == 1
extern "C"  __m128d  _mm_exp_pd      ( __m128d );
extern "C"  __m128d  _mm_sincos_pd   ( __m128d *, __m128d );
#endif

#if USE_LIBMVEC == 1
extern "C" __m128d  _ZGVbN2v_exp  ( __m128d  x );
extern "C" __m128d  _ZGVbN2v_sin  ( __m128d  x );
extern "C" __m128d  _ZGVbN2v_cos  ( __m128d  x );
#endif

namespace HLIB
{

//
// SSE2 functions
//
template <>
struct simd_traits< double, ISA_SSE2 >
{
    // SIMD base type
    using  value_t  = double;

    // SIMD vector type
    using  packed_t = __m128d;

    // SIMD instruction set
    enum { isa = ISA_SSE2 };

    // SIMD vector size
    enum { vector_size = 2 };

    //
    // SIMD functions
    //

    static packed_t  zero ()                   { return _mm_set1_pd( 0.0 ); }
    static packed_t  fill ( const value_t  f ) { return _mm_set1_pd( f ); }
    
    static packed_t  load  ( const value_t *  f ) { return _mm_loadu_pd( f ); }
    static void      store ( const packed_t   a,
                             value_t *        f ) { _mm_storeu_pd( f, a ); }
    
    static packed_t  add   ( const packed_t  x,
                             const packed_t  y ) { return _mm_add_pd( x, y ); }
    static packed_t  sub   ( const packed_t  x,
                             const packed_t  y ) { return _mm_sub_pd( x, y ); }
    static packed_t  mul   ( const packed_t  x,
                             const packed_t  y ) { return _mm_mul_pd( x, y ); }
    static packed_t  div   ( const packed_t  x,
                             const packed_t  y ) { return _mm_div_pd( x, y ); }

    static packed_t  muladd     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return add( mul( x, y ), z ); }
    static packed_t  negmuladd  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( z, mul( x, y ) ); }
    static packed_t  mulsub     ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( mul( x, y ), z ); }
    static packed_t  negmulsub  ( const packed_t  x,
                                  const packed_t  y,
                                  const packed_t  z ) { return sub( sub( zero(), mul( x, y ) ), z ); }

    static packed_t  sqrt   ( const packed_t  x ) { return _mm_sqrt_pd( x ); }
    
    static packed_t  rsqrt  ( const packed_t  x )
    {
        static const packed_t  vthree = fill( 3.0 );
        static const packed_t  vhalf  = fill( 0.5 );
        packed_t               res;
    
        // r = 1/√(x) in single precision
        res = _mm_cvtps_pd( _mm_rsqrt_ps( _mm_cvtpd_ps( x ) ) );

        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );
    
        // r = r · 1/2 (3.0 - r·r·x)
        res = mul( res,  mul( sub( vthree, mul( mul( res, res ), x  ) ), vhalf ) );

        return res;
    }
    
    static packed_t  exp    ( const packed_t  x )
    {
        #if USE_SVML == 1
    
        return _mm_exp_pd( x );

        #elif USE_LIBMVEC == 1
    
        return _ZGVbN2v_exp( x );
    
        #elif USE_AMDLIBM == 1
    
        return amd_vrd2_exp( x );

        #elif USE_ACML == 1
    
        return __vrd2_exp( x );

        #else
    
        //
        // fall back to standard floating point functions
        //
    
        double  sx[2];

        store( x, sx );

        sx[0] = std::exp( sx[0] );
        sx[1] = std::exp( sx[1] );

        return load( sx );
    
        #endif
    }

    static void  sincos ( const packed_t   a,
                          packed_t &       s,
                          packed_t &       c )
    {
        #if USE_SVML == 1

        s = _mm_sincos_pd( & c, a );
    
        #elif USE_LIBMVEC == 1

        s = _ZGVbN2v_sin( a );
        c = _ZGVbN2v_cos( a );
    
        #elif USE_AMDLIBM == 1
    
        amd_vrd2_sincos( a, & s, & c );
    
        #elif USE_ACML == 1
    
        __vrd2_sincos( a, & s, & c );
    
        #else

        //
        // fall back to standard floating point functions
        //
    
        double  sa[2], ss[2], sc[2];

        store( a, sa );

        #  if HAS_SINCOS == 1    

        ::sincos( sa[0], & ss[0], & sc[0] );
        ::sincos( sa[1], & ss[1], & sc[1] );

        #  else

        ss[0] = std::sin( sa[0] );
        ss[1] = std::sin( sa[1] );
        sc[0] = std::cos( sa[0] );
        sc[1] = std::cos( sa[1] );

        #  endif
    
        s = load( ss );
        c = load( sc );
    
        #endif
    }
};

template <>
struct packed< double, ISA_SSE2 >
{
    // SIMD base type
    using  value_t  = typename simd_traits< double, ISA_SSE2 >::value_t;

    // SIMD vector type
    using  packed_t = typename simd_traits< double, ISA_SSE2 >::packed_t;

    // SIMD vector size
    enum { isa = simd_traits< double, ISA_SSE2 >::isa };

    // SIMD vector size
    enum { vector_size = simd_traits< double, ISA_SSE2 >::vector_size };

    // SIMD data
    packed_t   x;

    // ctors
    packed ()                   : x( simd_traits< double, ISA_SSE2 >::zero() )   {}
    packed ( packed_t       y ) : x( y )          {}
    packed ( const value_t  f ) : x( simd_traits< double, ISA_SSE2 >::fill( f ) ) {}
    packed ( const value_t  a,
             const value_t  b ) : x( _mm_setr_pd( a, b ) ) {}
};

//
// yields { a[0], b[0] }
//
inline
packed< double, ISA_SSE2 >
unpacklo ( const packed< double, ISA_SSE2 >  a,
           const packed< double, ISA_SSE2 >  b )
{
    return _mm_unpacklo_pd( a.x, b.x );
}

//
// yields { a[1], b[1] }
//
inline
packed< double, ISA_SSE2 >
unpackhi ( const packed< double, ISA_SSE2 >  a,
           const packed< double, ISA_SSE2 >  b )
{
    return _mm_unpackhi_pd( a.x, b.x );
}

//
// a·b with complex a,b
//
inline
packed< double, ISA_SSE2 >
zmul ( const packed< double, ISA_SSE2 >  a,
       const packed< double, ISA_SSE2 >  b )
{
    const packed< double, ISA_SSE2 >  t1( _mm_mul_pd( a.x, b.x ) );
    const packed< double, ISA_SSE2 >  t2( _mm_mul_pd( a.x, _mm_shuffle_pd( b.x, b.x, 1 ) ) );
     
    return _mm_unpacklo_pd( _mm_sub_pd( t1.x, _mm_shuffle_pd( t1.x, t1.x, 1 ) ),
                            _mm_add_pd( t2.x, _mm_shuffle_pd( t2.x, t2.x, 1 ) ) );
}

//
// a/b with complex a,b
//
inline
packed< double, ISA_SSE2 >
zdiv ( const packed< double, ISA_SSE2 >  a,
       const packed< double, ISA_SSE2 >  b )
{
    const packed< double, ISA_SSE2 >  t1( _mm_mul_pd( a.x, b.x ) );
    const packed< double, ISA_SSE2 >  t2( _mm_mul_pd( _mm_shuffle_pd( a.x, a.x, 1 ), b.x ) );
    const packed< double, ISA_SSE2 >  t3( _mm_mul_pd( b.x, b.x ) );
        
    return _mm_div_pd( _mm_unpacklo_pd( _mm_add_pd( t1.x, _mm_shuffle_pd( t1.x, t1.x, 1 ) ),
                                        _mm_sub_pd( t2.x, _mm_shuffle_pd( t2.x, t2.x, 1 ) ) ),
                       _mm_add_pd( t3.x, _mm_shuffle_pd( t3.x, t3.x, 1 ) ) );
}

}// namespace HLIB

#endif // __SSE2__

#endif // __HLIB_PACKED_SSE2_HH

