#ifndef __HLIB_COMPLEX_HH
#define __HLIB_COMPLEX_HH
//
// Project     : HLib
// File        : complex.hh
// Description : class for a complex type
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <cmath>
#include <ostream>
#include <string>

#if defined(__SSE3__)
#  include "hpro/base/packed.hh"
#endif

#include "hpro/base/basetypes.hh"

namespace HLIB
{

//
// use fast versions of exp, log and sincos
//

#if USE_AMDLIBM == 1

extern "C"
{
void  amd_sincos   ( double, double *, double * );
void  amd_sincosf  ( float,  float *,  float *  );
}

#endif

#if USE_ACML == 1

extern "C"
{
void  fastsincos   ( double, double *, double * );
void  fastsincosf  ( float,  float *,  float *  );
}

#endif

//!
//! \class  complex
//! \brief  Class for a complex numerical type.
//!
//!         - Initially implemented due to incompatibilities for std::complex
//!           between various compilers.
//!         - Now used since std::complex is much slower.
//!
template <typename T>
class Complex
{
public:
    //
    // public types
    //

    using  real_t = T;

protected:
    //
    // real and imaginary part
    //

    real_t  _re, _im;

public:
    ////////////////////////////////////////////////////
    //
    // constructors
    //

    //! construct complex number 0
    Complex () noexcept
    {
        _re = real_t(0);
        _im = real_t(0);
    }
    
    //! construct complex number \a ar + i·\a ai
    Complex ( const real_t  ar,
              const real_t  ai = 0.0 ) noexcept
    {
        _re = ar;
        _im = ai;
    }

    //! copy ctor
    template <typename R>
    Complex ( const Complex<R> &  z ) noexcept
    {
        _re = real_t( z.re() );
        _im = real_t( z.im() );
    }

    ////////////////////////////////////////////////////
    //
    // access data fields
    //

    //! return real part
    real_t    re  () const noexcept { return _re; }

    //! return imaginary part
    real_t    im  () const noexcept { return _im; }

    //! set real part to \a f
    void      re  ( const T  f ) noexcept { _re = f; }

    //! set imaginary part to \a f
    void      im  ( const T  f ) noexcept { _im = f; }
    
    ////////////////////////////////////////////////////
    //
    // standard arithmetic operations
    //

    //! return *this + \a f
    Complex<T> &  operator += ( const T  f ) noexcept { _re += f; return *this; }

    //! return *this - \a f
    Complex<T> &  operator -= ( const T  f ) noexcept { _re -= f; return *this; }

    //! return *this * \a f
    Complex<T> &  operator *= ( const T  f ) noexcept { _re *= f; _im *= f; return *this; }

    //! return *this / \a f
    Complex<T> &  operator /= ( const T  f ) noexcept { _re /= f; _im /= f; return *this; }


    //! return *this + \a z
    template <typename R>
    Complex<T> &  operator += ( const Complex<R> &  z ) noexcept
    {
        _re += real_t(z.re());
        _im += real_t(z.im());

        return *this;
    }

    //! return *this - \a z
    template <typename R>
    Complex<T> &  operator -= ( const Complex<R> &  z ) noexcept
    {
        _re -= z.re();
        _im -= z.im();

        return *this;
    }

    //! return *this * \a z
    template <typename R>
    Complex<T> &  operator *= ( const Complex<R> &  z ) noexcept
    {
        const real_t  r = _re * z.re() - _im * z.im();
        
        _im = _re * z.im() + _im * z.re();
        _re = r;
        
        return *this;
    }


    //! return *this / \a z
    template <typename R>
    Complex<T> &  operator /= ( const Complex<R> &  z ) noexcept
    {
        const real_t  r =  _re * z.re() + _im * z.im();
        const real_t  n = z.norm();
        
        _im = (_im * z.re() - _re * z.im()) / n;
        _re = r / n;

        return *this;
    }

    ////////////////////////////////////////////////////
    //
    // non-standard arithmetic operations
    //

    //!
    //! return conjugated number
    //!
    const Complex<T>  conj () const noexcept { return Complex<T>( _re, -_im ); }

    //!
    //! conjugate *this
    //!
    Complex<T> &  conjugate () noexcept { _im = -_im; return *this; }

    //!
    //! returns the magnitude
    //!
    real_t  abs () const noexcept
    {
        const real_t  s = std::max( std::abs( _re ), std::abs( _im ) );
        
        if ( s == real_t(0) )
            return real_t(0);
        
        const real_t  x = _re / s; 
        const real_t  y = _im / s;
        
        return s * std::sqrt( x * x + y * y );
    }
    
    //!
    //! returns the phase angle 
    //!
    real_t  arg () const noexcept
    {
        return std::atan2( _im, _re );
    }
    
    //!
    //! return squared magnitude
    //!
    real_t  norm () const noexcept { return _re*_re + _im*_im; }

    //!
    //! return square root
    //!
    const Complex<T>  sqrt () const noexcept
    {
        const real_t x = re();
        const real_t y = im();

        if ( x == real_t(0) )
        {
            const real_t t = std::sqrt( std::abs( y ) / real_t(2) );
            
            return Complex( t, ( y < real_t(0) ? -t : t ));
        }// if
        else
        {
            const real_t t = std::sqrt( real_t(2) * ( abs() + std::abs(x) ) );
            const real_t u = t / real_t(2);
            
            return ( x > real_t(0)       ?
                     Complex( u, y / t ) :
                     Complex( std::abs(y) / t,
                              ( y < real_t(0) ? -u : u ) ) );
        }// else
    }
    
    //!
    //! returns sign : x/|x| with x = *this
    //!
    const Complex<T> sign () const noexcept
    {
        Complex<T> t( *this );

        t /= t.abs();
        return t;
    }
    
    ////////////////////////////////////////////////////
    //
    // transcendentals
    //

    //!
    //! return cosine of *this
    //!
    const Complex<T> cos () const noexcept
    {
        return Complex< T >(   std::cos( _re ) * std::cosh( _im ),
                             - std::sin( _re ) * std::sinh( _im ) );
    }
    
    //!
    //! return hyperbolic cosine of *this
    //!
    const Complex<T> cosh () const noexcept
    {
        return Complex< T >( std::cosh( _re ) * std::cos( _im ),
                             std::sinh( _re ) * std::sin( _im ) );
    }
    
    //!
    //! return sine of *this
    //!
    const Complex<T> sin () const noexcept
    {
        return Complex< T >( std::sin( _re ) * std::cosh( _im ),
                             std::cos( _re ) * std::sinh( _im ) );
    }
    
    //!
    //! return hyperbolic sine  of *this
    //!
    const Complex<T> sinh () const noexcept
    {
        return Complex< T >( std::sinh( _re ) * std::cos( _im ),
                             std::cosh( _re ) * std::sin( _im ) );
    }
    
    //!
    //! return tangent (sin/cos) of *this
    //!
    const Complex<T> tan () const noexcept
    {
        const T  s_re  = std::sin( _re );
        const T  c_re  = std::cos( _re );
        const T  sh_im = std::sinh( _im );
        const T  ch_im = std::cosh( _im );
        
        return ( Complex< T >( s_re * ch_im,   c_re * sh_im ) /
                 Complex< T >( c_re * ch_im, - s_re * sh_im ) );
    }
    
    //!
    //! return hyperbolic tangent (sinh/cosh) of *this
    //!
    const Complex<T> tanh () const noexcept
    {
        const T  sh_re = std::sinh( _re );
        const T  ch_re = std::cosh( _re );
        const T  s_im  = std::sin( _im );
        const T  c_im  = std::cos( _im );
        
        return ( Complex< T >( sh_re * c_im, ch_re * s_im ) /
                 Complex< T >( ch_re * c_im, sh_re * s_im ) );
    }
    
    ////////////////////////////////////////////////////
    //
    // misc.
    //

    //
    // comparison
    //

    //! return true if *this == \a f 
    bool operator == ( const real_t  f ) const noexcept { return (_re == f) && (_im == real_t(0)); }

    //! return true if *this != \a f 
    bool operator != ( const real_t  f ) const noexcept { return (_re != f) || (_im != real_t(0)); }

    //! return true if *this == \a z 
    template <typename R>
    bool operator == ( const Complex<R> & z ) const noexcept { return (_re == z.re()) && (_im == z.im()); }

    //! return true if *this != \a z 
    template <typename R>
    bool operator != ( const Complex<R> & z ) const noexcept { return (_re != z.re()) || (_im != z.im()); }

    //
    // assignment
    //

    //! set *this = f
    Complex<T> & operator = ( const T  f ) noexcept
    {
        _re = f;
        _im = real_t(0);
        
        return *this;
    }

    //! set *this = z
    template <typename R>
    Complex & operator = ( const Complex<R> &  z ) noexcept
    {
        _re = real_t( z.re() );
        _im = real_t( z.im() );
        
        return *this;
    }
    
    //
    // output
    //

    //! return string represantation of *this
    std::string  to_string () const;
};


////////////////////////////////////////////////////
//
// stream I/O
//

//! write complex to stream
template <typename T>
std::ostream &
operator << ( std::ostream &        os,
              const Complex< T > &  c );

//! read complex from stream
template <typename T>
std::istream &
operator >> ( std::istream &        is,
              Complex< T > &        c );

////////////////////////////////////////////////////
//
// specialised internal functions
//

//
//! return *this * \a z
//
#if defined(__SSE3__)
template <>
template <>
inline
Complex< double > &
Complex< double >::operator *= ( const Complex< double > &  z ) noexcept
{
    packed< double, ISA_SSE3 >  t( zmul( packed< double, ISA_SSE3 >( _re, _im ),
                                         packed< double, ISA_SSE3 >( z._re, z._im ) ) );
    
    store( t, reinterpret_cast< double * >( this ) );
        
    return *this;
}
#endif

//
//! return *this / \a z
//
#if defined(__SSE3__)
template <>
template <>
inline
Complex< double > &
Complex< double >::operator /= ( const Complex< double > &  z ) noexcept
{
    packed< double, ISA_SSE3 >  t( zdiv( packed< double, ISA_SSE3 >( _re, _im ),
                                         packed< double, ISA_SSE3 >( z._re, z._im ) ) );

    store( t, reinterpret_cast< double * >( this ) );
    
    return *this;
}
#endif   

////////////////////////////////////////////////////
//
// access data fields
//

//! return real part of \a z
template <typename T> inline T  re  ( const Complex<T> &  z ) noexcept { return z.re(); }

//! return imaginary part of \a z
template <typename T> inline T  im  ( const Complex<T> &  z ) noexcept { return z.im(); }

// special versions for float/double to allow real/complex template functions
inline float   re  ( const float   f ) noexcept { return f; }
inline double  re  ( const double  f ) noexcept { return f; }
inline float   im  ( const float     ) noexcept { return 0.f; }
inline double  im  ( const double    ) noexcept { return 0.0; }

////////////////////////////////////////////////////
//
// standard arithmetic operations
//

//! return -\a z
template <typename T> 
inline const Complex<T> operator - ( const Complex<T> &  z ) noexcept
{ return Complex<T>( -z.re(), -z.im() ); }

//! return \a z + \a f
template <typename T> 
inline const Complex<T> operator + ( const Complex<T> &  z,
                                     const T             f ) noexcept
{ Complex<T> t( z ); return t += f; }

//! return \a z - \a f
template <typename T> 
inline const Complex<T> operator - ( const Complex<T> &  z,
                                     const T             f ) noexcept
{ Complex<T> t( z ); return t -= f; }

//! return \a z * \a f
template <typename T> 
inline const Complex<T> operator * ( const Complex<T> &  z,
                                     const T             f ) noexcept
{ Complex<T> t( z ); return t *= f; }

//! return \a z / \a f
template <typename T> 
inline const Complex<T> operator / ( const Complex<T> &  z,
                                     const T             f ) noexcept
{ Complex<T> t( z ); return t /= f; }

//! return \a f + \a z
template <typename T> 
inline const Complex<T> operator + ( const T             f,
                                     const Complex<T> &  z ) noexcept
{ Complex<T> t( z ); return t += f; }

//! return \a f - \a z
template <typename T> 
inline const Complex<T> operator - ( const T             f,
                                     const Complex<T> &  z ) noexcept
{ Complex<T> t( f ); return t -= z; }

//! return \a f * \a z
template <typename T> 
inline const Complex<T> operator * ( const T             f,
                                     const Complex<T> &  z ) noexcept
{ Complex<T> t( z ); return t *= f; }

//! return \a f / \a z
template <typename T> 
inline const Complex<T> operator / ( const T             f,
                                     const Complex<T> &  z ) noexcept
{ Complex<T> t( f ); return t /= z; }

//! return \a z1 + \a z2
template <typename T> 
inline const Complex<T> operator + ( const Complex<T> &  z1,
                                     const Complex<T> &  z2 ) noexcept
{ Complex<T> t( z1 ); return t += z2; }

//! return \a z1 - \a z2
template <typename T> 
inline const Complex<T> operator - ( const Complex<T> &  z1,
                                     const Complex<T> &  z2 ) noexcept
{ Complex<T> t( z1 ); return t -= z2; }

//! return \a z1 * \a z2
template <typename T> 
inline const Complex<T> operator * ( const Complex<T> &  z1,
                                     const Complex<T> &  z2 ) noexcept
{ Complex<T> t( z1 ); return t *= z2; }

//! return \a z1 / \a z2
template <typename T> 
inline const Complex<T> operator / ( const Complex<T> &  z1,
                                     const Complex<T> &  z2 ) noexcept
{ Complex<T> t( z1 ); return t /= z2; }

////////////////////////////////////////////////////
//
// polar to arithmetic conversion
//

//! return complex number to polar coordinates \a rho and \a theta
template <typename T> 
inline
const Complex<T>
polar ( const T  rho,
        const T  theta ) noexcept
{
    T  s, c;

    s = std::sin( theta );
    c = std::cos( theta );

    return Complex<T>( rho * c, rho * s );
}

template <> 
inline
const Complex< float >
polar ( const float  rho,
        const float  theta ) noexcept
{
    float  s, c;

#if USE_SVML == 1
    
    sincosf( theta, & s, & c );
        
#elif USE_AMDLIBM == 1
    
    amd_sincosf( theta, & s, & c );
        
#elif USE_ACML == 1
    
    fastsincosf( theta, & s, & c );
        
#elif HAS_SINCOS == 1
        
    ::sincosf( theta, & s, & c );
        
#else
    
    s = std::sin( theta );
    c = std::cos( theta );
    
#endif

    return Complex< float >( rho * c, rho * s );
}

template <> 
inline
const Complex< double >
polar ( const double  rho,
        const double  theta ) noexcept
{
    double  s, c;

#if defined(__MIC__)
    
    // REMARK: sincos on Intel MIC is _way_ slower than calling sin/cos individually
    s = std::sin( theta );
    c = std::cos( theta );
    
#elif USE_SVML == 1
    
    ::sincos( theta, & s, & c );
        
#elif USE_AMDLIBM == 1
    
    amd_sincos( theta, & s, & c );
        
#elif USE_ACML == 1
    
    fastsincos( theta, & s, & c );
        
#elif HAS_SINCOS == 1
        
    ::sincos( theta, & s, & c );
        
#else
    
    s = std::sin( theta );
    c = std::cos( theta );
    
#endif

    return Complex< double >( rho * c, rho * s );
}

////////////////////////////////////////////////////
//
// non-standard arithmetic operations
//

inline float  conj ( const float  f ) noexcept { return f; }
inline double conj ( const double f ) noexcept { return f; }

//! return conjugate value of \a z
template <typename T> 
inline const Complex<T>  conj ( const Complex<T> & z ) noexcept { return z.conj(); }

//! return absolute value of \a z
template <typename T> 
inline T                 abs ( const Complex<T> & z ) noexcept { return z.abs(); } 

//! return phase angle of \a z
template <typename T> 
inline T                 arg ( const Complex<T> & z ) noexcept { return z.arg(); } 

//! return squared magnitude of \a z
template <typename T> 
inline T                 norm ( const Complex<T> & z ) noexcept { return z.norm(); }

//! return square root of \a z
template <typename T> 
inline const Complex<T>  sqrt ( const Complex<T> & z ) noexcept { return z.sqrt(); }

//! return sign of \a z
template <typename T> 
inline const Complex<T>  sign ( const Complex<T> & z ) noexcept { return z.sign(); }

//! return exp( \a z )
template <typename T> 
inline const Complex<T>  exp  ( const Complex<T> & z ) noexcept
{
    if ( z.re() == T(0) ) return polar<T>( T(1), z.im() );
    else                  return polar<T>( std::exp( z.re() ), z.im() );
}

//! return log( \a z )
template <typename T> 
inline const Complex<T>  log ( const Complex<T> & z ) noexcept
{
    return Complex<T>( std::log( abs(z) ), arg(z) );
}

//! return log( \a z ) w.r.t. base 10
template <typename T> 
inline const Complex<T>  log10 ( const Complex<T> & z ) noexcept
{
    const T  ONE_OVER_LN10 = T(0.43429448190325182765112891891660);

    return log( z ) * T( ONE_OVER_LN10 );
}

//! return \a x to the power of \a y
template <typename T> 
inline const Complex<T>  pow  ( const Complex< T > &  x,
                                const Complex< T > &  y ) noexcept
{
    return ( x == T(0) ? T(0) : exp( y * log( x ) ) );
}

//! return \a x to the power of \a y
template <typename T> 
inline const Complex<T>  pow ( const Complex< T > &  x,
                               const T               y ) noexcept
{
    return ( x == T(0) ? T(0) : exp( y * log( x ) ) );
}

//! return \a x to the power of \a y
template <typename T> 
inline const Complex<T>  pow ( const Complex< T > &  x,
                               const int             y ) noexcept
{
    Complex< T >  res( T(1), T(0) );

    if ( y >= 0 )
    {
        for ( int  i = 0; i < y; ++i )
            res *= x;
    }// if
    else
    {
        const int  abs_y = std::abs( y );
        
        for ( int  i = 0; i < abs_y; ++i )
            res *= x;

        res = T(1) / res;
    }// else

    return res;
}

//! return sine of \a x
template <typename T> 
inline const Complex<T>  sin  ( const Complex< T > &  x ) noexcept
{
    return x.sin();
}

//! return hyperbolic sine of \a x
template <typename T> 
inline const Complex<T>  sinh ( const Complex< T > &  x ) noexcept
{
    return x.sinh();
}

//! return cosine of \a x
template <typename T> 
inline const Complex<T>  cos  ( const Complex< T > &  x ) noexcept
{
    return x.cos();
}

//! return hyperbolic cosine of \a x
template <typename T> 
inline const Complex<T>  cosh ( const Complex< T > &  x ) noexcept
{
    return x.cosh();
}

//! return tangent of \a x
template <typename T> 
inline const Complex<T>  tan  ( const Complex< T > &  x ) noexcept
{
    return x.tan();
}

//! return hyperbolic tangent of \a x
template <typename T> 
inline const Complex<T>  tanh ( const Complex< T > &  x ) noexcept
{
    return x.tanh();
}

////////////////////////////////////////////////////
//
// default complex type based on "real"
//

using  complex = Complex<real>;

////////////////////////////////////////////////////
//
// imaginary number
//


}// namespace

#endif
