#ifndef __HLIB_SYSTEM_HH
#define __HLIB_SYSTEM_HH
//
// Project     : HLib
// File        : System.hh
// Description : module containing basic system routines
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include <cmath>
#include <string>
#include <limits>

#include "hpro/base/types.hh"
#include "hpro/base/String.hh"

namespace HLIB
{
    
///////////////////////////////////////////////////
//
// mathematical functions
//
///////////////////////////////////////////////////

#if USE_AMDLIBM == 1

extern "C" void  amd_sincos  ( double, double *, double * );

#endif

#if USE_ACML == 1

extern "C" void  fastsincos  ( double, double *, double * );

#endif

namespace Math
{
    //
    // constants
    //

    template < typename T >  constexpr T  pi () { return T(3.141592653589793238462643383279502884L); }

    //
    // minimum, maximum, intervals
    //

    // return minimum/maximum of \a f1 and \a f2
    template < typename T >  T  min     ( const T  f1, const T  f2 ) noexcept { return std::min( f1, f2 ); }
    template < typename T >  T  max     ( const T  f1, const T  f2 ) noexcept { return std::max( f1, f2 ); }

    // limit \a f to interval [f_min,f_max]
    template < typename T >  T  limit   ( const T  f_min,
                                          const T  f,
                                          const T  f_max ) noexcept { return std::min( f_max, std::max( f_min, f ) ); }

    // return true if lb < f < ub
    template < typename T >  bool  inside ( const T  lb,
                                            const T  f,
                                            const T  ub ) noexcept 
    {
        return (( lb < f ) && ( f < ub ));
    }

    // return true if lb ≤ f ≤ ub
    template < typename T >  bool  inside_eq ( const T  lb,
                                               const T  f,
                                               const T  ub ) noexcept 
    {
        return (( lb <= f ) && ( f <= ub ));
    }

    //
    // modulus, sign and roots
    //

    // return absolute value of argument
    inline int     abs     ( const int             f ) noexcept { return (f < 0 ? -f : f); }
    inline long    abs     ( const long            f ) noexcept { return (f < 0 ? -f : f); }
    inline float   abs     ( const float           f ) noexcept { return std::abs( f ); }
    inline double  abs     ( const double          f ) noexcept { return std::abs( f ); }
    inline float   abs     ( const Complex<float>  f ) noexcept { return f.abs(); }
    inline double  abs     ( const Complex<double> f ) noexcept { return f.abs(); }

    // return signum of argument
    inline float   sign    ( const float   f ) { return (f == 0.f ? 0.f : (f > 0.f ? 1.f : -1.f)); }
    inline double  sign    ( const double  f ) { return (f == 0.0 ? 0.0 : (f > 0.0 ? 1.0 : -1.0)); }

    inline Complex< float >   sign  ( const Complex< float >   f ) { return f.sign(); }
    inline Complex< double >  sign  ( const Complex< double >  f ) { return f.sign(); }

    // return square of argument
    inline float   square  ( const float   f ) { return f*f; }
    inline double  square  ( const double  f ) { return f*f; }

    inline Complex< float >   square  ( const Complex< float >   f ) { return f*f; }
    inline Complex< double >  square  ( const Complex< double >  f ) { return f*f; }

    // return square root of argument
    inline float   sqrt    ( const float   f ) { return std::sqrt( f ); }
    inline double  sqrt    ( const double  f ) { return std::sqrt( f ); }

    inline Complex< float >   sqrt  ( const Complex< float >   f ) { return f.sqrt(); }
    inline Complex< double >  sqrt  ( const Complex< double >  f ) { return f.sqrt(); }

    // return reciprocal square root of argument
    inline float   rsqrt   ( const float  f ) { return 1.0f / std::sqrt( f ); }
    inline double  rsqrt   ( const double f ) { return 1.0  / std::sqrt( f ); }

    //
    // power and logarithms
    //

    // compute x raised to the power of y
    inline float   pow     ( const float  x, const float  y ) { return std::pow( x, y ); }
    inline double  pow     ( const double x, const double y ) { return std::pow( x, y ); }

    // compute e raised to the power of x
    inline float   exp     ( const float   f ) { return std::exp( f ); }
    inline double  exp     ( const double  f ) { return std::exp( f ); }

    inline Complex< float >   exp  ( const Complex< float >   f ) { return HLIB::exp( f ); }
    inline Complex< double >  exp  ( const Complex< double >  f ) { return HLIB::exp( f ); }

    // compute base 2 logarithm of integers
    uint           log2    ( const uint n );

    // compute natural logarithm
    inline float   log     ( const float  f ) { return std::log( f ); }
    inline double  log     ( const double f ) { return std::log( f ); }

    // compute base 10 logarithm
    inline float   log10   ( const float  f ) { return std::log10( f ); }
    inline double  log10   ( const double f ) { return std::log10( f ); }

    //
    // rounding
    //

    // round argument down to the nearest integer
    inline float   floor   ( const float  f ) { return std::floor( f ); }
    inline double  floor   ( const double f ) { return std::floor( f ); }
    
    // round argument up to the nearest integer
    inline float   ceil    ( const float  f ) { return std::ceil( f ); }
    inline double  ceil    ( const double f ) { return std::ceil( f ); }
    
    //
    // trigonometry
    //

    // return sine and cosing of argument
    inline float   sin     ( const float  f ) { return std::sin( f ); }
    inline double  sin     ( const double f ) { return std::sin( f ); }
    inline float   cos     ( const float  f ) { return std::cos( f ); }
    inline double  cos     ( const double f ) { return std::cos( f ); }

    // return arc sine and cosing of argument
    inline float   asin    ( const float  f ) { return std::asin( f ); }
    inline double  asin    ( const double f ) { return std::asin( f ); }
    inline float   acos    ( const float  f ) { return std::acos( f ); }
    inline double  acos    ( const double f ) { return std::acos( f ); }

    // simultaneously compute sine and cosine
    inline void    sincos  ( const double f,
                             double &     s,
                             double &     c )
    {
#if USE_SVML == 1
        
        ::sincos( f, & s, & c );
        
#elif USE_AMDLIBM == 1
        
        amd_sincos( f, & s, & c );
        
#elif USE_ACML == 1
        
        fastsincos( f, & s, & c );
        
#elif HAS_SINCOS == 1
        
        ::sincos( f, & s, & c );
        
#else
        s = sin( f );
        c = cos( f );
#endif
    }

    //
    // check of values
    //

    // return true if given value contains Inf
    template <typename T> bool is_inf ( const T val );

    // return true if given value contains NaN
    template <typename T> bool is_nan ( const T val );

    //
    // misc.
    //

    // compute givens rotation (cs,sn) for given <a>, <b>
    void givens ( const float  a, const float  b, float &  cs, float &  sn );
    void givens ( const double a, const double b, double & cs, double & sn );

    void givens ( const Complex<float> &   a,
                  const Complex<float> &   b,
                  Complex<float> &         cs,
                  Complex<float> &         sn );
    void givens ( const Complex<double> &  a,
                  const Complex<double> &  b,
                  Complex<double> &        cs,
                  Complex<double> &        sn );

}// namespace Math

///////////////////////////////////////////////////
//
// information about types
//
///////////////////////////////////////////////////

namespace Limits
{
    
    // minimal representable value
    template <class T>  T  min      () { return std::numeric_limits< T >::min(); }
    
    // maximal representable value
    template <class T>  T  max      () { return std::numeric_limits< T >::max(); }
    
    // return difference between 1 and smallest value bigger than 1
    template <class T>  T  epsilon  () { return std::numeric_limits< T >::epsilon(); }

    // return NaN value
    template <class T>  T  nan      ();

}// namespace Limits

///////////////////////////////////////////////////
//
// Time related types and functions
//
///////////////////////////////////////////////////

namespace Time
{
    //
    // helper function to print timings
    //
    void  __autotimer_print ( const std::string &  format,
                              double               sec );

    //
    // stores time durations
    //
    struct TDuration
    {
        double  value;
    
        // default ctor
        TDuration ( const double  val = 0.0 ) : value(val) {}
    
        // convert to various units
        double  hour     () const { return value / 3600.0; }
        double  minute   () const { return value / 60.0; }
        double  seconds  () const { return value; }
        double  millisec () const { return value * 1000.0; }
        double  microsec () const { return value * 1000000.0; }
    
        // convert to double (seconds)
        operator double () const { return value; }
    
        // convert to string
        std::string  to_string     () const;

        // convert to string (extended version)
        std::string  to_string_ext () const;
    
        //! stream output
        friend std::ostream & operator << ( std::ostream & os, const TDuration & t )
        {
            return os << t.to_string();
        }
    };
    
    inline
    TDuration
    operator + ( const TDuration  d1,
                 const TDuration  d2 )
    {
        return d1.value + d2.value;
    }
    
    //
    // stores time points
    //
    template < int TIME_TYPE >
    struct TBaseTimePoint
    {
        double  value;
        
        // default ctor
        TBaseTimePoint ( const double  val = 0.0 ) : value(val) {}
    };
    
    template < int TIME_TYPE >
    inline TDuration
    operator - ( const TBaseTimePoint< TIME_TYPE >  t1,
                 const TBaseTimePoint< TIME_TYPE >  t2 )
    {
        return t1.value - t2.value;
    }
    
    enum TTimeType
    {
        PROCESS_TIME,
        THREAD_TIME,
        WALL_TIME
    };
    
    //! return CPU time since start of program
    double  cpu_time         ();

    //! return CPU time of current thread since thread start
    double  cpu_time_thread  ();

    //! return current wall time in seconds
    double  wall_time        ();

    //! measure and print time spent in current block
    #define DEF_AUTO_TIMER                                              \
    class TAutoTimer                                                    \
        {                                                               \
        private:                                                        \
            std::string  _format;                                       \
            TTimePoint   _tic;                                          \
                                                                        \
        public:                                                         \
            TAutoTimer ( const std::string &  aformat )                 \
                : _format( aformat )                                    \
                , _tic( now() )                                         \
            {}                                                          \
                                                                        \
            ~TAutoTimer ()                                              \
            {                                                           \
                auto  toc = since( _tic );                              \
                                                                        \
                __autotimer_print( _format, double( toc ) );            \
            }                                                           \
        }


namespace Process
{
    using TTimePoint = TBaseTimePoint< PROCESS_TIME >;

    //! return process time since start of program
    TTimePoint  now ();

    //! return process time since \a t
    TDuration   since ( const TTimePoint  t );

    DEF_AUTO_TIMER;

}// namespace CPU

namespace Thread
{
    using TTimePoint = TBaseTimePoint< THREAD_TIME >;

    //! return thread time since start of program
    TTimePoint  now ();

    //! return thread time since \a t
    TDuration   since ( const TTimePoint  t );

    DEF_AUTO_TIMER;

}// namespace Thread

namespace Wall
{
    using TTimePoint = TBaseTimePoint< WALL_TIME >;

    //! return current wall time
    TTimePoint  now ();

    //! return wall time since \a t
    TDuration   since ( const TTimePoint  t );

    DEF_AUTO_TIMER;

}// namespace Wall

}// namespace Time

///////////////////////////////////////////////////
//
// cryptographic functions
//
///////////////////////////////////////////////////

namespace Crypt
{

    //! compute CRC32 checksum of byte stream stored in \a data with \a size entries
    uint crc32 ( const size_t size, const void * data );

    //! template version if \see crc32 for automatic type handling
    template <typename T>
    uint crc32 ( const size_t nelem, const T * data )
    {
        return crc32( sizeof(T) * nelem, reinterpret_cast< const void * >( data ) );
    }

    //! compute Adler32 checksum of byte stream stored in \a data with \a size entries
    ulong adler32 ( const size_t size, const void * data );

}// namespace Crypt

///////////////////////////////////////////////////
//
// Support for RTTI (independent of C++)
//
///////////////////////////////////////////////////

namespace RTTI
{

    // convert a typename to a unique id
    uint         type_to_id        ( const std::string &  type );

    // return typename to a given id
    std::string  id_to_type        ( const uint           id );

    // register type in RTTI system (returns id)
    uint         register_type     ( const char *         type );

    // print all registered types with ids
    void         print_registered  ();

}// namespace RTTI

///////////////////////////////////////////////////
//
// memory management
//
///////////////////////////////////////////////////

namespace Mem
{

    //
    // return current memory usage in bytes
    //
    size_t       usage     ();

    //
    // return maximal memory usage of application so far in bytes
    //
    size_t       max_usage ();

    //
    // return string containing pretty printed current memory usage
    //
    std::string  to_string ();

    //
    // convert given number of bytes to human readable format
    //
    std::string  to_string ( const size_t bytes );

    //
    // wrapper around malloc for statistics (only for debugging)
    //
    void *       alloc     ( const size_t  size );
    
    //
    // wrapper around free for statistics (only for debugging)
    //
    void         free      ( void *        ptr );
    
}// namespace Mem

///////////////////////////////////////////////////
//
// machine properties
//
///////////////////////////////////////////////////

namespace Mach
{

    //
    // return string with list of CPUs associated to process
    //
    std::string  cpuset   ();

    //
    // return hostname
    //
    std::string  hostname ();

}// namespace Mach

}// namespace HLIB

#endif  // __HLIB_SYSTEM_HH
