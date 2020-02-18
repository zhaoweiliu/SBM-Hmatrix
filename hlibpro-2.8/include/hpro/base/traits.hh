#ifndef __HLIB_BLAS_TRAITS_HH
#define __HLIB_BLAS_TRAITS_HH
//
// Project     : HLib
// File        : traits.hh
// Description : provide type traits for BLAS functions
// Author      : Ronald Kriemann
// Copyright   : Max Planck Institute MIS 2004-2018. All Rights Reserved.
//

#include "hpro/base/types.hh"

namespace HLIB
{

//!
//! signals integer types
//!
template <typename T> struct is_integer                   { static constexpr bool value = false; };
template <>           struct is_integer< char >           { static constexpr bool value = true; };
template <>           struct is_integer< unsigned char >  { static constexpr bool value = true; };
template <>           struct is_integer< short >          { static constexpr bool value = true; };
template <>           struct is_integer< unsigned short > { static constexpr bool value = true; };
template <>           struct is_integer< int >            { static constexpr bool value = true; };
template <>           struct is_integer< unsigned int >   { static constexpr bool value = true; };
template <>           struct is_integer< long >           { static constexpr bool value = true; };
template <>           struct is_integer< unsigned long >  { static constexpr bool value = true; };

//!
//! signals floating point types
//!
template <typename T> struct is_float           { static constexpr bool value = false; };
template <>           struct is_float< float >  { static constexpr bool value = true; };
template <>           struct is_float< double > { static constexpr bool value = true; };

//!
//! signals complex valued types
//!
template <typename T> struct is_complex_type                       { static constexpr bool value = false; };
template <>           struct is_complex_type< Complex< float > >   { static constexpr bool value = true; };
template <>           struct is_complex_type< Complex< double > >  { static constexpr bool value = true; };

//!
//! signals single/double precision
//!
template <typename T> struct is_single_prec                       { static constexpr bool value = false; };
template <>           struct is_single_prec< float >              { static constexpr bool value = true; };
template <>           struct is_single_prec< Complex< float > >   { static constexpr bool value = true; };

template <typename T> struct is_double_prec                       { static constexpr bool value = false; };
template <>           struct is_double_prec< double >             { static constexpr bool value = true; };
template <>           struct is_double_prec< Complex< double > >  { static constexpr bool value = true; };

//!
//! provide real valued type forming base of T
//!
template <typename T> struct real_type                  { using  type_t = T; };
template <typename T> struct real_type< Complex< T > >  { using  type_t = T; };

//!
//! convert from C++ type to value_type_t
//!
template <typename T> struct value_type                  { static constexpr value_type_t value = real_valued; };
template <typename T> struct value_type< Complex< T > >  { static constexpr value_type_t value = complex_valued; };

//
// test if T1 and T2 are of the same type
//
template < typename T1, typename T2 > struct is_same_type           { static constexpr bool  value = false; };
template < typename T1 >              struct is_same_type< T1, T1 > { static constexpr bool  value = true;  };

//
// enable functions only, if template argument evaluates to true
//
template < bool >  struct enable_if         { };
template <>        struct enable_if< true > { using  result = void; };

// same, but with given return value
template < bool, typename T_result >  struct enable_if_res                   { };
template <       typename T_result >  struct enable_if_res< true, T_result > { using  result = T_result; };

}// namespace HLIB

#endif  // __HLIB_BLAS_TRAITS_HH
