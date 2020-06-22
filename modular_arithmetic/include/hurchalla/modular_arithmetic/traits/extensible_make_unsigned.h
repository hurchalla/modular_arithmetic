// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MODULAR_ARITHMETIC_EXTENSIBLE_MAKE_UNSIGNED_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_EXTENSIBLE_MAKE_UNSIGNED_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


// primary template
template <typename T>
struct extensible_make_unsigned {
    static_assert(std::is_integral<T>::value,
          "You'll need to specialize extensible_make_unsigned for your type T");
    using type = typename std::make_unsigned<T>::type;
};

// ---Specializations---
//
// You can use the following as examples of how to specialize for a particular
// integer type you wish to use.  Note that the primary template will pass its
// static assert and work fine if its type T is known to std::type_traits.
//
// ** If you specialize this type, you should probably also specialize
// extensible_make_signed (in extensible_make_signed.h) **

#if (HURCHALLA_COMPILER_HAS_UINT128_T())
template<>  // explicit specialization for T = __int128_t
struct extensible_make_unsigned<__int128_t> {
    using type = __uint128_t;
};
template<>  // explicit specialization for T = __uint128_t
struct extensible_make_unsigned<__uint128_t> {
    using type = __uint128_t;
};
#endif


}}  // end namespace

#endif
