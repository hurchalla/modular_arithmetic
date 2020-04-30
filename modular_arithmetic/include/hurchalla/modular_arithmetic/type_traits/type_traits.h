
#ifndef HURCHALLA_MODULAR_ARITHMETIC_TYPE_TRAITS_H_INCLUDED
#define HURCHALLA_MODULAR_ARITHMETIC_TYPE_TRAITS_H_INCLUDED


#include <type_traits>

namespace hurchalla { namespace modular_arithmetic {


// primary template
template <typename T>
struct extensible_make_unsigned {
    static_assert(std::is_integral<T>::value,
                     "You'll need to specialize this template for your type T");
    using type = typename std::make_unsigned<T>::type;
};

// primary template
template <typename T>
struct extensible_make_signed {
    static_assert(std::is_integral<T>::value,
                     "You'll need to specialize this template for your type T");
    using type = typename std::make_signed<T>::type;
};

// ---Specializations---
//
// You can use the following as examples of how to specialize for a particular
// integer type you wish to use.  Note that the primary template will pass its
// static assert and work fine if its type T is known to std::type_traits.

#if defined(__SIZEOF_INT128__)
// The macro __SIZEOF_INT128__ indicates if __int128 is supported.  See
//https://stackoverflow.com/questions/16088282/is-there-a-128-bit-integer-in-gcc

template<>  // explicit specialization for T = __int128_t
struct extensible_make_unsigned<__int128_t> {
    using type = __uint128_t;
};
template<>  // explicit specialization for T = __uint128_t
struct extensible_make_unsigned<__uint128_t> {
    using type = __uint128_t;
};

template<>  // explicit specialization for T = __int128_t
struct extensible_make_signed<__int128_t> {
    using type = __int128_t;
};
template<>  // explicit specialization for T = __uint128_t
struct extensible_make_signed<__uint128_t> {
    using type = __int128_t;
};
#endif


}}  // end namespace

#endif
