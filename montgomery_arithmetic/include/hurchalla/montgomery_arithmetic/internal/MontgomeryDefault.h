
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H__INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H__INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/internal/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/internal/MontySqrtRange.h"
#include <type_traits>
#include <limits>
#include <cstdint>

namespace hurchalla { namespace montgomery_arithmetic {


template <int BITS> struct sized_uint { using type = void; };
template<> struct sized_uint<8>  { using type = uint8_t; };
template<> struct sized_uint<16> { using type = uint16_t; };
template<> struct sized_uint<32> { using type = uint32_t; };
template<> struct sized_uint<64> { using type = uint64_t; };
// uint128_t would be nice, but it's unsupported on most compilers
//template<> struct sized_uint<128> { using type = uint128_t; 
//                                   static constexpr bool is_valid = true; };


// When this primary template is instantiated, useTypeTraits==true.
// In this case, we assume T is an integral type known to std::type_traits.
template <bool useTypeTraits, typename T>
class MontgomeryDefaultExtended {
    static_assert(std::is_integral<T>::value, "");
    using U = typename std::make_unsigned<T>::type;
    static constexpr int ubits = std::numeric_limits<U>::digits;
public:
    using type =
        typename std::conditional<
            !std::is_same<typename sized_uint<ubits*2>::type, void>::value
              && ubits*2 <= TARGET_BIT_WIDTH,
            MontySqrtRange<typename sized_uint<ubits*2>::type>,
            typename std::conditional<
                std::is_signed<T>::value,
                MontyHalfRange<U>,
                MontyFullRange<U>
            >::type
        >::type;
};

// When this specialization is instantiated, useTypeTraits==false.
// In this case, we assume U is not known to std::type_traits, and we require
// that U is an unsigned integral type. (Signed integral type for U is not
// possible because we would need to derive an unsigned integral type from it,
// and std::make_unsigned won't work for any type unknown to std::type_traits.)
template <typename U>
class MontgomeryDefaultExtended<false, U> {
    static_assert(std::numeric_limits<U>::is_integer, "");
    static_assert(!(std::numeric_limits<U>::is_signed), "");
    static constexpr int ubits = std::numeric_limits<U>::digits;
public:
    using type =
        typename std::conditional<
            !std::is_same<typename sized_uint<ubits*2>::type, void>::value
              && ubits*2 <= TARGET_BIT_WIDTH,
            MontySqrtRange<typename sized_uint<ubits*2>::type>,
            MontyFullRange<U>
        >::type;
};


template <typename T>
class MontgomeryDefault {
public:
    // we instantiate based on whether T is known to std::type_traits
    using type =
        typename MontgomeryDefaultExtended<std::is_integral<T>::value, T>::type;
};


}} // end namespace

#endif
