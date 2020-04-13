
#ifndef HURCHALLA_MODULAR_MONTGOMERY_MONTGOMERY_DEFAULT_H__INCLUDED
#define HURCHALLA_MODULAR_MONTGOMERY_MONTGOMERY_DEFAULT_H__INCLUDED


#include "hurchalla/modular_arithmetic/montgomery/internal/montyfullrange.h"
#include "hurchalla/modular_arithmetic/montgomery/internal/montyhalfrange.h"
#include "hurchalla/modular_arithmetic/montgomery/internal/montysqrtrange.h"
#include <type_traits>
#include <cstdint>

namespace hurchalla { namespace modular_montgomery {


template <int M> struct sized_uint { using type = void; };
template<> struct sized_uint<1> { using type = uint8_t; };
template<> struct sized_uint<2> { using type = uint16_t; };
template<> struct sized_uint<4> { using type = uint32_t; };
template<> struct sized_uint<8> { using type = uint64_t; };
// uint128_t would be nice, but it's unsupported on most compilers
//template<> struct sized_uint<16> { using type = uint128_t; 
//                                   static constexpr bool is_valid = true; };

template <typename T>
class MontgomeryDefault {
    static_assert(std::is_integral<T>::value, "");
    using U = typename std::make_unsigned<T>::type;
public:
    using type = typename std::conditional<
         !std::is_same<typename sized_uint<sizeof(T)*2>::type, void>::value
             && sizeof(T)*2*8 <= TARGET_BIT_WIDTH,
         MontySqrtRange<typename sized_uint<sizeof(T)*2>::type>,
         typename std::conditional<
             std::is_signed<T>::value,
             MontyHalfRange<U>,
             MontyFullRange<U>
             >::type
         >::type;
};


}} // end namespace

#endif
