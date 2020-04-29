
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/internal/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/internal/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/internal/sized_uint.h"
#include <type_traits>
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// If this primary template is instantiated, then T is not known to
// std::type_traits.  In this case we require that T is an unsigned integral
// type.  (Signed integral type for T is not possible because we would need to
// derive an unsigned integral type from it, and std::make_unsigned won't work
// for any type unknown to std::type_traits.)
template <typename T, typename Enable = void>
class MontgomeryDefault {
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static constexpr int ubits = std::numeric_limits<T>::digits;
public:
    using type =
        typename std::conditional<
            !std::is_same<typename sized_uint<ubits*2>::type, void>::value
              && ubits*2 <= TARGET_BIT_WIDTH,
            MontySqrtRange<typename sized_uint<ubits*2>::type>,
            MontyFullRange<T>
        >::type;
};

// If this partial specialization is instantiated, then T is an integral type
// known to std::type_traits (the class needs to safely use std::make_unsigned).
template <typename T>
class MontgomeryDefault<T,
                    typename std::enable_if<std::is_integral<T>::value>::type> {
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


}} // end namespace

#endif
