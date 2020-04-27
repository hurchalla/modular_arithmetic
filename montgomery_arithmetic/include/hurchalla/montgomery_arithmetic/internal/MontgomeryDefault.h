
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H__INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H__INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/internal/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/internal/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/internal/sized_uint.h"
#include <type_traits>
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


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
