
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/detail/sized_uint.h"
#include "hurchalla/modular_arithmetic/traits/type_traits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include <type_traits>
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// If this primary template is instantiated, then T is an unsigned integral type
template <typename T, typename Enable = void>
class MontgomeryDefault {
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static constexpr int ubits = std::numeric_limits<T>::digits;
public:
    using type =
        typename std::conditional<
            !(std::is_same<typename sized_uint<ubits*2>::type, void>::value)
                && (ubits*2 <= HURCHALLA_TARGET_BIT_WIDTH),
            MontySqrtRange<typename sized_uint<ubits*2>::type>,
            MontyFullRange<T>
        >::type;
};

// If this partial specialization is instantiated, T is a signed integral type.
template <typename T>
class MontgomeryDefault<T, typename std::enable_if<
std::numeric_limits<T>::is_integer && std::numeric_limits<T>::is_signed>::type>
{
    using U = typename modular_arithmetic::extensible_make_unsigned<T>::type;
    static constexpr int ubits = std::numeric_limits<U>::digits;
public:
    using type =
        typename std::conditional<
            !(std::is_same<typename sized_uint<ubits*2>::type, void>::value)
                && (ubits*2 <= HURCHALLA_TARGET_BIT_WIDTH),
            MontySqrtRange<typename sized_uint<ubits*2>::type>,
            MontyHalfRange<U>
        >::type;
};


}} // end namespace

#endif
