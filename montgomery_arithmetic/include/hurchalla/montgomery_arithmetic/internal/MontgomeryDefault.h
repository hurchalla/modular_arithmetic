
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/internal/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/internal/MontySqrtRange.h"
#include "hurchalla/montgomery_arithmetic/internal/sized_uint.h"
#include "hurchalla/modular_arithmetic/type_traits/type_traits.h"
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
                && (ubits*2 <= TARGET_BIT_WIDTH),
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
                && (ubits*2 <= TARGET_BIT_WIDTH),
            MontySqrtRange<typename sized_uint<ubits*2>::type>,
            MontyHalfRange<U>
        >::type;
};


}} // end namespace

#endif
