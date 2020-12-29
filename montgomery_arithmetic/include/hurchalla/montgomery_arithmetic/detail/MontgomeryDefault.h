// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontySixthRange.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>

namespace hurchalla { namespace montgomery_arithmetic { namespace detail {


// If this primary template is instantiated, then T is an unsigned integral type
template <typename T, typename Enable = void>
class MontgomeryDefault {
    static_assert(util::ut_numeric_limits<T>::is_integer, "");
    static_assert(!(util::ut_numeric_limits<T>::is_signed), "");
    static constexpr int ubits = util::ut_numeric_limits<T>::digits;
public:
    using type = typename std::conditional<
                      util::sized_uint<ubits*2>::is_valid
                                     && (ubits*2 <= HURCHALLA_TARGET_BIT_WIDTH),
                      MontySixthRange<typename util::sized_uint<ubits*2>::type>,
                      MontyFullRange<T>
                 >::type;
};

// If this partial specialization is instantiated, T is a signed integral type.
template <typename T>
class MontgomeryDefault<T, typename std::enable_if<
                                 util::ut_numeric_limits<T>::is_integer &&
                                 util::ut_numeric_limits<T>::is_signed>::type>
{
    using U = typename util::extensible_make_unsigned<T>::type;
    static constexpr int ubits = util::ut_numeric_limits<U>::digits;
public:
    using type = typename std::conditional<
                      util::sized_uint<ubits*2>::is_valid
                                     && (ubits*2 <= HURCHALLA_TARGET_BIT_WIDTH),
                      MontySixthRange<typename util::sized_uint<ubits*2>::type>,
                      MontyHalfRange<U>
                 >::type;
};


}}} // end namespace

#endif
