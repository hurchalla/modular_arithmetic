// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>

namespace hurchalla { namespace montgomery_arithmetic { namespace detail {


template <typename T>
class MontgomeryDefault {
    static_assert(util::ut_numeric_limits<T>::is_integer, "");
    using U = typename util::extensible_make_unsigned<T>::type;
    static constexpr int bitsT = util::ut_numeric_limits<T>::digits;
    static constexpr int target_bits = HURCHALLA_TARGET_BIT_WIDTH;
public:
    using type = typename std::conditional<
                (bitsT <= target_bits - 2),
                MontyQuarterRange<typename util::sized_uint<target_bits>::type>,
                MontyFullRange<U>
               >::type;
};


}}} // end namespace

#endif
