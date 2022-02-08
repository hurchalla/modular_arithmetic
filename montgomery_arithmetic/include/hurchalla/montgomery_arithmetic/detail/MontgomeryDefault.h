// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>

namespace hurchalla { namespace detail {


template <typename T>
class MontgomeryDefault final {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    using U = typename extensible_make_unsigned<T>::type;
    static constexpr int bitsT = ut_numeric_limits<T>::digits;
    static constexpr int target_bits = HURCHALLA_TARGET_BIT_WIDTH;
public:
    using type = typename std::conditional<
                     (bitsT <= target_bits - 2),
                     MontyQuarterRange<typename sized_uint<target_bits>::type>,
                     MontyFullRange<U>
                 >::type;
};


}} // end namespace

#endif
