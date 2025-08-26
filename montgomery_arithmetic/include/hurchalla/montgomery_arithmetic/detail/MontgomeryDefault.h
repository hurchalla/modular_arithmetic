// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTGOMERY_DEFAULT_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyFullRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyHalfRange.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyQuarterRange.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/sized_uint.h"
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
                     typename std::conditional<
                         (bitsT <= target_bits - 1),
                         MontyHalfRange<typename sized_uint<target_bits>::type>,
                         MontyFullRange<U>
                     >::type
                 >::type;
};

// Implementation note: when bitsT > target_bits (e.g. T == __int128_t on a 64
// bit system), we purposely never use MontyHalfRange above and instead default
// to MontyFullRange, because MontyFullRange uses unsigned hi_lo mults, whereas
// MontyHalfRange uses signed hi_lo multiplications...
// When bitsT > target_bits we're forced to use a 'slow' hi_lo mult routine,
// since there's no simple asm instruction that's applicable- e.g. on x86_64,
// we need far more than a single MUL or IMUL.  And unfortunately we don't have
// a signed routine that's as good as unsigned when bitsT > target_bits.  For
// details see the comments for slow_signed_multiply_to_hilo_product() in
// hurchalla/util/detail/platform_specific/impl_signed_multiply_to_hilo_product.h


}} // end namespace

#endif
