// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_BASE_MONTGOMERY_VALUE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_BASE_MONTGOMERY_VALUE_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>

namespace hurchalla { namespace detail {


template <typename T>
class BaseMontgomeryValue {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    T value;
protected:
    HURCHALLA_FORCE_INLINE explicit BaseMontgomeryValue(T a) : value(a) {}
    HURCHALLA_FORCE_INLINE T get() const { return value; }
public:
    // This next constructor purposely does not initialize 'value' - the
    // contents are undefined until the object is assigned to.
    HURCHALLA_FORCE_INLINE BaseMontgomeryValue() = default;

    template <class PerfTag = CSelectDefaultTag>
    HURCHALLA_FORCE_INLINE void cmov(bool cond, BaseMontgomeryValue v)
    {
          // value = cond ? v.value : value
        value = conditional_select<T, PerfTag>(cond, v.value, value);
    }
};


}} // end namespace

#endif
