// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_BASE_MONTGOMERY_VALUE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_BASE_MONTGOMERY_VALUE_H_INCLUDED


#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_unsigned.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
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

    template <typename T1 = T> HURCHALLA_FORCE_INLINE
    typename std::enable_if<(ut_numeric_limits<T1>::digits <=
                             HURCHALLA_TARGET_BIT_WIDTH)>::type 
    cmov(bool cond, BaseMontgomeryValue v)
    {
        HURCHALLA_CMOV(cond, value, v.value);
    }

    template <typename T1 = T> HURCHALLA_FORCE_INLINE
    typename std::enable_if<(ut_numeric_limits<T1>::digits >
                             HURCHALLA_TARGET_BIT_WIDTH)>::type 
    cmov(bool cond, BaseMontgomeryValue v)
    {
        value = cond ? v.value : value;
    }

    HURCHALLA_FORCE_INLINE void cmov_masked(bool cond, BaseMontgomeryValue v)
    {
        using U = typename extensible_make_unsigned<T>::type;
        using P = typename safely_promote_unsigned<U>::type;
        P mask = static_cast<P>(static_cast<P>(0) - static_cast<P>(cond));
        P maskflip = static_cast<P>(static_cast<P>(cond) - static_cast<P>(1));
        value = static_cast<T>((mask & static_cast<P>(v.value)) |
                               (maskflip & static_cast<P>(value)));
    }
};


}} // end namespace

#endif
