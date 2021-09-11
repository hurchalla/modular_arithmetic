// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// Let the theoretical constant R = 2^(ut_numeric_limits<T>::digits).
template <typename T>
class MontyQuarterRange : public MontyCommonBase<MontyQuarterRange, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    static_assert(ut_numeric_limits<T>::digits >= 2, "");
    using BC = MontyCommonBase<::hurchalla::detail::MontyQuarterRange, T>;
    using BC::n_;
    using BC::inv_n_;
    using V = typename BC::MontgomeryValue;
public:
    using montvalue_type = V;
    using uint_type = T;
    using MontyTag = QuarterrangeTag;

    explicit MontyQuarterRange(T modulus) : BC(modulus)
    {
        // MontyQuarterRange requires  modulus < R/4
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(modulus < Rdiv4);
    }
    MontyQuarterRange(const MontyQuarterRange&) = delete;
    MontyQuarterRange& operator=(const MontyQuarterRange&) = delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return static_cast<T>((static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 2)) - 1);
    }

    HURCHALLA_FORCE_INLINE T getExtendedModulus() const
    {
        return static_cast<T>(2*n_);
    }

    HURCHALLA_FORCE_INLINE V getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        T cv = static_cast<T>(x.get() - n_);
#if 0
        cv = (x.get() < n_) ? x.get() : cv;
#else
        HURCHALLA_CMOV(x.get() < n_, cv, x.get());
#endif
        HPBC_POSTCONDITION2(cv < n_);
        return V(cv);
    }
};


}} // end namespace

#endif
