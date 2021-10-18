// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/monty_tag_structs.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// Let the theoretical constant R = 2^(ut_numeric_limits<T>::digits).
template <typename T>
class MontyFullRange final : public MontyCommonBase<MontyFullRange, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    using BC = MontyCommonBase<::hurchalla::detail::MontyFullRange, T>;
    using BC::n_;
    using typename BC::W;
    using typename BC::V;
    using typename BC::C;
 public:
    using widevalue_type = W;
    using montvalue_type = V;
    using canonvalue_type = C;
//    using squaringvalue_type = SQ;
    using uint_type = T;
    using MontyTag = FullrangeTag;

    explicit MontyFullRange(T modulus) : BC(modulus) {}
    MontyFullRange(const MontyFullRange&) = delete;
    MontyFullRange& operator=(const MontyFullRange&) = delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return (ut_numeric_limits<T>::max() % 2 == 0) ?
                static_cast<T>(ut_numeric_limits<T>::max() - 1) :
                ut_numeric_limits<T>::max();
    }

    HURCHALLA_FORCE_INLINE T getExtendedModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(W x) const
    {
        // this static_assert guarantees 0 <= x.get()
        static_assert(!(ut_numeric_limits<decltype(x.get())>::is_signed), "");
        HPBC_PRECONDITION2(x.get() < n_);
        return C(x.get());
    }
};


}} // end namespace

#endif
