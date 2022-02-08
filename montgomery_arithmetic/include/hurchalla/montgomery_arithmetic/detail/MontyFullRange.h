// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// The name "Fullrange" signifies that there are essentially no preconditions on
// the value of the modulus used in the Montgomery representation.


struct TagMontyFullrange final {};


// struct used internally by MontyFullRange
template <typename T>
struct MontyFRValueTypes {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    // regular montgomery value type
    struct V : public BaseMontgomeryValue<T> {
        HURCHALLA_FORCE_INLINE V() = default;
     protected:
        template <typename> friend class MontyFullRange;
        HURCHALLA_FORCE_INLINE explicit V(T a) : BaseMontgomeryValue<T>(a) {}
    };
    // canonical montgomery value type
    struct C : public V {
        HURCHALLA_FORCE_INLINE C() = default;
        HURCHALLA_FORCE_INLINE friend bool operator==(const C& x, const C& y)
            { return x.get() == y.get(); }
        HURCHALLA_FORCE_INLINE friend bool operator!=(const C& x, const C& y)
            { return !(x == y); }
     protected:
        template <template<class> class, template<class> class, typename>
          friend class MontyCommonBase;
        template <typename> friend class MontyFullRange;
        HURCHALLA_FORCE_INLINE explicit C(T a) : V(a) {}
    };
    // fusing montgomery value (addend/subtrahend for fmadd/fmsub)
    using FV = C;
};


// Let the theoretical constant R = 2^(ut_numeric_limits<T>::digits).
template <typename T>
class MontyFullRange final :
                  public MontyCommonBase<MontyFullRange, MontyFRValueTypes, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    using BC = MontyCommonBase<::hurchalla::detail::MontyFullRange,
                               ::hurchalla::detail::MontyFRValueTypes, T>;
    using BC::n_;
    using typename BC::V;
    using typename BC::C;
    using FV = typename MontyFRValueTypes<T>::FV;
 public:
    using MontyTag = TagMontyFullrange;
    using uint_type = T;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;

    using BC::fmadd;
    using BC::fmsub;

    explicit MontyFullRange(T modulus) : BC(modulus) {}
    MontyFullRange(const MontyFullRange&) = delete;
    MontyFullRange& operator=(const MontyFullRange&) = delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return (ut_numeric_limits<T>::max() % 2 == 0) ?
                static_cast<T>(ut_numeric_limits<T>::max() - 1) :
                ut_numeric_limits<T>::max();
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        // this static_assert guarantees 0 <= x.get()
        static_assert(!(ut_numeric_limits<decltype(x.get())>::is_signed), "");
        HPBC_PRECONDITION2(x.get() < n_);
        return C(x.get());
    }

    static_assert(std::is_same<FV, C>::value, "");
    // Note: fmsub and fmadd with FusingValue (FV) arguments will match to
    // fmsub and fmadd with CanonicalValue args, since C is_same as FV.

    HURCHALLA_FORCE_INLINE FV getFusingValue(V x) const
    {
        return getCanonicalValue(x);
    }

    using BC::add;
    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T result = modular_addition_prereduced_inputs(x.get(), y.get(), n_);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: add(V, C) will match to add(V x, V y) above

    using BC::subtract;
    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T result = modular_subtraction_prereduced_inputs(x.get(), y.get(), n_);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: subtract(V, C) will match to subtract(V x, V y) above
    // Note: subtract(C, V) will match to subtract(V x, V y) above

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T result = absolute_value_difference(x.get(), y.get());
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: unordered_subtract(C, V) and unordered_subtract(V, C) will match to
    // unordered_subtract(V x, V y) above

private:
    // functions called by the 'curiously recurring template pattern' base (BC).
    friend BC;

    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(bool& resultIsZero, T u_hi, T u_lo, PTAG) const
    {
        HPBC_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        T result = REDC_standard(u_hi, u_lo, n_, BC::inv_n_, PTAG());
        resultIsZero = (result == 0);
        HPBC_ASSERT2(result < n_);
        return V(result);
    }
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(T u_hi, T u_lo, PTAG) const
    {
        bool resultIsZero;  // ignored
        return montyREDC(resultIsZero, u_hi, u_lo, PTAG());
    }

    // return the high word of the product, and write the low word of the
    // product to u_lo.
    HURCHALLA_FORCE_INLINE T multiplyToHiLo(T& u_lo, V x, V y) const
    {
        return unsigned_multiply_to_hilo_product(u_lo, x.get(), y.get());
    }
    HURCHALLA_FORCE_INLINE T squareToHiLo(T& u_lo, V x) const
    {
        return unsigned_multiply_to_hilo_product(u_lo, x.get(), x.get());
    }
    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        return (x.get() < n_);
    }
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        return (0 <= x.get() && x.get() < n_);
    }
    // get a Natural number (i.e. number >= 0) congruent to x (mod n)
    HURCHALLA_FORCE_INLINE T getNaturalEquivalence(V x) const
    {
        return x.get();
    }
};


}} // end namespace

#endif
