// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyTags.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/absolute_value_difference.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// The name "Fullrange" signifies that there are essentially no preconditions on
// the value of the modulus used in the Montgomery representation.


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
    struct FV : public V {
        HURCHALLA_FORCE_INLINE FV() = default;
     protected:
        template <typename> friend class MontyFullRange;
        HURCHALLA_FORCE_INLINE explicit FV(T a) : V(a) {}
    };
    // squaring value type - used for square() optimizations (fyi, those
    // optimizations wouldn't help much or at all for monty types other than
    // MontyFullRange).
    struct SV {
        HURCHALLA_FORCE_INLINE SV() = default;
     protected:
        template <typename> friend class MontyFullRange;
        HURCHALLA_FORCE_INLINE T getbits() const { return bits; }
        HURCHALLA_FORCE_INLINE T get_subtrahend() const { return subtrahend; }
        HURCHALLA_FORCE_INLINE SV(T bbits, T subt) :
                                                bits(bbits), subtrahend(subt) {}
     private:
        T bits;
        T subtrahend;
    };
};


// Let the theoretical constant R = 1<<(ut_numeric_limits<T>::digits).
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
    using SV = typename MontyFRValueTypes<T>::SV;
 public:
    using MontyTag = TagMontyFullrange;
    using uint_type = T;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;
    using squaringvalue_type = SV;

    explicit MontyFullRange(T modulus) : BC(modulus) {}

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return (ut_numeric_limits<T>::max() % 2 == 0) ?
                static_cast<T>(ut_numeric_limits<T>::max() - 1) :
                ut_numeric_limits<T>::max();
    }

    HURCHALLA_FORCE_INLINE V negate(V x) const
    {
        return subtract(BC::getZeroValue(), x, LowuopsTag());
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        // this static_assert guarantees 0 <= x.get()
        static_assert(!(ut_numeric_limits<decltype(x.get())>::is_signed), "");
        HPBC_CLOCKWORK_PRECONDITION2(x.get() < n_);
        return C(x.get());
    }

    // Note: internal to MontyFullRange, the contents of FusingValue (FV) and
    // CanonicalValue (C) variables are interchangeable.  Other Monty types
    // use FV and C as completely distinct types, and so for genericity we
    // always present C and FV to the outside world as being unrelated.
    HURCHALLA_FORCE_INLINE FV getFusingValue(V x) const
    {
        C cv = getCanonicalValue(x);
        return FV(cv.get());
    }
    using BC::fmadd;
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, FV fv, PTAG) const
    {
        C cv = C(fv.get());
        return fmadd(x, y, cv, PTAG());
    }
    using BC::fmsub;
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, FV fv, PTAG) const
    {
        C cv = C(fv.get());
        return fmsub(x, y, cv, PTAG());
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isValid(x));
        HPBC_CLOCKWORK_PRECONDITION2(isValid(y));
        namespace hc = ::hurchalla;
        T result = hc::modular_addition_prereduced_inputs(x.get(), y.get(), n_);
        HPBC_CLOCKWORK_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: add(V, C) will match to add(V x, V y) above
    HURCHALLA_FORCE_INLINE C add(C cx, C cy) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(cx.get() < n_);
        HPBC_CLOCKWORK_PRECONDITION2(cy.get() < n_);
        namespace hc = ::hurchalla;
        T result = hc::modular_addition_prereduced_inputs(cx.get(),cy.get(),n_);
        HPBC_CLOCKWORK_POSTCONDITION2(result < n_);
        return C(result);
    }

    template <class PTAG>
    HURCHALLA_FORCE_INLINE V subtract(V x, V y, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isValid(x));
        HPBC_CLOCKWORK_PRECONDITION2(isValid(y));
        T result = ::hurchalla::modular_subtraction_prereduced_inputs
                                      <decltype(n_),PTAG>(x.get(), y.get(), n_);
        HPBC_CLOCKWORK_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: subtract(V, C, PTAG) will match to subtract(V x, V y, PTAG) above
    // Note: subtract(C, V, PTAG) will match to subtract(V x, V y, PTAG) above
    template <class PTAG>
    HURCHALLA_FORCE_INLINE C subtract(C cx, C cy, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(cx.get() < n_);
        HPBC_CLOCKWORK_PRECONDITION2(cy.get() < n_);
        T result = ::hurchalla::modular_subtraction_prereduced_inputs
                                    <decltype(n_),PTAG>(cx.get(), cy.get(), n_);
        HPBC_CLOCKWORK_POSTCONDITION2(result < n_);
        return C(result);
    }

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isValid(x));
        HPBC_CLOCKWORK_PRECONDITION2(isValid(y));
        T result = ::hurchalla::absolute_value_difference(x.get(), y.get());
        HPBC_CLOCKWORK_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: unordered_subtract(C, V) and unordered_subtract(V, C) will match to
    // unordered_subtract(V x, V y) above

    HURCHALLA_FORCE_INLINE V two_times(V x) const
    {
        return add(x, x);
    }
    HURCHALLA_FORCE_INLINE C two_times(C x) const
    {
        return add(x, x);
    }


    HURCHALLA_FORCE_INLINE SV getSquaringValue(V x) const
    {
        return SV(x.get(), 0);
    }

    template <class PTAG> HURCHALLA_FORCE_INLINE
    SV squareSV(SV sv, PTAG) const
    {
        // see squareToHiLo in MontyFullRangeMasked.h for basic ideas of
        // proof for why u_hi and u_lo are correct
        namespace hc = ::hurchalla;
        T a = sv.getbits();
        T sqlo;
        T sqhi = hc::unsigned_multiply_to_hilo_product(sqlo, a, a);
        T a_or_zero = sv.get_subtrahend();
        T u_hi = static_cast<T>(sqhi - a_or_zero - a_or_zero);
        T u_lo = sqlo;
        HPBC_CLOCKWORK_ASSERT2(u_hi < n_);

        T minu, subt;
        hc::REDC_incomplete(minu, subt, u_hi, u_lo, n_, BC::inv_n_, PTAG());
        T res = static_cast<T>(minu - subt);
        T subtrahend = (minu < subt) ? res : static_cast<T>(0);
        SV result(res, subtrahend);
        return result;
    }

    template <class PTAG> HURCHALLA_FORCE_INLINE
    V squareToMontgomeryValue(SV sv, PTAG) const
    {
        // see squareToHiLo in MontyFullRangeMasked.h for basic ideas of
        // proof for why u_hi and u_lo are correct
        namespace hc = ::hurchalla;
        T a = sv.getbits();
        T sqlo;
        T sqhi = hc::unsigned_multiply_to_hilo_product(sqlo, a, a);
        T a_or_zero = sv.get_subtrahend();
        T u_hi = static_cast<T>(sqhi - a_or_zero - a_or_zero);
        T u_lo = sqlo;
        HPBC_CLOCKWORK_ASSERT2(u_hi < n_);

        T res = hc::REDC_standard(u_hi, u_lo, n_, BC::inv_n_, PTAG());
        V result(res);
        return result;
    }

    // probably I would not want to use this, instead preferring to get a SV
    // via squareToMontgomeryValue
    HURCHALLA_FORCE_INLINE V getMontgomeryValue(SV sv) const
    {
        T nonneg_value = sv.get_subtrahend() != 0 ? sv.get() + n_ : sv.get();
        return V(nonneg_value);
    }

private:
    // functions called by the 'curiously recurring template pattern' base (BC).
    friend BC;

    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(bool& resultIsZero, T u_hi, T u_lo, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        namespace hc = ::hurchalla;
        T result = hc::REDC_standard(u_hi, u_lo, n_, BC::inv_n_, PTAG());
        resultIsZero = (result == 0);
        HPBC_CLOCKWORK_ASSERT2(result < n_);
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
        namespace hc = ::hurchalla;
        return hc::unsigned_multiply_to_hilo_product(u_lo, x.get(), y.get());
    }
    HURCHALLA_FORCE_INLINE T squareToHiLo(T& u_lo, V x) const
    {
        namespace hc = ::hurchalla;
        return hc::unsigned_multiply_to_hilo_product(u_lo, x.get(), x.get());
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
