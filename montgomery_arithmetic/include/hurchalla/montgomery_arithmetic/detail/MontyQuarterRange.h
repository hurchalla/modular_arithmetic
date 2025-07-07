// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/quarterrange_get_canonical.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/absolute_value_difference.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_signed.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// The name "Quarterrange" signifies that the modulus must be less than R/4,
// where  R = 1<<(ut_numeric_limits<T>::digits).  For example, if T is uint64_t
// then R = 1<<64 and R/4 == 1<<62, and thus it would require  modulus < 1<<62.

// The MontyQuarterRange class functions require/allow an unusual input range:
// for an input x, they allow  0 <= x < 2*n, where n is the modulus.  Similarly,
// the return value range will be  0 <= returnValue < 2*n.  Obviously neither
// inputs nor outputs necessarily belong to the minimal residue class modulo n -
// i.e. they might not be fully reduced, modulo n.  Note that the algorithm for
// montgomery REDC requires that  u = x*y < n*R;  this will always be satisfied
// for any multiplication x*y of Quarterrange montgomery values.  To see why,
// keep in mind that Quarterrange requires n < R/4  and that all inputs are less
// than 2*n.  Thus the multiplication
// u = x*y < (2*n)*(2*n) == (4*n)*n < (4*n)*(R/4) == n*R,  which means u < n*R,
// as required.  For more details, see also section 5 from
// "Montgomery's Multiplication Technique: How to Make It Smaller and Faster"
// https://www.comodo.com/resources/research/cryptography/CDW_CHES_99.ps


struct TagMontyQuarterrange final {};  // IDs MontyQuarterRange independent of T


// struct used internally by MontyQuarterRange
template <typename T>
struct MontyQRValueTypes {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    // regular montgomery value type
    struct V : public BaseMontgomeryValue<T> {
        HURCHALLA_FORCE_INLINE V() = default;
     protected:
        template <typename> friend class MontyQuarterRange;
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
        template <typename> friend class MontyQuarterRange;
        HURCHALLA_FORCE_INLINE explicit C(T a) : V(a) {}
    };
    // fusing montgomery value (addend/subtrahend for fmadd/fmsub)
    struct FV : public V {
        HURCHALLA_FORCE_INLINE FV() = default;
     protected:
        template <typename> friend class MontyQuarterRange;
        HURCHALLA_FORCE_INLINE explicit FV(T a) : V(a) {}
    };
};


// Let the theoretical constant R = 1<<(ut_numeric_limits<T>::digits).
template <typename T>
class MontyQuarterRange final : public
                      MontyCommonBase<MontyQuarterRange, MontyQRValueTypes, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    static_assert(ut_numeric_limits<T>::digits >= 2, "");
    using BC = MontyCommonBase<::hurchalla::detail::MontyQuarterRange,
                               ::hurchalla::detail::MontyQRValueTypes, T>;
    using BC::n_;
    using typename BC::V;
    using typename BC::C;
    using FV = typename MontyQRValueTypes<T>::FV;

    using S = typename extensible_make_signed<T>::type;
    static_assert(static_cast<S>(-1) == ~(static_cast<S>(0)),
                                  "S must use two's complement representation");
    static_assert(static_cast<S>(static_cast<T>(static_cast<S>(-1))) ==
                  static_cast<S>(-1), "Casting a signed S value to unsigned and"
                               " back again must result in the original value");
 public:
    using MontyTag = TagMontyQuarterrange;
    using uint_type = T;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;

    explicit MontyQuarterRange(T modulus) : BC(modulus)
    {
        // MontyQuarterRange requires  modulus < R/4
        constexpr T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(modulus < Rdiv4);
    }

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return static_cast<T>((static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 2)) - 1);
    }

    HURCHALLA_FORCE_INLINE V negate(V x) const
    {
        return subtract(BC::getZeroValue(), x);
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(0 <= x.get() && x.get() < static_cast<T>(2*n_));
        T result = quarterrange_get_canonical<T>::call(x.get(), n_);
        HPBC_POSTCONDITION2(0 <= result && result < n_);
        return C(result);
    }

    // Note: internal to MontyQuarterRange, the contents of FusingValue (FV) and
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
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        // a type V value can be up to 2*n_, which requires digits-1 bits.
        constexpr int max_bits_needed = ut_numeric_limits<T>::digits - 1;
        constexpr T limit = static_cast<T>(static_cast<T>(1)<<max_bits_needed);
        T n2 = static_cast<T>(2*n_);
        HPBC_ASSERT2(x.get() < limit);
        HPBC_ASSERT2(y.get() < limit);
        HPBC_ASSERT2(n2 < limit);
        S sx = static_cast<S>(x.get());
        S sy = static_cast<S>(y.get());
        S sn2 = static_cast<S>(n2);
        S modsum = ::hurchalla::modular_addition_prereduced_inputs(sx, sy, sn2);
        HPBC_ASSERT2(modsum >= 0);
        T result = static_cast<T>(modsum);
        HPBC_POSTCONDITION2(result < n2);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V add(V x, C cy) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(cy.get() < n_);
        C cx = getCanonicalValue(x);
        T result = static_cast<T>(cx.get() + cy.get());
        HPBC_ASSERT2(result < 2*n_);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE C add(C cx, C cy) const
    {
        HPBC_PRECONDITION2(cx.get() < n_);
        HPBC_PRECONDITION2(cy.get() < n_);
        // a type C value has range [0,n), which requires digits-2 bits.
        constexpr int max_bits_needed = ut_numeric_limits<T>::digits - 2;
        constexpr T limit = static_cast<T>(static_cast<T>(1)<<max_bits_needed);
        HPBC_ASSERT2(cx.get() < limit);
        HPBC_ASSERT2(cy.get() < limit);
        HPBC_ASSERT2(n_ < limit);
        S sx = static_cast<S>(cx.get());
        S sy = static_cast<S>(cy.get());
        S sn = static_cast<S>(n_);
        S modsum = ::hurchalla::modular_addition_prereduced_inputs(sx, sy, sn);
        HPBC_ASSERT2(modsum >= 0);
        T result = static_cast<T>(modsum);
        HPBC_POSTCONDITION2(result < n_);
        return C(result);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        // a type V value can be up to 2*n, which requires digits-1 bits.
        constexpr int max_bits_needed = ut_numeric_limits<T>::digits - 1;
        constexpr T limit = static_cast<T>(static_cast<T>(1)<<max_bits_needed);
        T n2 = static_cast<T>(2*n_);
        HPBC_ASSERT2(x.get() < limit);
        HPBC_ASSERT2(y.get() < limit);
        HPBC_ASSERT2(n2 < limit);
        S sx = static_cast<S>(x.get());
        S sy = static_cast<S>(y.get());
        S sn2 = static_cast<S>(n2);
        namespace hc = ::hurchalla;
        S moddiff = hc::modular_subtraction_prereduced_inputs(sx, sy, sn2);
        HPBC_ASSERT2(moddiff >= 0);
        T result = static_cast<T>(moddiff);
        HPBC_POSTCONDITION2(result < n2);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE C subtract(C cx, C cy) const
    {
        HPBC_PRECONDITION2(cx.get() < n_);
        HPBC_PRECONDITION2(cy.get() < n_);
        // a type C value has range [0,n), which requires digits-2 bits.
        constexpr int max_bits_needed = ut_numeric_limits<T>::digits - 2;
        constexpr T limit = static_cast<T>(static_cast<T>(1)<<max_bits_needed);
        HPBC_ASSERT2(cx.get() < limit);
        HPBC_ASSERT2(cy.get() < limit);
        HPBC_ASSERT2(n_ < limit);
        S sx = static_cast<S>(cx.get());
        S sy = static_cast<S>(cy.get());
        S sn = static_cast<S>(n_);
        namespace hc = ::hurchalla;
        S moddiff = hc::modular_subtraction_prereduced_inputs(sx, sy, sn);
        HPBC_ASSERT2(moddiff >= 0);
        T result = static_cast<T>(moddiff);
        HPBC_POSTCONDITION2(result < n_);
        return C(result);
    }
    // Note: subtract(C,V) and subtract(V,C) will match to subtract(V,V) above

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        // a type V value can be up to 2*n_, which requires digits-1 bits.
        constexpr int max_bits_needed = ut_numeric_limits<T>::digits - 1;
        constexpr T limit = static_cast<T>(static_cast<T>(1)<<max_bits_needed);
        HPBC_ASSERT2(x.get() < limit);
        HPBC_ASSERT2(y.get() < limit);
        S sx = static_cast<S>(x.get());
        S sy = static_cast<S>(y.get());
        S absdiff = ::hurchalla::absolute_value_difference(sx, sy);
        HPBC_ASSERT2(absdiff >= 0);
        T result = static_cast<T>(absdiff);
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    // Note: unordered_subtract(C, V) and unordered_subtract(V, C) and
    // unordered_subtract(C, C) all match to unordered_subtract(V x, V y) above.

private:
    // functions called by the 'curiously recurring template pattern' base (BC).
    friend BC;

    // note: PTAG gets ignored here
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(bool& resultIsZero, T u_hi, T u_lo, PTAG) const
    {
        HPBC_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        namespace hc = ::hurchalla;
        bool isNegative;  // ignored
        T result = hc::REDC_incomplete(isNegative, u_hi, u_lo, n_, BC::inv_n_);
        resultIsZero = (result == 0);
        T sum = static_cast<T>(result + n_);
        HPBC_POSTCONDITION2(0 < sum && sum < static_cast<T>(2*n_));
        return V(sum);
    }
#if 1
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(T u_hi, T u_lo, PTAG) const
    {
        bool resultIsZero;  // ignored
        return montyREDC(resultIsZero, u_hi, u_lo, PTAG());
    }
#else
    HURCHALLA_FORCE_INLINE
    V montyREDC(bool& resultIsZero, T u_hi, T u_lo, LowlatencyTag) const
    {
        HPBC_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        namespace hc = ::hurchalla;
        bool isNegative;  // ignored
#if 0
// Enabling this section would result in the same code as the template version
// of this function, above.  But we can reduce latency via an optimization
// compilers don't always find, in the #else section.
        T result = hc::REDC_incomplete(isNegative, u_hi, u_lo, n_, BC::inv_n_);
        resultIsZero = (result == 0);
        result = static_cast<T>(result + n_);
#else
        u_hi = static_cast<T>(u_hi + n_);
        // REDC_incomplete uses u_hi only in a single subtract which sets the
        // result, so it makes no difference for correctness in this function if
        // we move the addition of n_ + u_hi to instead be prior to REDC.  But
        // it will lower latency to do the add before REDC.
        T result = hc::REDC_incomplete(isNegative, u_hi, u_lo, n_, BC::inv_n_);
        resultIsZero = (result == n_);
#endif
        HPBC_POSTCONDITION2(0 < result && result < static_cast<T>(2*n_));
        return V(result);
    }
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(T u_hi, T u_lo, PTAG) const
    {
        HPBC_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        namespace hc = ::hurchalla;
        bool isNegative;  // ignored
#if 0
// This is the obvious code to use, and the #else is an optimization.
// Compilers in theory should find the optimization (latest clang and gcc both
// do), but we enable the optimized version to be certain we get it.
        T result = hc::REDC_incomplete(isNegative, u_hi, u_lo, n_, BC::inv_n_);
        result = static_cast<T>(result + n_);
#else
        u_hi = static_cast<T>(u_hi + n_);
        T result = hc::REDC_incomplete(isNegative, u_hi, u_lo, n_, BC::inv_n_);
#endif
        HPBC_POSTCONDITION2(0 < result && result < static_cast<T>(2*n_));
        return V(result);
    }
#endif

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
        return (x.get() < static_cast<T>(2*n_));
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
