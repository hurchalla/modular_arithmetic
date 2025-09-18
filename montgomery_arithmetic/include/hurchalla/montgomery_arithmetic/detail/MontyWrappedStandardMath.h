// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/BaseMontgomeryValue.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/get_R_mod_n.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyTags.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_multiplicative_inverse.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/absolute_value_difference.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// This class provides a standard modular arithmetic implementation, wrapped
// inside a Monty template.  This allows standard modular arithmetic to be used
// with a generic MontgomeryForm interface.


template <typename T>
class MontyWrappedStandardMath final {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    T modulus_;

    struct V : public BaseMontgomeryValue<T> {  // regular montgomery value type
        HURCHALLA_FORCE_INLINE V() = default;
     protected:
        friend MontyWrappedStandardMath;
        HURCHALLA_FORCE_INLINE explicit V(T a) : BaseMontgomeryValue<T>(a) {}
    };
    struct C : public V {                     // canonical montgomery value type
        HURCHALLA_FORCE_INLINE C() = default;
        HURCHALLA_FORCE_INLINE friend bool operator==(const C& x, const C& y)
            { return x.get() == y.get(); }
        HURCHALLA_FORCE_INLINE friend bool operator!=(const C& x, const C& y)
            { return !(x == y); }
     protected:
        friend MontyWrappedStandardMath;
        HURCHALLA_FORCE_INLINE explicit C(T a) : V(a) {}
    };
    struct FV : public V {                     // fusing montgomery value type
        HURCHALLA_FORCE_INLINE FV() = default;
     protected:
        friend MontyWrappedStandardMath;
        HURCHALLA_FORCE_INLINE explicit FV(T a) : V(a) {}
    };

    // intended for use in postconditions/preconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        // this static_assert guarantees 0 <= x.get()
        static_assert(!(ut_numeric_limits<T>::is_signed), "");
        return (x.get() < modulus_);
    }

    using SV = V;

 public:
    using MontyTag = TagMontyWrappedmath;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;
    using squaringvalue_type = SV;
    using uint_type = T;

    explicit MontyWrappedStandardMath(T modulus) : modulus_(modulus)
    {
        HPBC_CLOCKWORK_PRECONDITION2(modulus > 0);
    }

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return ut_numeric_limits<T>::max();
    }

    HURCHALLA_FORCE_INLINE T getModulus() const
    {
        return modulus_;
    }

    template <class PTAG>   // PTAG is ignored by this class
    HURCHALLA_FORCE_INLINE V convertIn(T a, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(0 <= a);
        if HURCHALLA_LIKELY(a < modulus_)
            return V(a);
        else
            return V(static_cast<T>(a % modulus_));
    }
    template <class PTAG>   // PTAG is ignored by this class
    HURCHALLA_FORCE_INLINE T convertOut(V x, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        T ret = x.get();
        HPBC_CLOCKWORK_POSTCONDITION2(0 <= ret && ret < modulus_);
        return ret;
    }

    template <class PTAG>   // PTAG is ignored by this class
    HURCHALLA_FORCE_INLINE T remainder(T a, PTAG) const
    {
        return static_cast<T>(a % modulus_);
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        return C(x.get());
    }

    HURCHALLA_FORCE_INLINE C getUnityValue() const
    {
        HPBC_CLOCKWORK_INVARIANT2(isCanonical(V(static_cast<T>(1))));
        return C(static_cast<T>(1));
    }
    HURCHALLA_FORCE_INLINE C getZeroValue() const
    {
        HPBC_CLOCKWORK_INVARIANT2(isCanonical(V(static_cast<T>(0))));
        return C(static_cast<T>(0));
    }
    HURCHALLA_FORCE_INLINE C getNegativeOneValue() const
    {
        HPBC_CLOCKWORK_INVARIANT2(modulus_ > 0);
        T negOne = static_cast<T>(modulus_ - static_cast<T>(1));
        HPBC_CLOCKWORK_INVARIANT2(isCanonical(V(negOne)));
        return C(negOne);
    }

    HURCHALLA_FORCE_INLINE V negate(V x) const
    {
        return subtract(getZeroValue(), x, 0);  // 0 is arbitrary, for PTAG
    }

    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, bool& isZero, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(y));
        T result = ::hurchalla::modular_multiplication_prereduced_inputs(
                                                    x.get(), y.get(), modulus_);
        isZero = (getCanonicalValue(V(result)).get() == getZeroValue().get());
        HPBC_CLOCKWORK_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, C z, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(y));
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(z));
        bool isZero;
        V product = multiply(x, y, isZero, PTAG());
        V result = subtract(product, z, 0);  // 0 is arbitrary, for PTAG
        HPBC_CLOCKWORK_POSTCONDITION2(isCanonical(result));
        return result;
    }
    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, C z, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(y));
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(z));
        bool isZero;
        V product = multiply(x, y, isZero, PTAG());
        V result = add(product, z);
        HPBC_CLOCKWORK_POSTCONDITION2(isCanonical(result));
        return result;
    }

    // Note: internal to MontyWrappedStandardMath the contents of FusingValue
    // (FV) and CanonicalValue (C) variables are interchangeable.  Other Monty
    // types use FV and C as completely distinct types, and so for genericity we
    // always present C and FV to the outside world as being unrelated.
    HURCHALLA_FORCE_INLINE FV getFusingValue(V x) const
    {
        C cv = getCanonicalValue(x);
        return FV(cv.get());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, FV fv, PTAG) const
    {
        C cv = C(fv.get());
        return fmadd(x, y, cv, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, FV fv, PTAG) const
    {
        C cv = C(fv.get());
        return fmsub(x, y, cv, PTAG());
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(y));
        T result = ::hurchalla::modular_addition_prereduced_inputs(
                                                    x.get(), y.get(), modulus_);
        HPBC_CLOCKWORK_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    // Note: add(V, C) and add(C, V) will match to add(V x, V y) above.
    HURCHALLA_FORCE_INLINE C add(C x, C y) const
    {
        V v = add(V(x), V(y));
        return C(v.get());
    }

    template <class PTAG>
    HURCHALLA_FORCE_INLINE V subtract(V x, V y, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(y));
        T result = ::hurchalla::modular_subtraction_prereduced_inputs(
                                                    x.get(), y.get(), modulus_);
        HPBC_CLOCKWORK_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    // Note: subtract(V, C, PTAG) and subtract(C, V, PTAG) will match to
    // subtract(V x, V y, PTAG) above.
    template <class PTAG>
    HURCHALLA_FORCE_INLINE C subtract(C x, C y, PTAG) const
    {
        V v = subtract(V(x), V(y), 0);  // 0 is arbitrary, for PTAG
        return C(v.get());
    }

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(y));
        T result = ::hurchalla::absolute_value_difference(x.get(), y.get());
        HPBC_CLOCKWORK_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    // Note: unordered_subtract(V, C) and unordered_subtract(C, V) will match
    // to unordered_subtract(V x, V y) above.

    HURCHALLA_FORCE_INLINE V two_times(V x) const
    {
        return add(x, x);
    }
    HURCHALLA_FORCE_INLINE C two_times(C cx) const
    {
        return add(cx, cx);
    }


    HURCHALLA_FORCE_INLINE SV getSquaringValue(V x) const
    {
        static_assert(std::is_same<V, SV>::value, "");
        return x;
    }
    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE SV squareSV(SV sv, PTAG) const
    {
        static_assert(std::is_same<V, SV>::value, "");
        return square(sv, PTAG());
    }
    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V squareToMontgomeryValue(SV sv, PTAG) const
    {
        static_assert(std::is_same<V, SV>::value, "");
        return square(sv, PTAG());
    }
    HURCHALLA_FORCE_INLINE V getMontgomeryValue(SV sv) const
    {
        static_assert(std::is_same<V, SV>::value, "");
        return sv;
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE C inverse(V x, PTAG) const
    {
        namespace hc = ::hurchalla;
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        T gcd;  // ignored
        T inv = hc::modular_multiplicative_inverse(x.get(), modulus_, gcd);

        HPBC_CLOCKWORK_POSTCONDITION2(inv < modulus_);
        //POSTCONDITION: Return 0 if the inverse does not exist. Otherwise
        //   return the value of the inverse (which would never be 0, given that
        //   modulus_ > 1).
        HPBC_CLOCKWORK_POSTCONDITION2(inv == 0 || 1 ==
            hc::modular_multiplication_prereduced_inputs(inv,x.get(),modulus_));
        return C(inv);
    }

    template <class PTAG> HURCHALLA_FORCE_INLINE
    V divideBySmallPowerOf2(C cx, int exponent, PTAG) const
    {
        V pow_of_two = twoPowLimited(static_cast<size_t>(exponent), PTAG());
        C inv_pow_of_two = inverse(pow_of_two, PTAG());
        C zero = getZeroValue();
        HPBC_CLOCKWORK_ASSERT2(inv_pow_of_two != zero);
        bool isZero;
        V product = multiply(inv_pow_of_two, cx, isZero, PTAG());
        HPBC_CLOCKWORK_ASSERT2((cx == zero) == isZero);
        C result = getCanonicalValue(product);
        return result;
    }

    // Returns the greatest common divisor of the standard representations
    // (non-montgomery) of both x and the modulus, using the supplied functor.
    // The functor must take two integral arguments of the same type and return
    // the gcd of those two arguments.  Usually you would make the functor's
    // operator() a templated function, where the template parameter is the
    // unknown type of the integral arguments.  Or more simply, you can just use
    // a lambda, with 'auto' type for the function parameters.
    template <class F>
    HURCHALLA_FORCE_INLINE T gcd_with_modulus(V x, const F& gcd_functor) const
    {
        HPBC_CLOCKWORK_INVARIANT2(modulus_ > 0);
        // We want to return the value  q = gcd(convertOut(x), modulus_).  Since
        // this class simply wraps standard integer domain values within a
        // MontgomeryForm interface, x.get() == convertOut(x).
        T p = gcd_functor(x.get(), modulus_);
        // Our postconditions assume the Functor implementation is correct.
        HPBC_CLOCKWORK_POSTCONDITION2(0 < p && p <= modulus_ &&
                            (x.get() == 0 || p <= x.get()));
        HPBC_CLOCKWORK_POSTCONDITION2(modulus_ % p == 0);
        HPBC_CLOCKWORK_POSTCONDITION2(x.get() % p == 0);
        return p;
    }

    // This class doesn't do anything special for square functions.
    // It just delegates to the functions above.
    // --------
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V square(V x, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        bool isZero;
        return multiply(x, x, isZero, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fusedSquareSub(V x, C cv, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        return fmsub(x, x, cv, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fusedSquareAdd(V x, C cv, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isCanonical(x));
        return fmadd(x, x, cv, PTAG());
    }


    // returns R mod N
    HURCHALLA_FORCE_INLINE C getMontvalueR() const
    {
        T result = ::hurchalla::get_R_mod_n(modulus_);
        HPBC_CLOCKWORK_POSTCONDITION2(result < modulus_);
        return C(result);
    }
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V twoPowLimited_times_x(size_t exponent, C cx, PTAG) const
    {
        static constexpr int digitsT = ut_numeric_limits<T>::digits;
        int power = static_cast<int>(exponent);
        HPBC_CLOCKWORK_PRECONDITION2(0 <= power && power < digitsT);

        T tmp = cx.get();
        HPBC_CLOCKWORK_INVARIANT2(tmp < modulus_);
        T u_lo = static_cast<T>(tmp << power);
        int rshift = digitsT - power;
        HPBC_CLOCKWORK_ASSERT2(rshift > 0);
        T u_hi = static_cast<T>(tmp >> 1) >> (rshift - 1);
        HPBC_CLOCKWORK_ASSERT2(u_hi < modulus_);
        // It's very strange to use REDC when this class is meant to wrap
        // standard arithmetic within the monty interface and not actually
        // use mont arith.  But we need REDC here, due to the extra R factor
        // that is expected to be in cx whenever this function is called.
        T inv_modulus = ::hurchalla::inverse_mod_R(modulus_);
        T result = ::hurchalla::REDC_standard(u_hi, u_lo, modulus_, inv_modulus, PTAG());

        HPBC_CLOCKWORK_POSTCONDITION2(result < modulus_);
        return V(result);
    }
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V twoPowLimited_times_x_v2(size_t exponent, C cx, PTAG) const
    {
        static constexpr int digitsT = ut_numeric_limits<T>::digits;
        int power = static_cast<int>(exponent);
        HPBC_CLOCKWORK_PRECONDITION2(0 < power && power <= digitsT);

        T tmp = cx.get();
        HPBC_CLOCKWORK_INVARIANT2(tmp < modulus_);
        T u_lo = static_cast<T>(static_cast<T>(tmp << 1) << (power - 1));
        int rshift = digitsT - power;
        HPBC_CLOCKWORK_ASSERT2(0 <= rshift && rshift < digitsT);
        T u_hi = static_cast<T>(tmp >> rshift);

        HPBC_CLOCKWORK_ASSERT2(u_hi < modulus_);
        // It's very strange to use REDC when this class is meant to wrap
        // standard arithmetic within the monty interface and not actually
        // use mont arith.  But we need REDC here, due to the extra R factor
        // that is expected to be in cx whenever this function is called.
        T inv_modulus = ::hurchalla::inverse_mod_R(modulus_);
        T result = ::hurchalla::REDC_standard(u_hi, u_lo, modulus_, inv_modulus, PTAG());

        HPBC_CLOCKWORK_POSTCONDITION2(result < modulus_);
        return V(result);
    }

    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE T getMagicValue(PTAG) const
    {
        T result = ::hurchalla::get_R_mod_n(modulus_);
        HPBC_CLOCKWORK_POSTCONDITION2(result < modulus_);
        return result;
    }
    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V convertInExtended_aTimesR(T a, T RmodN, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(RmodN == getMagicValue(PTAG()));
        T tmp = a;
        if (tmp >= modulus_)
            tmp = static_cast<T>(tmp % modulus_);
        HPBC_CLOCKWORK_ASSERT2(tmp < modulus_);
        HPBC_CLOCKWORK_ASSERT2(RmodN < modulus_);
        T result = ::hurchalla::modular_multiplication_prereduced_inputs(
                                                          tmp, RmodN, modulus_);
        HPBC_CLOCKWORK_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }

    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V twoPowLimited(size_t exponent, PTAG) const
    {
        static constexpr int digitsT = ut_numeric_limits<T>::digits;
        int power = static_cast<int>(exponent);
        HPBC_CLOCKWORK_PRECONDITION2(0 <= power && power < digitsT);
        T tmp = static_cast<T>(static_cast<T>(1) << power);
        if (tmp >= modulus_)
            tmp = static_cast<T>(tmp % modulus_);
        HPBC_CLOCKWORK_POSTCONDITION2(tmp < modulus_);
        return V(tmp);
    }
    // PTAG Performance TAG - ignored by this class
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V RTimesTwoPowLimited(size_t exponent, T RmodN, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(RmodN == getMagicValue(PTAG()));
        static constexpr int digitsT = ut_numeric_limits<T>::digits;
        int power = static_cast<int>(exponent);
        HPBC_CLOCKWORK_PRECONDITION2(0 <= power && power < digitsT);

        T tmp = static_cast<T>(static_cast<T>(1) << power);
        if (tmp >= modulus_)
            tmp = static_cast<T>(tmp % modulus_);

        HPBC_CLOCKWORK_ASSERT2(tmp < modulus_);
        HPBC_CLOCKWORK_ASSERT2(RmodN < modulus_);
        T result = ::hurchalla::modular_multiplication_prereduced_inputs(
                                                          tmp, RmodN, modulus_);
        HPBC_CLOCKWORK_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
};


}} // end namespace

#endif
