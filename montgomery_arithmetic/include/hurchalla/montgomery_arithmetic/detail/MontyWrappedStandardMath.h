// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/BaseMontgomeryValue.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/absolute_value_difference.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/programming_by_contract.h"
#include "hurchalla/util/compiler_macros.h"

namespace hurchalla { namespace detail {


// This class provides a standard modular arithmetic implementation, wrapped
// inside a Monty template.  This allows standard modular arithmetic to be used
// with a generic MontgomeryForm interface.


struct TagMontyWrappedmath final {};  //IDs MontyWrappedStandardMath without T


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

 public:
    using MontyTag = TagMontyWrappedmath;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;
    using uint_type = T;

    explicit MontyWrappedStandardMath(T modulus) : modulus_(modulus)
    {
        HPBC_PRECONDITION2(modulus > 0);
    }
    MontyWrappedStandardMath(const MontyWrappedStandardMath&) = delete;
    MontyWrappedStandardMath& operator=(const MontyWrappedStandardMath&)=delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return ut_numeric_limits<T>::max();
    }

    HURCHALLA_FORCE_INLINE T getModulus() const
    {
        return modulus_;
    }

    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        HPBC_PRECONDITION2(0 <= a);
        if HURCHALLA_LIKELY(a < modulus_)
            return V(a);
        else
            return V(static_cast<T>(a % modulus_));
    }
    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        T ret = x.get();
        HPBC_POSTCONDITION2(0 <= ret && ret < modulus_);
        return ret;
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        return C(x.get());
    }

    HURCHALLA_FORCE_INLINE C getUnityValue() const
    {
        HPBC_INVARIANT2(isCanonical(V(static_cast<T>(1))));
        return C(static_cast<T>(1));
    }
    HURCHALLA_FORCE_INLINE C getZeroValue() const
    {
        HPBC_INVARIANT2(isCanonical(V(static_cast<T>(0))));
        return C(static_cast<T>(0));
    }
    HURCHALLA_FORCE_INLINE C getNegativeOneValue() const
    {
        HPBC_INVARIANT2(modulus_ > 0);
        T negOne = static_cast<T>(modulus_ - static_cast<T>(1));
        HPBC_INVARIANT2(isCanonical(V(negOne)));
        return C(negOne);
    }

    HURCHALLA_FORCE_INLINE V negate(V x) const
    {
        return subtract(getZeroValue(), x);
    }

    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, bool& isZero, PTAG) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result = modular_multiplication_prereduced_inputs(x.get(),
                                                             y.get(), modulus_);
        isZero = (getCanonicalValue(V(result)).get() == getZeroValue().get());
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, C z, PTAG) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        HPBC_PRECONDITION2(isCanonical(z));
        bool isZero;
        V product = multiply(x, y, isZero, PTAG());
        V result = subtract(product, z);
        HPBC_POSTCONDITION2(isCanonical(result));
        return result;
    }
    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, C z, PTAG) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        HPBC_PRECONDITION2(isCanonical(z));
        bool isZero;
        V product = multiply(x, y, isZero, PTAG());
        V result = add(product, z);
        HPBC_POSTCONDITION2(isCanonical(result));
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
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result=modular_addition_prereduced_inputs(x.get(), y.get(), modulus_);
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    // Note: add(V, C) and add(C, V) will match to add(V x, V y) above.
    HURCHALLA_FORCE_INLINE C add(C x, C y) const
    {
        V v = add(V(x), V(y));
        return C(v.get());
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result =
              modular_subtraction_prereduced_inputs(x.get(), y.get(), modulus_);
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    // Note: subtract(V, C) and subtract(C, V) will match to subtract(V x, V y)
    // above.
    HURCHALLA_FORCE_INLINE C subtract(C x, C y) const
    {
        V v = subtract(V(x), V(y));
        return C(v.get());
    }

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result = absolute_value_difference(x.get(), y.get());
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    // Note: unordered_subtract(V, C) and unordered_subtract(C, V) will match
    // to unordered_subtract(V x, V y) above.


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
        HPBC_INVARIANT2(modulus_ > 0);
        // We want to return the value  q = gcd(convertOut(x), modulus_).  Since
        // this class simply wraps standard integer domain values within a
        // MontgomeryForm interface, x.get() == convertOut(x).
        T p = gcd_functor(x.get(), modulus_);
        // Our postconditions assume the Functor implementation is correct.
        HPBC_POSTCONDITION2(0 < p && p <= modulus_ &&
                            (x.get() == 0 || p <= x.get()));
        HPBC_POSTCONDITION2(modulus_ % p == 0);
        HPBC_POSTCONDITION2(x.get() % p == 0);
        return p;
    }

    // This class doesn't do anything special for square functions.
    // It just delegates to the functions above.
    // --------
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V square(V x, PTAG) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        bool isZero;
        return multiply(x, x, isZero, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fusedSquareSub(V x, C cv, PTAG) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        return fmsub(x, x, cv, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fusedSquareAdd(V x, C cv, PTAG) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        return fmadd(x, x, cv, PTAG());
    }
};


}} // end namespace

#endif
