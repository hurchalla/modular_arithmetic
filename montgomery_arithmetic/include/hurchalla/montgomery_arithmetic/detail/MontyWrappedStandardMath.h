// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

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
template <typename T>
class MontyWrappedStandardMath final {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    T modulus_;

    struct W : public BaseMontgomeryValue<T> { // wide montgomery value type
        using BaseMontgomeryValue<T>::BaseMontgomeryValue;
    };
    struct SQ : public BaseMontgomeryValue<T> { //squaring montgomery value type
        using BaseMontgomeryValue<T>::BaseMontgomeryValue;
    };
    struct V : public W { using W::W; };  // regular montomery value type
    struct C : public V { using V::V; };  // canonical montgomery value type

    // intended for use in postconditions/preconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(W x) const
    {
        // this static_assert guarantees 0 <= x.get()
        static_assert(!(ut_numeric_limits<T>::is_signed), "");
        return (x.get() < modulus_);
    }
 public:
    using widevalue_type = W;
    using montvalue_type = V;
    using canonvalue_type = C;
    using squaringvalue_type = SQ;
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
    HURCHALLA_FORCE_INLINE T convertOut(W x) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        T ret = x.get();
        HPBC_POSTCONDITION2(0 <= ret && ret < modulus_);
        return ret;
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(W x) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        return C(x.get());
    }

    HURCHALLA_FORCE_INLINE C getUnityValue() const
    {
        HPBC_INVARIANT2(isCanonical(W(static_cast<T>(1))));
        return C(static_cast<T>(1));
    }
    HURCHALLA_FORCE_INLINE C getZeroValue() const
    {
        HPBC_INVARIANT2(isCanonical(W(static_cast<T>(0))));
        return C(static_cast<T>(0));
    }
    HURCHALLA_FORCE_INLINE C getNegativeOneValue() const
    {
        HPBC_INVARIANT2(modulus_ > 0);
        T negOne = static_cast<T>(modulus_ - static_cast<T>(1));
        HPBC_INVARIANT2(isCanonical(W(negOne)));
        return C(negOne);
    }

    template <class PTAG>   // Performance TAG (ignored by this class)
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, bool& isZero, PTAG) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result = modular_multiplication_prereduced_inputs(x.get(),
                                                             y.get(), modulus_);
        isZero = (getCanonicalValue(W(result)).get() == getZeroValue().get());
        HPBC_POSTCONDITION2(isCanonical(W(result)));
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
    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result=modular_addition_prereduced_inputs(x.get(), y.get(), modulus_);
        HPBC_POSTCONDITION2(isCanonical(W(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V add_canonical_value(V x, C y) const
    {
        V z = add(x, y);
        HPBC_POSTCONDITION2(isCanonical(z));  // add() guarantees this
        return z;
    }
    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result =
              modular_subtraction_prereduced_inputs(x.get(), y.get(), modulus_);
        HPBC_POSTCONDITION2(isCanonical(W(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V subtract_canonical_value(V x, C y) const
    {
        V z = subtract(x, y);
        HPBC_POSTCONDITION2(isCanonical(z));  // subtract() guarantees this
        return z;
    }
    HURCHALLA_FORCE_INLINE C subtract_dual_canonical_values(C x, C y) const
    {
        V z = subtract(x, y);
        HPBC_POSTCONDITION2(isCanonical(z));  // subtract() guarantees this
        return C(z.get());
    }
    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result = absolute_value_difference(x.get(), y.get());
        HPBC_POSTCONDITION2(isCanonical(W(result)));
        return V(result);
    }

    // Returns the greatest common divisor of the standard representations
    // (non-montgomery) of both x and the modulus, using the supplied functor.
    // The functor must take two integral arguments of the same type and return
    // the gcd of those two arguments.
    template <class F>
    HURCHALLA_FORCE_INLINE T gcd_with_modulus(W x, const F& gcd_functor) const
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
};


}} // end namespace

#endif
