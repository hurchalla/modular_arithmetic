// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_WRAPPED_STANDARD_MATH_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
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
public:
    class MontgomeryValue {
        friend MontyWrappedStandardMath;
        template <class> friend struct MontPowImpl;
        HURCHALLA_FORCE_INLINE explicit MontgomeryValue(T val) : value(val) {}
    public:
#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Weffc++"
#endif
        // This next constructor purposely does not initialize 'value' - the
        // contents are undefined until the object is assigned to.
        HURCHALLA_FORCE_INLINE MontgomeryValue() {}
#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif
    protected:
        HURCHALLA_FORCE_INLINE T get() const { return value; }
        T value;
    };
private:
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    using V = MontgomeryValue;
    T modulus_;
public:
    using uint_type = T;
    using montvalue_type = V;

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

    // intended for use in postconditions/preconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        // this static_assert guarantees 0 <= x.get()
        static_assert(!(ut_numeric_limits<T>::is_signed), "");
        return (x.get() < modulus_);
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

    HURCHALLA_FORCE_INLINE V getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        return x;
    }

    HURCHALLA_FORCE_INLINE V getUnityValue() const
    {
        return V(static_cast<T>(1));
    }
    HURCHALLA_FORCE_INLINE V getZeroValue() const
    {
        return V(static_cast<T>(0));
    }
    HURCHALLA_FORCE_INLINE V getNegativeOneValue() const
    {
        return V(static_cast<T>(modulus_ - static_cast<T>(1)));
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
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, V z, PTAG) const
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
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, V z, PTAG) const
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
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V add_canonical_value(V x, V y) const
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
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V subtract_canonical_value(V x, V y) const
    {
        V z = subtract(x, y);
        HPBC_POSTCONDITION2(isCanonical(z));  // subtract() guarantees this
        return z;
    }
    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(x));
        HPBC_PRECONDITION2(isCanonical(y));
        T result = absolute_value_difference(x.get(), y.get());
        HPBC_POSTCONDITION2(isCanonical(V(result)));
        return V(result);
    }

    // Returns the greatest common denominator of the standard representations
    // (non-montgomery) of both x and the modulus, using the supplied functor.
    // The functor must take two integral arguments of the same type and return
    // the gcd of those two arguments.
    template <template<class> class Functor>
    HURCHALLA_FORCE_INLINE T gcd_with_modulus(MontgomeryValue x) const
    {
        // We want to return the value  q = gcd(convertOut(x), modulus_).  Since
        // this class simply wraps standard integer domain values within a
        // MontgomeryForm interface, x.get() == convertOut(x).
        Functor<T> gcd_func;
        T p = gcd_func(x.get(), modulus_);
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
