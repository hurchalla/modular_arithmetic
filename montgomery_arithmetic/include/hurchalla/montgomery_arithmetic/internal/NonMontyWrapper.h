
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_NON_MONTY_WRAPPER_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_NON_MONTY_WRAPPER_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


template <typename T>
struct NonMontyWrapper final {
private:
    static_assert(std::numeric_limits<T>::is_integer, "");
    using V = MontgomeryValue<T>;
    T modulus_;
public:
    using converting_type = T;
    using montvalue_type = V;

    explicit NonMontyWrapper(T modulus) : modulus_(modulus)
    {
        precondition(modulus > 0);
    }
    NonMontyWrapper(const NonMontyWrapper&) = delete;
    NonMontyWrapper& operator=(const NonMontyWrapper&) = delete;

    V convertIn(T a) const
    {
        precondition(a >= 0 && a < modulus_);
        return V(a);
    }
    T convertOut(V x) const
    {
        T ret = x.get();
        postcondition(ret >= 0 && ret < modulus_);
        return ret;
    }

    V getCanonicalForm(V x) const  { return x; }
    V getUnityValue() const        { return V(static_cast<T>(1)); }
    V getZeroValue() const         { return V(static_cast<T>(0)); }
    V getNegativeOneValue() const  { return V(modulus_ - static_cast<T>(1)); }

    V multiply(V x, V y) const
    {
        precondition(x.get() >= 0 && x.get() < modulus_);
        precondition(y.get() >= 0 && y.get() < modulus_);
        T result =
          modular_multiplication_prereduced_inputs(x.get(), y.get(), modulus_);
        postcondition(result >= 0 && result < modulus_);
        return V(result);
    }
    V add(V x, V y) const
    {
        precondition(x.get() >= 0 && x.get() < modulus_);
        precondition(y.get() >= 0 && y.get() < modulus_);
        T result =
          modular_addition_prereduced_inputs(x.get(), y.get(), modulus_);
        postcondition(result >= 0 && result < modulus_);
        return V(result);
    }
    V subtract(V x, V y) const
    {
        precondition(x.get() >= 0 && x.get() < modulus_);
        precondition(y.get() >= 0 && y.get() < modulus_);
        T result =
          modular_subtraction_prereduced_inputs(x.get(), y.get(), modulus_);
        postcondition(result >= 0 && result < modulus_);
        return V(result);
    }

//    bool isValid(V x) const { return (x.get() < modulus_); }
};
/*
    // Returns the modular sum of (the montgomery values) x and y.  The return
    // value is in montgomery form but might not be canonical - call
    // getCanonicalForm() to use it in comparisons.
    V add(V x, V y) const { return impl.add(x, y); }

    // Returns the modular difference of (the montgomery values) x and y.  More
    // precisely, x minus y.  The return value is in montgomery form but might
    // not be canonical - call getCanonicalForm() to use it in comparisons.
    V subtract(V x, V y) const { return impl.subtract(x, y); }
*/


}} // end namespace

#endif
