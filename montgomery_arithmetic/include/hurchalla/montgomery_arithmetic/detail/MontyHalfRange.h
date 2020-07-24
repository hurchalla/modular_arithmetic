// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_HALF_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_HALF_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/montmul_half_range.h"
#include "hurchalla/montgomery_arithmetic/detail/monty_common.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


// MontyHalfRange is exactly the same as MontyFullRange, except that the
// constructor has the precondition that modulus < R/2, and in that multiply()
// calls montmul_half_range() rather than montmul_full_range().
// [The theoretical constant R = 2^(ma_numeric_limits<T>::digits).]
template <typename T>
class MontyHalfRange final : public MontyCommonBase<MontyHalfRange, T> {
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    using MontyCommonBase<
        ::hurchalla::montgomery_arithmetic::MontyHalfRange, T>::n_;
    using MontyCommonBase<
        ::hurchalla::montgomery_arithmetic::MontyHalfRange, T>::neg_inv_n_;
    using typename MontyCommonBase<
        ::hurchalla::montgomery_arithmetic::MontyHalfRange, T>::CanonicalValue;
    using V = typename MontyCommonBase<
        ::hurchalla::montgomery_arithmetic::MontyHalfRange, T>::MontgomeryValue;
public:
    using montvalue_type = V;
    using template_param_type = T;
    using canonical_value_type = CanonicalValue;

    explicit MontyHalfRange(T modulus) : MontyCommonBase<
        ::hurchalla::montgomery_arithmetic::MontyHalfRange, T>(modulus)
    {
        // MontyHalfRange requires  modulus < R/2
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(modulus < Rdiv2);
    }

    MontyHalfRange(const MontyHalfRange&) = delete;
    MontyHalfRange& operator=(const MontyHalfRange&) = delete;

    static constexpr T max_modulus()
    {
        return static_cast<T>((static_cast<T>(1) <<
                   (modular_arithmetic::ma_numeric_limits<T>::digits - 1)) - 1);
    }

    HURCHALLA_FORCE_INLINE bool isValid(V x) const { return (x.get() < n_); }

    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        T a = montout_non_minimized(x.get(), n_, neg_inv_n_);
        // montout_non_minimized() postconditions guarantee that since x < n_,
        // a < n_.  Thus 'a' is already minimized.
        T minimized_result = a;
        HPBC_POSTCONDITION2(minimized_result < n_);
        return minimized_result;
    }

    HURCHALLA_FORCE_INLINE CanonicalValue getCanonicalForm(V x) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        return CanonicalValue(x.get());
    }

    // Aside from the constructor's precondition of modulus < R/2, this
    // function is the only thing that differs from class MontyFullRange.
    HURCHALLA_FORCE_INLINE V multiply(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        HPBC_PRECONDITION2(y.get() < n_);

        T result = montmul_half_range(x.get(), y.get(), n_, neg_inv_n_);
        // montmul_half_range()'s postcondition guarantees the following
        HPBC_POSTCONDITION2(result < n_);

        return V(result);
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        HPBC_PRECONDITION2(y.get() < n_);
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_addition_prereduced_inputs(x.get(), y.get(), n_);
        HPBC_POSTCONDITION2(z < n_);
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        HPBC_PRECONDITION2(y.get() < n_);
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_subtraction_prereduced_inputs(x.get(), y.get(), n_);
        HPBC_POSTCONDITION2(z < n_);
        return V(z);
    }
};


}} // end namespace

#endif
