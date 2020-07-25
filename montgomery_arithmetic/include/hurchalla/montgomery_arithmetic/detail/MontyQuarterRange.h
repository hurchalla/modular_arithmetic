// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/montmul_quarter_range.h"
#include "hurchalla/montgomery_arithmetic/detail/monty_common.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

#if defined(_MSC_VER)
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace montgomery_arithmetic {


// Let the theoretical constant R = 2^(ma_numeric_limits<T>::digits).
template <typename T>
class MontyQuarterRange final : public MontyCommonBase<MontyQuarterRange, T> {
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    using MontyCommonBase<
      ::hurchalla::montgomery_arithmetic::MontyQuarterRange, T>::n_;
    using MontyCommonBase<
      ::hurchalla::montgomery_arithmetic::MontyQuarterRange, T>::neg_inv_n_;
    using V = typename MontyCommonBase<
      ::hurchalla::montgomery_arithmetic::MontyQuarterRange,T>::MontgomeryValue;
public:
    using montvalue_type = V;
    using template_param_type = T;

    explicit MontyQuarterRange(T modulus) : MontyCommonBase<
              ::hurchalla::montgomery_arithmetic::MontyQuarterRange, T>(modulus)
    {
        // MontyQuarterRange requires  modulus < R/4
        T Rdiv4 = static_cast<T>(static_cast<T>(1) <<
                        (modular_arithmetic::ma_numeric_limits<T>::digits - 2));
        HPBC_PRECONDITION2(modulus < Rdiv4);
    }

    MontyQuarterRange(const MontyQuarterRange&) = delete;
    MontyQuarterRange& operator=(const MontyQuarterRange&) = delete;

    static constexpr T max_modulus()
    {
        return static_cast<T>((static_cast<T>(1) <<
                   (modular_arithmetic::ma_numeric_limits<T>::digits - 2)) - 1);
    }

    HURCHALLA_FORCE_INLINE bool isValid(V x) const { return (x.get() < 2*n_); }

    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        T a = montout_non_minimized(x.get(), n_, neg_inv_n_);
        // montout_non_minimized() postconditions guarantee the following
        T minimized_result = (a >= n_) ? static_cast<T>(a - n_) : a;
        HPBC_POSTCONDITION2(minimized_result < n_);
        return minimized_result;
    }

    HURCHALLA_FORCE_INLINE V getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        T cf = (x.get() < n_) ? x.get() : static_cast<T>(x.get() - n_);
        HPBC_POSTCONDITION2(cf < n_);
        return V(cf);
    }

    HURCHALLA_FORCE_INLINE V multiply(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < 2*n_);
        // Since we have a class precondition that our modulus n_ < R/4,  we are
        // able to satisfy the preconditions to call montmul_quarter_range().

        T result = montmul_quarter_range(x.get(), y.get(), n_, neg_inv_n_);
        // montmul_quarter_range()'s postcondition guarantees the following
        HPBC_POSTCONDITION2(result < 2*n_);

        return V(result);
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < 2*n_);
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_addition_prereduced_inputs(x.get(), y.get(), 
                                                          static_cast<T>(2*n_));
        HPBC_POSTCONDITION2(z < 2*n_);
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < 2*n_);
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_subtraction_prereduced_inputs(x.get(), y.get(),
                                                          static_cast<T>(2*n_));
        HPBC_POSTCONDITION2(z < 2*n_);
        return V(z);
    }
};


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
