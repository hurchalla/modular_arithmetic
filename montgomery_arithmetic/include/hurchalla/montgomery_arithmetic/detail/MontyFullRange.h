// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


template <typename T>
class MontyFullRange : public MontyCommonBase<MontyFullRange, T> {
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    using BC = MontyCommonBase<
                         ::hurchalla::montgomery_arithmetic::MontyFullRange, T>;
    using BC::n_;
    using BC::neg_inv_n_;
    using V = typename BC::MontgomeryValue;
public:
    using montvalue_type = V;
    using template_param_type = T;
    using MontyTag = FullrangeTag;

    explicit MontyFullRange(T modulus) : BC(modulus) {}
    MontyFullRange(const MontyFullRange&) = delete;
    MontyFullRange& operator=(const MontyFullRange&) = delete;

    static constexpr T max_modulus()
    {
        return (modular_arithmetic::ma_numeric_limits<T>::max() % 2 == 0) ?
           static_cast<T>(modular_arithmetic::ma_numeric_limits<T>::max() - 1) :
           modular_arithmetic::ma_numeric_limits<T>::max();
    }

    HURCHALLA_FORCE_INLINE T getExtendedModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE V getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        return x;
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V famul(V x, V y, V z, PTAG) const
    {
        HPBC_PRECONDITION2(x.get() < n_);
        HPBC_PRECONDITION2(y.get() < n_);
        HPBC_PRECONDITION2(z.get() < n_);

        // Unfortunately for MontyFullRange, it's not possible to do a simple
        // non-modular add (x+y) prior to the multiply by z, since the sum might
        // be greater than n_.  Having a sum > n_ would break a precondition
        // for calling multiply, which requires that both input parameters are
        // less than n_.  More precisely, multiply() requires that its input
        // parameters a and b satisfy  a*b < n_*R.  If we let  a = x+y > n_  and
        // b = z, some combinations of a, z and n_ would break this requirement:
        // For example,  a=n_+3, z=n_-1, and n_=R-1  would result in 
        // a*b == (n_+3)*(n_-1) == n_*n_ + 2*n_ - 3 == n_*(R-1) + 2*n_ - 3 ==
        //        n_*R + n_ - 3 == n_*R + R-1 - 3 == n_*R + (R-4).  We can
        // safely assume R is always larger than 4, so we would have
        // a*b == n_*R + (R-4) > n_*R,  which breaks the requirement.
        // We will instead use modular addition to get the sum, which results in
        // sum < n_.  Due to the need for modular addition, famul() just wraps
        // this class's add_canonical_value and multiply functions.
        V sum = BC::add_canonical_value(x, y);
        V result = BC::multiply(sum, z, PTAG());

        HPBC_POSTCONDITION2(result.get() < n_);
        return result;
    }
};


}} // end namespace

#endif
