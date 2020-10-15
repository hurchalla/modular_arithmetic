// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_QUARTER_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/monty_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


// Let the theoretical constant R = 2^(ma_numeric_limits<T>::digits).
template <typename T>
class MontyQuarterRange : public MontyCommonBase<MontyQuarterRange, T> {
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    using BC = MontyCommonBase<
                      ::hurchalla::montgomery_arithmetic::MontyQuarterRange, T>;
    using BC::n_;
    using V = typename BC::MontgomeryValue;
public:
    using montvalue_type = V;
    using template_param_type = T;
    using MontyTag = QuarterrangeTag;

    explicit MontyQuarterRange(T modulus) : BC(modulus)
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

    HURCHALLA_FORCE_INLINE T getExtendedModulus() const
    {
        return static_cast<T>(2*n_);
    }

    HURCHALLA_FORCE_INLINE V getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        T cv = (x.get() < n_) ? x.get() : static_cast<T>(x.get() - n_);
        HPBC_POSTCONDITION2(cv < n_);
        return V(cv);
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V famul(V x, V y, V z, PTAG) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < n_);   // y must be canonical
        HPBC_PRECONDITION2(z.get() < 2*n_);

        // Unfortunately for MontyQuarterRange, it's not possible to do a simple
        // non-modular add (x+y) prior to the multiply by z, since the sum might
        // be greater than 2*n_.  Having a sum > 2*n_ would break a precondition
        // for calling multiply, which requires  a*b < n_*R.  To see why, let
        // a = sum = x+y, and let b = z.  This gives us
        // a*b == (x+y)*z < (2*n_ + n_)*2*n_ == 6*n_*n_.  Our class
        // constructor requires modulus < R/4, and so we know n_ < R/4.  Thus we
        // would have  a*b == 6*n_*n_ < 6*n_*(R/4) == (3/2)*n_*R.  This result
        // of  a*b < (3/2)*n_*R  is insufficient to guarantee that we would
        // satisfy the multiply requirement of  a*b < n_*R.
        // Instead we will use modular addition to get the sum, which results in
        // sum < 2*n_ (this guarantees we satisfy  a*b < n_*R).
        // Due to the need for modular addition, famul() just wraps this class's
        // add_canonical_value and multiply functions.
        V sum = BC::add_canonical_value(x, y);
        V result = BC::multiply(sum, z, PTAG());

        // multiply's postcondition guarantees the following for this monty type
        HPBC_POSTCONDITION2(0 < result.get() && result.get() < 2*n_);
        return result;
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V famulIsZero(V x, V y, V z, bool& isZero,PTAG) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < n_);   // y must be canonical
        HPBC_PRECONDITION2(z.get() < 2*n_);

        V result = famul(x, y, z, PTAG());
        // The canonical zeroValue == (0*R)%n_ == 0.  The equivalence class for
        // zeroValue therefore is composed of values that satisfy  0 + m*n_,
        // where m is any integer.  Since famul()'s postcondition guarantees
        // 0 < result < 2*n, the only result that can belong to zeroValue's
        // equivalence class is result == n_.
        // Note: other functions like add() or sub() have return results in the
        // range of 0 <= result < 2*n, and so our postcondition of 0 < result is
        // particular to multiply().  It is not a general invariant property of
        // montgomery values for this monty type.
        isZero = (result.get() == n_);

        HPBC_POSTCONDITION2(0 < result.get() && result.get() < 2*n_);
        return result;
    }
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V multiplyIsZero(V x, V y, bool& isZero, PTAG) const
    {
        HPBC_PRECONDITION2(x.get() < 2*n_);
        HPBC_PRECONDITION2(y.get() < 2*n_);

        V result = BC::multiply(x, y, PTAG());
        // multiply()'s postcondition guarantees 0 < result < 2*n.  Thus as
        // shown in famulIsZero()'s comments, the only result value that can
        // belong to the zero value equivalence class is result == n_.
        isZero = (result.get() == n_);

        HPBC_POSTCONDITION2(0 < result.get() && result.get() < 2*n_);
        return result;
    }
};


}} // end namespace

#endif
