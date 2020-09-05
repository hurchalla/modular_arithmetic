// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/detail/monty_common.h"
#include "hurchalla/montgomery_arithmetic/detail/negative_inverse_mod_r.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/MontHelper.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/modular_arithmetic/absolute_value_difference.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"

namespace hurchalla { namespace montgomery_arithmetic {


// For discussion purposes throughout this file, given an unsigned integral type
// T, let R = 2^(ma_numeric_limits<T>::digits).  For example: if T is uint64_t
// then R = 2^64.  The name 'R' is based on the wikipedia presentation.
//
// This base class uses the CRTP idiom
// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
// This is the base class shared by most montgomery forms (MontySqrtRange is an
// exception).
template <template <typename> class Derived, typename T>
class MontyCommonBase {
public:
    class MontgomeryValue {
        friend Derived<T>; friend MontyCommonBase;
        explicit MontgomeryValue(T val) : value(val) {}
    public:
#ifdef __GNUC__
#  pragma GCC diagnostic push
#  pragma GCC diagnostic ignored "-Weffc++"
#endif
        MontgomeryValue() {} // This constructor purposely does not initialize
        // 'value' - the contents are undefined until the object is assigned to.
#ifdef __GNUC__
#  pragma GCC diagnostic pop
#endif
    protected:
        T get() const { return value; }
        T value;
    };
private:
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
    static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
    using D = Derived<T>;
    using V = MontgomeryValue;
protected:
    const T n_;   // the modulus
    const T neg_inv_n_;
    const T r_mod_n_;
    const T r_squared_mod_n_;

    explicit MontyCommonBase(T modulus) : n_(modulus),
                              neg_inv_n_(negative_inverse_mod_r(n_)),
                              r_mod_n_(getRModN(n_)),
                              r_squared_mod_n_( modular_arithmetic::
                                  modular_multiplication_prereduced_inputs(
                                                       r_mod_n_, r_mod_n_, n_) )
    {
        HPBC_PRECONDITION2(modulus % 2 == 1);
        HPBC_PRECONDITION2(modulus > 1);
        // Note: unityValue == (the montgomery form of 1)==(1*R)%n_ == r_mod_n_.
        //
        // getRModN() guarantees the below.  getUnityValue() and
        // getNegativeOneValue() both rely on it.
        HPBC_INVARIANT2(0 < r_mod_n_ && r_mod_n_ < n_);
        HPBC_INVARIANT2(r_squared_mod_n_ < n_);
    }
    MontyCommonBase(const MontyCommonBase&) = delete;
    MontyCommonBase& operator=(const MontyCommonBase&) = delete;

private:
    static T getRModN(T n)
    {
        HPBC_PRECONDITION2(n % 2 == 1);
        HPBC_PRECONDITION2(n > 1);
        // Assign a tmp T variable rather than directly using the intermediate
        // expression, in order to avoid a negative value (and a wrong answer)
        // in cases where 'n' would be promoted to type 'int'.
        T tmp = static_cast<T>(static_cast<T>(0) - n);
        // Compute R%n.  For example, if R==2^64, arithmetic wraparound behavior
        // of the unsigned integral type T results in (0 - n) representing
        // (2^64 - n).  Thus, rModN = R%n == (2^64)%n == (2^64 - n)%n == (0-n)%n
        T rModN = static_cast<T>(tmp % n);
        // Since n is odd and > 1, n does not divide R==2^x.  Thus, rModN != 0
        HPBC_POSTCONDITION2(0 < rModN && rModN < n);
        return rModN;
    }

protected:
    // intended for use in preconditions/postconditions
    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        T em = static_cast<const D*>(this)->getExtendedModulus();
        return (x.get() < em);
    }

public:
    // intended for use in preconditions/postconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        V cfx = static_cast<const D*>(this)->getCanonicalValue(x);
        // Any fully reduced value (0 <= value < n_) must be canonical.  Class
        // Derived must be implemented to respect this.
        HPBC_INVARIANT2((0 <= x.get() && x.get() < n_) ? x.get() == cfx.get() :
                                                         x.get() != cfx.get());
        bool good = isValid(x);
        return (x.get() == cfx.get() && good);
    }

    HURCHALLA_FORCE_INLINE T getModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        return multiply(V(a), V(r_squared_mod_n_), OutofplaceLowlatencyTag());
    }

    HURCHALLA_FORCE_INLINE T convertOut(V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        T result = montout(x.get(), n_, neg_inv_n_, typename D::MontyTag());
        // montout() postcondition guarantees result < n_
        HPBC_POSTCONDITION2(result < n_);
        return result;
    }

    HURCHALLA_FORCE_INLINE V getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_
        HPBC_POSTCONDITION2(isCanonical(V(r_mod_n_)));
        return V(r_mod_n_);
    }

    HURCHALLA_FORCE_INLINE V getZeroValue() const
    {
        // zeroValue == (0*R)%n_
        HPBC_POSTCONDITION2(isCanonical(V(0)));
        return V(0);
    }

    HURCHALLA_FORCE_INLINE V getNegativeOneValue() const
    {
        // We want to get returnVal = getCanonicalValue(subtract(getZeroValue(),
        //                                               getUnityValue())).
        //   getZeroValue() returns a value belonging to the equivalence class
        //   0*R (mod n_).  This equivalence class can equally be represented by
        //   the value n_ (mod n_).  getUnityValue() returns a value belonging
        //   to the equivalence class 1*R (mod n_), which is the same as 
        //   r_mod_n_ (mod n_).  Therefore the subtraction results in the
        //   equivalence class (n_ - r_mod_n_) (mod n_).
        //   The constructor established the invariant  0 < r_mod_n_ < n_
        //   Thus we also know  0 < n_ - r_mod_n_ < n_.  This means
        //   (n_ - r_mod_n_)  is fully reduced, and thus canonical.
        HPBC_INVARIANT2(n_ > r_mod_n_);
        T ret = static_cast<T>(n_ - r_mod_n_);
        HPBC_ASSERT2(0 < ret && ret < n_);

        HPBC_POSTCONDITION2(isCanonical(V(ret)));
        return V(ret);
    }

    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T em = static_cast<const D*>(this)->getExtendedModulus();
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_addition_prereduced_inputs(x.get(), y.get(), em);
        HPBC_PRECONDITION2(isValid(V(z)));
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V add_canonical_value(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isCanonical(y));
        HPBC_PRECONDITION2(y.get() < n_); // isCanonical() should guarantee this
        T z = MontHelper<T>::modadd_canonical_second_addend(x.get(), y.get(),
                                                                            n_);
        HPBC_PRECONDITION2(isValid(V(z)));
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        T em = static_cast<const D*>(this)->getExtendedModulus();
        namespace ma = hurchalla::modular_arithmetic;
        T z = ma::modular_subtraction_prereduced_inputs(x.get(), y.get(), em);
        HPBC_PRECONDITION2(isValid(V(z)));
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V subtract_canonical_value(V x, V y) const
    {
        HPBC_PRECONDITION2(isCanonical(y));
        HPBC_PRECONDITION2(y.get() < n_); // isCanonical() should guarantee this
        HPBC_PRECONDITION2(isValid(x));
        T z = MontHelper<T>::modsub_canonical_subtrahend(x.get(), y.get(), n_);
        HPBC_POSTCONDITION2(isValid(V(z)));
        return V(z);
    }

    HURCHALLA_FORCE_INLINE V unordered_subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        namespace ma = hurchalla::modular_arithmetic;
        T result = ma::absolute_value_difference(x.get(), y.get());
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V multiply(V x, V y, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        // As a precondition, montmul requires  x*y < n*R.  This will always be
        // satisfied for all Monty classes known to derive from this class --
        // For MontyFullRange: its constructor requires modulus < R, which
        //   means n < R.  Since MontyFullRange's isValid(a) returns (a < n),
        //   we know by this function's preconditions that x < n and y < n.
        //   Therefore  x*y < n*n < n*R.
        // For MontyHalfRange: its constructor requires modulus < R/2, which
        //   means n < R/2.  Its isValid(a) returns (a < n), so we know by this
        //   function's preconditions that x < n and y < n.  Thus
        //   x*y < n*n < n*R/2 < n*R.
        // For MontyQuarterRange: its constructor requires modulus < R/4, which
        //   means n < R/4.  Its isValid(a) returns (a < 2*n), so we know by
        //   this function's preconditions that x < 2*n and y < 2*n.  Thus
        //   x*y < (2*n)*(2*n) == 4*n*n < 4*n*R/4 == n*R.
        // For MontySixthRange: its constructor requires modulus < R/6, which
        //   means n < R/6.  Its isValid(a) returns (a < 2*n), so we know by
        //   this function's preconditions that x < 2*n and y < 2*n.  Thus
        //   x*y < (2*n)*(2*n) == 4*n*n < 4*n*R/6 == (2/3)*n*R < n*R.
        T result = montmul(x.get(), y.get(), n_, neg_inv_n_,
                                                typename D::MontyTag(), PTAG());
        HPBC_PRECONDITION2(isValid(V(result)));
        return V(result);
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, V z, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        HPBC_PRECONDITION2(isCanonical(z));
        HPBC_PRECONDITION2(z.get() < n_); // isCanonical() should guarantee this
        // As a precondition, montfmsub requires  x*y < n*R.  See multiply() for
        // why this will always be satisfied.
        T result = montfmsub(x.get(), y.get(), z.get(), n_, neg_inv_n_,
                                                typename D::MontyTag(), PTAG());
        HPBC_PRECONDITION2(isValid(V(result)));
        return V(result);
    }

    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, V z, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        HPBC_PRECONDITION2(isCanonical(z));
        HPBC_PRECONDITION2(z.get() < n_); // isCanonical() should guarantee this
        // As a precondition, montfmadd requires  x*y < n*R.  See multiply() for
        // why this will always be satisfied.
        T result = montfmadd(x.get(), y.get(), z.get(), n_, neg_inv_n_,
                                                typename D::MontyTag(), PTAG());
        HPBC_PRECONDITION2(isValid(V(result)));
        return V(result);
    }
};


}} // end namespace

#endif
