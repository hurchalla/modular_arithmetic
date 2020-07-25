// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/negative_inverse_mod_r.h"
#include "hurchalla/modular_arithmetic/modular_multiplication.h"
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

public:
    // intended for use in postconditions/preconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        V cfx = static_cast<const D*>(this)->getCanonicalValue(x);
        // Any fully reduced value (0 <= value < n_) must be canonical.  Class
        // Derived must be implemented to respect this.
        HPBC_INVARIANT2((0 <= x.get() && x.get() < n_) ? x.get() == cfx.get() :
                                                         x.get() != cfx.get());
        bool good = static_cast<const D*>(this)->isValid(x);
        return (x.get() == cfx.get() && good);
    }

    HURCHALLA_FORCE_INLINE T getModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        return static_cast<const D*>(this)->multiply(V(a), V(r_squared_mod_n_));
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
};


}} // end namespace

#endif
