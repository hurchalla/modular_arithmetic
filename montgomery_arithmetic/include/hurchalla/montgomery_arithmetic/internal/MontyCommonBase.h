
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/internal/negative_inverse_mod_r.h"
#include "hurchalla/montgomery_arithmetic/internal/monty_common.h"
#include "hurchalla/montgomery_arithmetic/internal/MontgomeryValue.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "hurchalla/montgomery_arithmetic/internal/compiler_macros.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// For discussion purposes throughout this file, given an unsigned integral type
// T, let R = 2^(std::numeric_limits<T>::digits).  For example: if T is uint64_t
// then R = 2^64.  The name 'R' is based on the wikipedia presentation.
//
// This base class uses the CRTP idiom
// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
// This is the base class shared by most montgomery forms (MontySqrtRange is an
//exception).
template <template <typename> class Derived, typename T>
class MontyCommonBase {
    using D = Derived<T>;
public:
    using V = MontgomeryValue<T>;
protected:
    const T n_;   // the modulus
    const T neg_inv_n_;
    const T r_mod_n_;
    const T r_squared_mod_n_;

    explicit MontyCommonBase(T modulus) : n_(modulus),
                            neg_inv_n_(negative_inverse_mod_r(modulus)),
                            r_mod_n_(getRModN(modulus)),
                            r_squared_mod_n_(getRSquaredModN(r_mod_n_, modulus))
    {
        // Note: unityValue == (the montgomery form of 1)==(1*R)%n_ == r_mod_n_.
        //
        // getRModN() guarantees the below.  getUnityValue() and
        // getNegativeOneValue() both rely on it.
        HPBC_INVARIANT2(0 < r_mod_n_ && r_mod_n_ < modulus);
        HPBC_INVARIANT2(0 < r_squared_mod_n_ && r_squared_mod_n_ < modulus);
    }
    MontyCommonBase(const MontyCommonBase&) = delete;
    MontyCommonBase& operator=(const MontyCommonBase&) = delete;

//    bool isReduced(T a) const { return (0 <= a && a < getModulus()); }

public:
    // intended for use in postconditions/preconditions
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        V cfx = static_cast<const D*>(this)->getCanonicalForm(x);
        // Any fully reduced value (0 <= value < n_) must be canonical.  Class
        // Derived must be implemented to respect this.
        HPBC_INVARIANT2((0 <= x.get() && x.get() < n_) ? x == cfx : x != cfx);
        bool good = static_cast<const D*>(this)->isValid(x);
        return (x == cfx && good); 
    }

    HURCHALLA_FORCE_INLINE T getModulus() const { return n_; }

    HURCHALLA_FORCE_INLINE V convertIn(T a) const
    {
        return static_cast<const D*>(this)->multiply(V(a), V(r_squared_mod_n_));
    }

    HURCHALLA_FORCE_INLINE V getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_
        HPBC_INVARIANT2(isCanonical(V(r_mod_n_)));
        return V(r_mod_n_);
    }

    HURCHALLA_FORCE_INLINE V getZeroValue() const
    {
        V zero(0); // zeroValue == (0*R)%n_
        HPBC_INVARIANT2(isCanonical(zero));
        return zero;
    } 

    HURCHALLA_FORCE_INLINE V getNegativeOneValue() const
    {
        // We want to get  returnVal = getCanonicalForm(subtract(getZeroValue(),
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
        T ret = n_ - r_mod_n_;
        HPBC_ASSERT2(0 < ret && ret < n_);
        HPBC_INVARIANT2(isCanonical(V(ret)));

        return V(ret);
    }
};


}} // end namespace

#endif
