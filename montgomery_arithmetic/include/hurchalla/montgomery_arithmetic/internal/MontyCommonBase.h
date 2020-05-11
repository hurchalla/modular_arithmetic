
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/montgomery_arithmetic/internal/negative_inverse_mod_r.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include "hurchalla/montgomery_arithmetic/internal/compiler_macros.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// Except for where stated otherwise, the algorithms and variable names in this
// file are based on the webpage (in April 2020)
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication

// For discussion purposes throughout this file, given an unsigned integral type
// T, let R = 2^(std::numeric_limits<T>::digits).  For example: if T is uint64_t
// then R = 2^64.  The name 'R' is based on the wikipedia presentation.


// Returns rModN == R%n.  R is described above.
template <typename T>
T getRModN(T n)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");
    precondition2(n % 2 == 1);
    precondition2(n > 1);

    // Assign a tmp T variable rather than directly using the intermediate
    // expression, in order to avoid a negative value (and a wrong answer) in
    // cases where 'n' would be promoted to type 'int'
    T tmp = static_cast<T>(0) - n;

    // Compute R%n.  For example, if R==2^64, arithmetic wraparound behavior of
    // the unsigned integral type T results in (0 - n) representing (2^64 - n).
    // Thus, rModN = R%n == (2^64)%n == (2^64 - n)%n == (0-n)%n
    T rModN = tmp % n;

    // Since n is odd and > 1, n does not divide R==2^x.  Thus, rModN != 0
    postcondition2(0 < rModN && rModN < n);
    return rModN;
}


// Returns rSquaredModN == (R*R)%n.
// The input parameter rModN must equal R%n (call getRModN(n) to get this value)
template <typename T>
T getRSquaredModN(T rModN, T n)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");
    precondition2(n % 2 == 1);
    precondition2(n > 1);

    // rSqrModN == (R*R)%n == ((R%n)*(R%n))%n == (rModN*rModN)%n
    T rSqrModN = modular_multiplication_prereduced_inputs(rModN, rModN, n);

    // Since n is odd and > 1, n does not divide R*R==2^y.  Thus, rSqrModN != 0
    postcondition2(0 < rSqrModN && rSqrModN < n);
    return rSqrModN;
}




// Make sure MontySqrtRange doesn't derive from this class - it needs different
// functions for convertIn() and getZeroValue().



// The class member variable names are based on the webpage
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
//
// For discussion purposes, let R = 2^(std::numeric_limits<T>::digits).  For
// example if T is uint64_t, then R = 2^64.
//
// This base class uses the CRTP idiom
// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <template <typename> class Derived, typename T>
class MontyCommonBase {
    using D = Derived<T>;
    using V = D::montvalue_type;
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
        // Note: unityValue == (the montgomery form of 1)==(1*R)%N == rModN
        //
        // getRModN() guarantees the below.  getUnityValue() and
        // getNegativeOneValue() both rely on it.
        invariant2(0 < r_mod_n_ && r_mod_n_ < modulus);
        invariant2(0 < r_squared_mod_n_ && r_squared_mod_n_ < modulus);
    }
    MontyCommonBase(const MontyCommonBase&) = delete;
    MontyCommonBase& operator=(const MontyCommonBase&) = delete;

//    bool isReduced(T a) const { return (0 <= a && a < getModulus()); }

public:
    // intended for use in postconditions/preconditions
    FORCE_INLINE bool isCanonical(V x) const
    {
        V cfx = static_cast<D*>(this)->getCanonicalForm(x);
        // Any fully reduced value (0 <= value < N) must be canonical.  Class
        // Derived must be implemented to respect this.
        invariant2((0 <= x.get() && x.get() < n_) ? x == cfx : x != cfx);
        bool good = static_cast<D*>(this)->isValid(x);
        return (x == cfx && good); 
    }

    FORCE_INLINE T getModulus() const { return n_; }

    FORCE_INLINE V convertIn(T a) const
    {
        return static_cast<D*>(this)->multiply(V(a), V(r_squared_mod_n_));
    }

    FORCE_INLINE V getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_
        invariant2(isCanonical(V(r_mod_n_)));
        return V(r_mod_n_);
    }

    FORCE_INLINE V getZeroValue() const
    {
        V zero(0); // zeroValue == (0*R)%N
        invariant2(isCanonical(zero));
        return zero;
    } 

    FORCE_INLINE V getNegativeOneValue() const
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
        invariant2(n_ > r_mod_n_);
        T ret = n_ - r_mod_n_;
        assert_body2(0 < ret && ret < n_);
        invariant2(isCanonical(V(ret)));

        return V(ret);
    }

    FORCE_INLINE V square(V x) const
    {
        return static_cast<D*>(this)->multiply(x, x);
    }
};


}} // end namespace

#endif
