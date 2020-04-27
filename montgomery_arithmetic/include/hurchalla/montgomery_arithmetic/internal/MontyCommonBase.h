
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_BASE_H_INCLUDED


#include "hurchalla/modular_arithmetic/modular_multiplication.h"
#include "hurchalla/montgomery_arithmetic/internal/negative_inverse_mod_r.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


// Except for where stated otherwise, the algorithms and variable names in this
// file are based on the webpage (in April 2020)
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication

// For discussion purposes throughout this file, given an unsigned integral type
// T, let R = 2^(sizeof(T)*8).  For example: if T is uint64_t, then R = 2^64.
// The name 'R' is based on the wikipedia presentation.


// Returns rModN == R%N.  N is modulus of the monty form.  R described above.
template <typename T>
T getRModN(T N)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");
    precondition2(N % 2 == 1);
    precondition2(N > 1);

    // Assign a tmp T variable rather than directly using the intermediate
    // expression, in order to avoid a negative value (and a wrong answer) in
    // cases where 'N' would be promoted to type 'int'
    T tmp = static_cast<T>(0) - N;

    // Compute R%N.  For example, if R==2^64, arithmetic wraparound behavior of
    // the unsigned integral type T results in (0 - N) representing (2^64 - N).
    // Thus, rModN = R%N == (2^64)%N == (2^64 - N)%N == (0-N)%N
    T rModN = tmp % N;

    // Since N is odd and > 1, N does not divide R==2^x.  Thus, rModN != 0
    postcondition2(0 < rModN && rModN < N);
    return rModN;
}


// Returns rSquaredModN == (R*R)%N.  rModN == R%N  (as computed by getRModN(N))
template <typename T>
T getRSquaredModN(T rModN, T N)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T>::is_modulo, "");
    precondition2(N % 2 == 1);
    precondition2(N > 1);

    // rSqrModN == (R*R)%N == ((R%N)*(R%N))%N == (rModN*rModN)%N
    T rSqrModN = modular_multiplication_prereduced_inputs(rModN, rModN, N);

    // Since N is odd and > 1, N does not divide R*R==2^y.  Thus, rSqrModN != 0
    postcondition2(0 < rSqrModN && rSqrModN < N);
    return rSqrModN;
}


// The class member variable names are based on the webpage
// https://en.wikipedia.org/wiki/Montgomery_modular_multiplication
//
// For discussion purposes, let R = 2^(sizeof(T)*8).  For example: if T is
// uint64_t, then R = 2^64.  This is based on the above presentation.
//
// This base class uses the CRTP idiom
// https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
template <template <typename> class Derived, typename T>
class MontyCommonBase {
    using V = Derived<T>::montvalue_type;
protected:
    const T n_;   // the modulus
    const T nprime_;
    const T r_mod_n_;
    const T r_squared_mod_n_;

    explicit MontyCommonBase(T modulus) : n_(modulus),
                            nprime_(negative_inverse_mod_r(modulus)),
                            r_mod_n_(getRModN(modulus)),
                            r_squared_mod_n_(getRSquaredModN(r_mod_n_, modulus))
    {
        // Note: unityValue == (the montgomery form of 1)==(1*R)%N == rModN
        //
        // getRModN() guarantees the below.  getUnityValue() and
        // getNegativeOneValue() both rely on it.
        invariant(0 < r_mod_n_ && r_mod_n_ < modulus);
        invariant(0 < r_squared_mod_n_ && r_squared_mod_n_ < modulus);
    }
    MontyCommonBase(const MontyCommonBase&) = delete;
    MontyCommonBase& operator=(const MontyCommonBase&) = delete;

    bool isReduced(T a) const { return (0 <= a && a < getModulus()); }

public:
    T getModulus() const { return n_; }

    // intended for use in postconditions/preconditions
    bool isCanonical(V x) const
    {
        V cfx = static_cast<Derived<T>*>(this)->getCanonicalForm(x);
        // Aside from values belonging to special equivalence class 0 (mod N),
        // any fully reduced value (0 < value < N) must be canonical.  Class
        // Derived must be implemented to respect this.  Note: Derived class has
        // freedom to decide whether the canonical form of the equivalence class
        // 0 (mod N)  is equal to 0 or N.
        invariant2((0 < x.get() && x.get() < n_) ? x == cfx :
                                     x.get() == 0 || x.get() == n_ || x != cfx);
        bool good = static_cast<Derived<T>*>(this)->isValid(x);
        return (x == cfx && good); 
    }

    V getUnityValue() const
    {
        // as noted in constructor, unityValue == (1*R)%n_ == r_mod_n_
        postcondition(isCanonical(V(r_mod_n_)));
        return V(r_mod_n_);
    }

    V getNegativeOneValue() const
    {
        // The constructor established the invariant  0 < r_mod_n_ < n_
        // Thus we also know  0 < n_ - r_mod_n_ < n_.
        T ret = n_ - r_mod_n_;
        postcondition(0 < ret && ret < n_);
        // Note: We want returnVal = getCanonicalForm(subtract(getZeroValue(),
        //                                               getUnityValue())).
        //   getZeroValue() returns a value belonging to the equivalence class
        //   0*R (mod N).  This equivalence class can equally be represented by
        //   the value N (mod N).  getUnityValue() returns a value belonging to
        //   the equivalence class 1*R (mod N), which is the same as 
        //   r_mod_n_ (mod N).  The subtract() would return a value belonging to 
        //   the equivalence class (getZeroValue() - r_mod_n_) (mod N).
        //   Substituting representation of the equivalence class for the zero
        //   value, this is equal to (N - r_mod_n_) (mod N).  We know from the
        //   postcondition above that 0 < N - r_mod_n_ < N, which means that
        //   N - r_mod_n_  is fully reduced.  Aside from values belonging to the
        //   special equivalence class 0 (mod N), we always know that any value
        //   that is fully reduced is canonical. Since (N - r_mod_n_) is greater
        //   than zero and fully reduced, it therefore must be canonical.
        postcondition(isCanonical(V(ret)));
        return V(ret);
    }
};


}} // end namespace

#endif
