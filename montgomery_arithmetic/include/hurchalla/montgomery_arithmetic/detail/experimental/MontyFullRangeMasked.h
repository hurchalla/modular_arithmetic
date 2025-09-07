// Copyright (c) 2020-2025 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_MASKED_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_FULL_RANGE_MASKED_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyTags.h"
#include "hurchalla/modular_arithmetic/detail/optimization_tag_structs.h"
#include "hurchalla/modular_arithmetic/modular_addition.h"
#include "hurchalla/modular_arithmetic/modular_subtraction.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/safely_promote_unsigned.h"
#include "hurchalla/util/conditional_select.h"
#include "hurchalla/util/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/modular_arithmetic/detail/clockwork_programming_by_contract.h"
#include <type_traits>

namespace hurchalla { namespace detail {


// The name "Fullrange" signifies that there are essentially no preconditions on
// the value of the modulus used in the Montgomery representation.


// struct used internally by MontyFullRangeMasked
template <typename T>
struct MfrmValueTypes {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    struct C; // forward declare C so that V can friend it
    struct FV; // forward declare FV so that V can friend it
    // regular montgomery value type
    struct V {
        // V is a signed integer abstract data type, internally represented by
        // an unusual extended two's complement scheme.  The data member lowbits
        // contains the low bits of the two's complement binary digits, and this
        // is augmented in V by a single implicit extra-high bit.  This extra-
        // high bit is 0 if getmask() == 0, and 1 if getmask() != 0.
        // Here is an example:
        // Note that V depends on T (its full name is MfrmValueTypes<T>::V).
        // So when T is uint64_t, class V is a 65 bit signed integer type, with
        // getbits() giving the low 64 bits and getmask() implicitly providing
        // the high 65th bit.
        HURCHALLA_FORCE_INLINE V() = default;
        template <class PerfTag = CSelectDefaultTag>
        HURCHALLA_FORCE_INLINE void cmov(bool cond, V v)
        {
            // lowbits = (cond) ? v.lowbits : lowbits
            lowbits = ::hurchalla::conditional_select<T, PerfTag>(
                                                      cond, v.lowbits, lowbits);
            // signmask = (cond) ? v.signmask : signmask
            signmask = ::hurchalla::conditional_select<T, PerfTag>(
                                                    cond, v.signmask, signmask);
        }
     protected:
        friend struct C;
        friend struct FV;
        template <typename> friend class MontyFullRangeMasked;
        HURCHALLA_FORCE_INLINE T getbits() const { return lowbits; }
        HURCHALLA_FORCE_INLINE T getmask() const { return signmask; }
        HURCHALLA_FORCE_INLINE V(T lbits, T smask):
                                              lowbits(lbits), signmask(smask) {}
     private:
        T lowbits;
        T signmask;
    };
    // canonical montgomery value type
    struct C : public BaseMontgomeryValue<T> {
        HURCHALLA_FORCE_INLINE C() = default;
        // implicitly convert from canonical value C to montgomery value V
        HURCHALLA_FORCE_INLINE operator V() const { return V(this->get(), 0); }
        HURCHALLA_FORCE_INLINE friend bool operator==(const C& x, const C& y)
            { return x.get() == y.get(); }
        HURCHALLA_FORCE_INLINE friend bool operator!=(const C& x, const C& y)
            { return !(x == y); }
     protected:
        template <typename> friend class MontyFullRangeMasked;
        template <template<class> class, template<class> class, typename>
          friend class MontyCommonBase;
        HURCHALLA_FORCE_INLINE explicit C(T a) : BaseMontgomeryValue<T>(a) {}
    };
    // fusing montgomery value (addend/subtrahend for fmadd/fmsub)
    struct FV : public BaseMontgomeryValue<T> {
        HURCHALLA_FORCE_INLINE FV() = default;
        // implicitly convert from fusing value FV to montgomery value V
        HURCHALLA_FORCE_INLINE operator V() const { return V(this->get(), 0); }
     protected:
        template <typename> friend class MontyFullRangeMasked;
        HURCHALLA_FORCE_INLINE explicit FV(T a) : BaseMontgomeryValue<T>(a) {}
    };
};


// Let the theoretical constant R = 1<<(ut_numeric_limits<T>::digits).
template <typename T>
class MontyFullRangeMasked final :
           public MontyCommonBase<MontyFullRangeMasked, MfrmValueTypes, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    using BC = MontyCommonBase<::hurchalla::detail::MontyFullRangeMasked,
                               ::hurchalla::detail::MfrmValueTypes, T>;
    using BC::n_;
    using typename BC::V;
    using typename BC::C;
    using FV = typename MfrmValueTypes<T>::FV;
    using SV = V;
 public:
    using MontyTag = TagMontyFullrangeMasked;
    using uint_type = T;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;
    using squaringvalue_type = SV;

    explicit MontyFullRangeMasked(T modulus) : BC(modulus) {}

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return (ut_numeric_limits<T>::max() % 2 == 0) ?
                static_cast<T>(ut_numeric_limits<T>::max() - 1) :
                ut_numeric_limits<T>::max();
    }

    HURCHALLA_FORCE_INLINE V negate(V x) const
    {
        return subtract(BC::getZeroValue(), x, LowuopsTag());
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isValid(x));
        T tmpn = static_cast<T>(x.getmask() & n_);
        C result = C(static_cast<T>(x.getbits() + tmpn));
        HPBC_CLOCKWORK_POSTCONDITION2(result.get() < n_);
        return result;
    }

    // Note: internal to MontyFullRangeMasked, the contents of FusingValue (FV)
    // and CanonicalValue (C) variables are interchangeable.  Other Monty types
    // use FV and C as completely distinct types, and so for genericity we
    // always present C and FV to the outside world as being unrelated.
    HURCHALLA_FORCE_INLINE FV getFusingValue(V x) const
    {
        C cv = getCanonicalValue(x);
        return FV(cv.get());
    }
    using BC::fmadd;
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, FV fv, PTAG) const
    {
        C cv = C(fv.get());
        return fmadd(x, y, cv, PTAG());
    }
    using BC::fmsub;
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, FV fv, PTAG) const
    {
        C cv = C(fv.get());
        return fmsub(x, y, cv, PTAG());
    }

    HURCHALLA_FORCE_INLINE V add(V x, C cy) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isValid(x));
        HPBC_CLOCKWORK_PRECONDITION2(0 <= cy.get() && cy.get() < n_);
        HPBC_CLOCKWORK_ASSERT2(0 <= cy.get() && cy.get() < n_);
#if defined(_MSC_VER)
// MSVC only section-- for fewer uops and lower latency (x64 and arm).  Either
// section would be correct, but msvc generates better machine code here
        C cx = getCanonicalValue(x);
        HPBC_CLOCKWORK_ASSERT2(0 <= cx.get() && cx.get() < n_);
        T resultval = ::hurchalla::modular_addition_prereduced_inputs(
                                                        cx.get(), cy.get(), n_);
        T result_smask = 0;
#else
        // if x is negative, set tmpn = n, otherwise set tmpn = 0.
        T tmpn = static_cast<T>(x.getmask() & n_);
        // x.getbits() has range [-n, n),  so x.getbits()-n has range [-2n, 0)
        T tmpx = static_cast<T>(x.getbits() - n_);       // has range [-2n, 0)
        tmpx = static_cast<T>(tmpx + tmpn);              // has range [-n, 0)
        T resultval = static_cast<T>(tmpx + cy.get());   // has range [-n, n)
        bool b = (resultval < tmpx);
        T result_smask = static_cast<T>(static_cast<T>(b) - static_cast<T>(1));
#endif
        V result = V(resultval, result_smask);
        HPBC_CLOCKWORK_POSTCONDITION2(isValid(result));
        return result;
    }
    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        C cy = getCanonicalValue(y);
        return add(x, cy);
    }
    HURCHALLA_FORCE_INLINE C add(C cx, C cy) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(cx.get() < n_);
        HPBC_CLOCKWORK_PRECONDITION2(cy.get() < n_);
        T result = ::hurchalla::modular_addition_prereduced_inputs(
                                                        cx.get(), cy.get(), n_);
        HPBC_CLOCKWORK_POSTCONDITION2(result < n_);
        return C(result);
    }

    template <class PTAG>
    HURCHALLA_FORCE_INLINE V subtract(V x, C cy, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isValid(x));
        HPBC_CLOCKWORK_PRECONDITION2(0 <= cy.get() && cy.get() < n_);
        C cx = getCanonicalValue(x);
        HPBC_CLOCKWORK_ASSERT2(0 <= cx.get() && cx.get() < n_);
        T resultval = static_cast<T>(cx.get() - cy.get());
        bool b = (cx.get() < cy.get());
        T result_smask = static_cast<T>(0 - static_cast<T>(b));
        V result = V(resultval, result_smask);
        HPBC_CLOCKWORK_POSTCONDITION2(isValid(result));
        return result;
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V subtract(C cx, V y, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(0 <= cx.get() && cx.get() < n_);
        HPBC_CLOCKWORK_PRECONDITION2(isValid(y));
        C cy = getCanonicalValue(y);
        HPBC_CLOCKWORK_ASSERT2(0 <= cy.get() && cy.get() < n_);
        T resultval = static_cast<T>(cx.get() - cy.get());
        bool b = (cx.get() < cy.get());
        T result_smask = static_cast<T>(0 - static_cast<T>(b));
        V result = V(resultval, result_smask);
        HPBC_CLOCKWORK_POSTCONDITION2(isValid(result));
        return result;
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE V subtract(V x, V y, PTAG) const
    {
        C cy = getCanonicalValue(y);
        return subtract(x, cy, PTAG());
    }
    template <class PTAG>
    HURCHALLA_FORCE_INLINE C subtract(C cx, C cy, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(cx.get() < n_);
        HPBC_CLOCKWORK_PRECONDITION2(cy.get() < n_);
        T result = ::hurchalla::modular_subtraction_prereduced_inputs
                                    <decltype(n_),PTAG>(cx.get(), cy.get(), n_);
        HPBC_CLOCKWORK_POSTCONDITION2(result < n_);
        return C(result);
    }

    template <typename J, typename K>
    HURCHALLA_FORCE_INLINE V unordered_subtract(J x, K y) const
    {
        return subtract(x, y, LowuopsTag());
    }

    HURCHALLA_FORCE_INLINE V two_times(V x) const
    {
        return add(x, x);
    }
    HURCHALLA_FORCE_INLINE C two_times(C cx) const
    {
        return add(cx, cx);
    }


    HURCHALLA_FORCE_INLINE SV getSquaringValue(V x) const
    {
        static_assert(std::is_same<V, SV>::value, "");
        return x;
    }
    HURCHALLA_FORCE_INLINE SV squareSV(SV sv) const
    {
        static_assert(std::is_same<V, SV>::value, "");
        return BC::square(sv, LowlatencyTag());
    }
    HURCHALLA_FORCE_INLINE V squareToMontgomeryValue(SV sv) const
    {
        static_assert(std::is_same<V, SV>::value, "");
        return BC::square(sv, LowlatencyTag());
    }
    HURCHALLA_FORCE_INLINE V getMontgomeryValue(SV sv) const
    {
        static_assert(std::is_same<V, SV>::value, "");
        return sv;
    }

private:
    // functions called by the 'curiously recurring template pattern' base (BC).
    friend BC;

    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(bool& resultIsZero, T u_hi, T u_lo, PTAG) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        bool isNegative;
        T resultval = ::hurchalla::REDC_incomplete(
                                        isNegative, u_hi, u_lo, n_, BC::inv_n_);
        T result_smask = static_cast<T>(0 - static_cast<T>(isNegative));
        resultIsZero = (resultval == 0);
        V result = V(resultval, result_smask);
        HPBC_CLOCKWORK_POSTCONDITION2(isValid(result));
        return result;
    }
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(T u_hi, T u_lo, PTAG) const
    {
        bool resultIsZero;  // ignored
        return montyREDC(resultIsZero, u_hi, u_lo, PTAG());
    }

    // This function computes the full two-word product of x*x.  It returns the
    // high word of x*x, and writes the low word of x*x to u_lo.
    HURCHALLA_FORCE_INLINE T squareToHiLo(T& u_lo, V x) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isValid(x));
        T a = x.getbits();
        T umlo;
        T umhi = ::hurchalla::unsigned_multiply_to_hilo_product(umlo, a, a);
        T masked_a = static_cast<T>(x.getmask() & a);
        T result_hi = static_cast<T>(umhi - static_cast<T>(2) * masked_a);
        u_lo = umlo;
        // Complete details are in the proof below, but roughly what we do here
        // is get a*a as a two-word product (umhi, umlo).  We let s == 1 if x is
        // negative (in which case x.getmask() is all ones) and s == 0 otherwise
        // (in which case x.getmask() is all zeros).  Note that s*s == s.  We
        // then use the identity  x == a - s*R  (due to x's two's complement
        // representation) to get  x*x == a*a - 2*s*a*R + s*s*R*R  and simplify
        // to  x*x == a*a + s*(R-2*a)*R.  The term  s*(R-2*a)*R  has a low word
        // of 0, and a high word (which might borrow) of  s*(R-2*a).  The term
        // a*a  has a low word and high word respectively of umlo and umhi.
        // When we add the low words of the two terms, we simply get umlo (there
        // is no carry).  When we add the high words we get  umhi + s*(R-2*a).
        // The proof below shows that this result neither borrows nor carries
        // (i.e. 0 <= umhi + s*(R-2*a) < R).  Since math in C++ is always
        // performed mod R on (unsigned) type T expressions, we can calculate
        // umhi + s*(R-2*a)  using  umhi - s*2*a.  Simplifying further, we get
        // high_word_result = umhi - 2*(x.getmask() & a).

        HPBC_CLOCKWORK_POSTCONDITION2(result_hi < n_);
        // Since we have class invariant n < R, and our isValid(x) precondition
        // requires -n < x < n,  we know  x*x < n*n < n*R.
        // Thus  result_hi*R + umlo == x*x < n*R, and since  0 <= umlo, we know
        // result_hi*R <= result_hi*R + umlo < n*R.  Since
        // result_hi*R < n*R,  we know  result_hi < n.
        return result_hi;

        // Proof that the algorithm expressed in code above is correct:
        // (First read the comments in MfrmValueTypes::V at top of this file)
        //
        // In the discussion that follows we'll mostly ignore C++ types, by
        // treating the variables as mathematical integers belonging to the set
        // Z of all possible integers (we can note that all C++ integer types
        // are subsets of Z).
        // Let the variable s equal 1 if x is negative, and otherwise equal 0.
        // Recall that the theoretical constant
        // R = pow(2, ut_numeric_limits<T>::digits).
        // Let the variable 'a' equal x.getbits().  Since x.getbits() is an
        // unsigned (type T) two's complement binary value that represents a
        // signed integer, we can write x as
        // x == a - s*R.
        // Thus  x*x == (a - s*R)*(a - s*R) == a*a - 2*s*R*a + s*s*R*R.
        // Since s is either 0 or 1, s*s == s, and thus
        // x*x == a*a - 2*s*R*a + s*R*R
        // x*x == a*a + s*(R-2*a)*R
        //
        // Since the square of any integer is >= 0, we know  x*x >= 0.
        // Since we have class invariant n < R, and our isValid(x) precondition
        // requires -n <= x < n,  we know  -R < x < R.  Thus x*x < R*R, and thus
        // 0 <= x*x < R*R.
        //
        // Because 'a' equals x.getbits() and x.getbits() returns type T, we
        // know that  0 <= a < R.  And thus
        // 0 <= a*a < R*R.
        // Let asqrHi and asqrLo be values satisfying
        // a*a == asqrHi*R + asqrLo,  with
        // 0 <= asqrHi < R  and  0 <= asqrLo < R.
        // Since 0 <= a*a < R*R, we know that asqrHi and asqrLo exist.
        // We can replace a*a in the x*x equality above, giving us
        // x*x == (asqrHi*R + asqrLo) + s*(R-2*a)*R
        // x*x == (asqrHi + s*(R-2*a))*R + asqrLo
        //
        // Assume temporarily that  asqrHi + s*(R-2*a) <= -1.
        //   Then (asqrHi + s*(R-2*a))*R <= -R
        //   and (asqrHi + s*(R-2*a))*R + asqrLo <= -R + asqrLo
        //   Since  x*x == (asqrHi + s*(R-2*a))*R + asqrLo,  x*x <= -R + asqrLo.
        //   Since asqrLo < R,  asqrLo <= R - 1, and thus
        //   x*x <= -R + asqrLo <= -R + (R - 1) == -1.
        //   But x*x <= -1 is a contradiction because we know  0 <= x*x < R*R.
        // Assume temporarily that  asqrHi + s*(R-2*a) >= R.
        //   Then (asqrHi + s*(R-2*a))*R >= R*R.  And
        //   (asqrHi + s*(R-2*a))*R + asqrLo >= R*R + asqrLo.
        //   Since  x*x == (asqrHi + s*(R-2*a))*R + asqrLo,  x*x >= R*R + asqrLo
        //   Since we know  asqrLo >= 0,  x*x >= R*R + asqrLo >= R*R.
        //   But x*x >= R*R is a contradiction because we know  0 <= x*x < R*R.
        // Therefore both assumptions are false, and so
        // 0 <= asqrHi + s*(R-2*a) < R.
        //
        // Recall that asqrHi and asqrLo are simply any values that satisfy
        // asqrHi*R + asqrLo == a*a, with 0 <= asqrHi < R  and  0 <= asqrLo < R.
        // When we compute
        // T umlo; T umhi = unsigned_multiply_to_hilo_product(umlo, a, a);
        // umhi and umlo satisfy these requirements, since all type T variables
        // have bounds [0, R), and unsigned_multiply_to_hilo_product() computes
        // the full two-word product of a*a.
        // Therefore, since
        // x*x == (asqrHi + s*(R-2*a))*R + asqrLo,  we can substitute and get
        // x*x == (umhi + s*(R-2*a))*R + umlo.  Likewise, substitution gives us
        // 0 <= umhi + s*(R-2*a) < R.
        //
        // Let "%%" be a conceptual modulo operator that always produces the
        // non-negative remainder for its result.  For example,  8 %% 5 == 3,
        // and  -8 %% 5 == 2.
        // Since  0 <= umhi + s*(R-2*a) < R,  we know
        // umhi + s*(R-2*a) == (umhi + s*(R-2*a)) %% R.  Furthermore,
        // umhi + s*(R-2*a) == (umhi + (s*(R-2*a))%%R) %% R
        //                  == (umhi + (s*((R-2*a)%%R))%%R) %% R
        //                  == (umhi + (s*((-2*a)%%R))%%R) %% R
        // In C and C++, mathematical operations (e.g. addition and subtraction
        // and multiplication) on type T (an unsigned integer type) values are
        // always performed modulo R, and always satisfy  0 <= finalresult < R.
        // This is the same as applying the "%%"" operator to every math result.
        // Therefore, we can calculate the expression
        // (umhi + (s*((-2*a)%%R))%%R) %% R
        // on a computer in C++ via the expression
        // T result_hi = static_cast<T>(umhi) +
        //                static_cast<T>(s)*static_cast<T>(-2)*static_cast<T>(a)
        // and therefore result_hi is equal to the mathematical integer (set Z)
        // result of  umhi + s*(R-2*a).
        //
        // Since static_cast<T>(a) is equal to x.getbits(), we can calculate
        // static_cast<T>(-2)*static_cast<T>(a)  using
        // T neg2a = static_cast<T>(-2) * x.getbits();
        // And umhi is already type T so there is no need to cast it.  Therefore
        // we can simplify result_hi by using
        // T result_hi = umhi + static_cast<T>(s) * neg2a;
        // 
        // All that remains is s.  As specified earlier, if x < 0 then s == 1,
        // and if x >= 0 then s == 0.  x.getmask() returns a type T bitmask
        // consisting of all ones when x < 0, and returns a bitmask consisting
        // of all zeros when x >= 0.
        // Therefore  x.getmask() & neg2a == static_cast<T>(s) * neg2a.
        // And therefore we can simplify result_hi further, using
        // T result_hi = umhi + (x.getmask() & neg2a);
        //
        // As shown above, result_hi == umhi + s*(R-2*a),  and
        // x*x == (umhi + s*(R-2*a))*R + umlo,  and so
        // x*x == result_hi*R + umlo.
        // Our function (squareToHiLo()) needs to provide the unsigned two-word
        // full product of x*x, and we can see that result_hi and umlo are this
        // desired two-word product.
        //
        // We can express all of this in code via
        // T a = x.getbits();
        // T umlo;
        // T umhi = ::hurchalla::unsigned_multiply_to_hilo_product(umlo, a, a);
        // T neg2a = static_cast<T>(-2) * a;
        // T result_hi = umhi + (x.getmask() & neg2a);
        // u_lo = umlo;
        // return result_hi;
        //
        // This function uses a slightly optimized version of the code above.
    }


    // Let u be an arbitrary double-word value that is congruent (mod n_) to the
    // product of x and y, and that satisfies 0 <= u < n_*R.  This function
    // returns the high word of u, and writes the low word of u to u_lo.
    // Important performance note: parameter 'y' adds extra latency to this
    // function if it is part of a loop-carried (or other) dependency chain, and
    // so when you need to multiply one variable that is part of a dependency
    // chain with another variable that is not, you should pass the dependency
    // variable as parameter 'x' and the non-dependency variable as 'y'.
    HURCHALLA_FORCE_INLINE T multiplyToHiLo(
                                     T& HURCHALLA_RESTRICT u_lo, V x, V y) const
    {
        HPBC_CLOCKWORK_PRECONDITION2(isValid(x));  // x has range [-n, n)
        HPBC_CLOCKWORK_PRECONDITION2(isValid(y));
        C cy = getCanonicalValue(y);     // cy has range [0, n)
        T a = x.getbits();
        T b = cy.get();                  // b has range [0, n)
        HPBC_CLOCKWORK_ASSERT2(0 <= b && b < n_);

        T u_hi = ::hurchalla::unsigned_multiply_to_hilo_product(u_lo, a, b);
#if 1
// This section follows the algorithm described in the proof.  It will likely
// perform well on all architectures due to its predictable conditional branch.
// And it's well suited to RISC-V, since RISC-V has no conditional move/select.
        if (b != 0) {    // usually predictable
           T tmp = static_cast<T>(n_ - b);
           T maskedtmp = static_cast<T>(x.getmask() & tmp);
           u_hi = static_cast<T>(u_hi + maskedtmp);
        }
#else
// This section is an optimization(?) of the code above, relying on conditional
// move/select; perf measurements would be needed to know if it has any benefit.
        T sumn = static_cast<T>(u_hi + n_);
        T tmp = static_cast<T>(x.getmask() & b);
        u_hi = conditional_select(tmp != 0, sumn, u_hi); //uhi=(tmp!=0)?sumn:uhi
        u_hi = static_cast<T>(u_hi - tmp);
#endif
        HPBC_CLOCKWORK_POSTCONDITION2(u_hi < n_);
        return u_hi;

        // Proof that the algorithm expressed in code above is correct:
        // (First read the comments in MfrmValueTypes::V at top of this file)
        //
        // In the discussion that follows we'll mostly ignore C++ types, by
        // treating the variables as mathematical integers belonging to the set
        // Z of all possible integers (we can note that all C++ integer types
        // are subsets of Z).
        // Let b = getCanonicalValue(y).get().  Due to getCanonicalValue's post-
        // conditions, we know that 0 <= b < n, with b ≡ y (mod n); and thus
        // also  x*y ≡ x*b (mod n).
        // Let s = 1 if x is negative, and otherwise let s = 0.
        // Recall that the theoretical constant
        // R = pow(2, ut_numeric_limits<T>::digits).
        // Let the variable  a = x.getbits().  Since x.getbits() is an
        // unsigned (type T) two's complement binary value that represents a
        // signed integer, we can write x as
        // x == a - s*R.
        // Thus  x*b == (a - s*R)*b == a*b - s*b*R.
        // Because 'a' equals x.getbits() and x.getbits() returns type T, we
        // know that  0 <= a < R.  And since 0 <= b < n < R, we know that
        // 0 <= a*b < R*R.
        // 
        // We'll separately handle the two possible cases, b==0 and b>0.
        // Case 1, b > 0:
        //   Since x*b == a*b - s*b*R  and  x*y ≡ x*b (mod n),  we know
        //   x*y ≡ (a*b - s*b*R) ≡ (a*b + (n-b)*s*R)  (mod n).
        //   Let u = a*b + (n-b)*s*R.  We can see  u ≡ x*y  (mod n).
        //      Since we have b < n, we know (n-b) > 0.  Furthermore since
        //      s >= 0 and R > 0, we know  (n-b)*s*R >= 0.  Finally, since
        //      a*b >= 0, we know  a*b + (n-b)*s*R >= 0.
        //   Thus  u = a*b + (n-b)*s*R >= 0.
        //      Since a is type T, a < R; and our case specifies b > 0.  Thus
        //      a*b < b*R.  Therefore  u = a*b + (n-b)*s*R  <  b*R + (n-b)*s*R.
        //      Thus  u < (b - s*b + s*n)*R.
        //      Assume temporarily that  s == 1:  Then u < n*R.
        //      Assume temporarily that  s == 0:  Then u < b*R.  Since b < n
        //         and R > 0,  b*R < n*R.  And so  u < n*R.
        //      s can only be 0 or 1, and for both assumptions, we have u < n*R.
        //   Thus u < n*R.
        //   Put together, we know  0 <= u < n*R.
        //   Since also n < R and R > 0,  0 <= u < n*R < R*R.  Thus u will fit
        //   within a two-word (type T words) representation, and can be
        //   computed with double-word (multiprecision with a fixed size of two
        //   words) computer arithmetic performed modulo R*R.
        //   We showed that u is congruent to x*y (mod n), and that
        //   0 <= u < n*R.  Thus the value u (specifically, its high word and
        //   low word) will satisfy this function's postconditions when b > 0.
        // Case 2, b == 0:
        //   x*b == a*b - s*b*R == 0 == a*b.  Thus x*y ≡ x*b ≡ a*b (mod n).
        //   Let v = a*b.  We can see that v is congruent to x*y (mod n), and
        //   since v == 0 and n > 0 and R > 0, v satisfies 0 <= v < n*R.  Thus
        //   the value v (specifically, its high word and low word) will satisfy
        //   this function's postconditions when b == 0.
        //   We can use unsigned_multiply_to_hilo_product() to compute the two-
        //   word product v = a*b.  Note: we could simply set u = 0 without
        //   using the multiply, but case b==0 is a rare/unusual case, and we
        //   get smaller code by always performing the multiply a*b, and using
        //   it for both case 1 and case 2.
    }

    // the valid range for x:  -n <= x < n.
    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        return ((x.getmask() != 0) ? static_cast<T>(0 - x.getbits()) <= n_
                                   : x.getbits() < n_);
    }
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        return (x.getmask() == 0);
    }
    // get a Natural number (i.e. number >= 0) congruent to x (mod n)
    HURCHALLA_FORCE_INLINE T getNaturalEquivalence(V x) const
    {
        C cx = getCanonicalValue(x);
        return cx.get();
    }
};


}} // end namespace

#endif
