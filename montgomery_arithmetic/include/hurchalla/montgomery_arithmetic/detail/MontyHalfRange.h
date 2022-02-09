// Copyright (c) 2020-2022 Jeffrey Hurchalla.
/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 */

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_HALF_RANGE_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_HALF_RANGE_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/low_level_api/optimization_tag_structs.h"
#include "hurchalla/montgomery_arithmetic/low_level_api/REDC.h"
#include "hurchalla/montgomery_arithmetic/detail/MontyCommonBase.h"
#include "hurchalla/util/traits/ut_numeric_limits.h"
#include "hurchalla/util/traits/extensible_make_signed.h"
#include "hurchalla/util/signed_multiply_to_hilo_product.h"
#include "hurchalla/util/compiler_macros.h"
#include "hurchalla/util/programming_by_contract.h"
#include <type_traits>


namespace hurchalla { namespace detail {


// The name "Halfrange" signifies that the modulus must be less than R/2, where
// R = 1<<(ut_numeric_limits<T>::digits).  For example, if T is uint64_t then
// R = 1<<64 and R/2 == 1<<63, and thus it would require  modulus < (1<<63).


struct TagMontyHalfrange final {};


// struct used internally by MontyHalfRange
template <typename T>
struct MontyHRValueTypes {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    using SignedT = typename extensible_make_signed<T>::type;
    // Assert that casting a SignedT value to an (unsigned) T results in
    // (SignedT_value + R) % R.
    static_assert(static_cast<T>(1 + static_cast<T>(static_cast<SignedT>(-1)))
                  == 0, "");
    // Assert SignedT uses two's complement representation, maybe redundant
    static_assert(static_cast<SignedT>(-1) == ~(static_cast<SignedT>(0)), "");

    struct C; // forward declare C so that V can friend it
    // regular montgomery value type
    struct V : public BaseMontgomeryValue<SignedT> {
        HURCHALLA_FORCE_INLINE V() = default;
     protected:
        friend struct C;
        template <typename> friend class MontyHalfRange;
        HURCHALLA_FORCE_INLINE
        explicit V(SignedT a) : BaseMontgomeryValue<SignedT>(a) {}
    };
    // canonical montgomery value type
    struct C : public BaseMontgomeryValue<T> {
        HURCHALLA_FORCE_INLINE C() = default;
        // implicitly convert from canonical value C to montgomery value V
        HURCHALLA_FORCE_INLINE operator V() const
        {
            T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
            HPBC_PRECONDITION2(0 <= this->get() && this->get() < Rdiv2);
            return V(static_cast<SignedT>(this->get()));
        }
        HURCHALLA_FORCE_INLINE friend bool operator==(const C& x, const C& y)
            { return x.get() == y.get(); }
        HURCHALLA_FORCE_INLINE friend bool operator!=(const C& x, const C& y)
            { return !(x == y); }
     protected:
        template <typename> friend class MontyHalfRange;
        template <template<class> class, template<class> class, typename>
          friend class MontyCommonBase;
        HURCHALLA_FORCE_INLINE explicit C(T a) : BaseMontgomeryValue<T>(a) {}
    };
    // fusing montgomery value (addend/subtrahend for fmadd/fmsub)
    struct FV : public V {
        HURCHALLA_FORCE_INLINE FV() = default;
     protected:
        template <typename> friend class MontyHalfRange;
        HURCHALLA_FORCE_INLINE explicit FV(SignedT a) : V(a) {}
    };
};


// Let the theoretical constant R = 1<<(ut_numeric_limits<T>::digits).
template <typename T>
class MontyHalfRange final :
                  public MontyCommonBase<MontyHalfRange, MontyHRValueTypes, T> {
    static_assert(ut_numeric_limits<T>::is_integer, "");
    static_assert(!(ut_numeric_limits<T>::is_signed), "");
    static_assert(ut_numeric_limits<T>::is_modulo, "");
    static_assert(ut_numeric_limits<T>::digits >= 1, "");

    using S = typename MontyHRValueTypes<T>::SignedT;
    static_assert(ut_numeric_limits<S>::is_signed, "");

    using BC = MontyCommonBase<::hurchalla::detail::MontyHalfRange,
                               ::hurchalla::detail::MontyHRValueTypes, T>;
    using BC::n_;
    using typename BC::V;
    using typename BC::C;
    using FV = typename MontyHRValueTypes<T>::FV;
 public:
    using MontyTag = TagMontyHalfrange;
    using uint_type = T;
    using montvalue_type = V;
    using canonvalue_type = C;
    using fusingvalue_type = FV;

    explicit MontyHalfRange(T modulus) : BC(modulus)
    {
        // MontyHalfRange requires  modulus < R/2
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_PRECONDITION2(modulus < Rdiv2);
    }
    MontyHalfRange(const MontyHalfRange&) = delete;
    MontyHalfRange& operator=(const MontyHalfRange&) = delete;

    static HURCHALLA_FORCE_INLINE constexpr T max_modulus()
    {
        return static_cast<T>((static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1)) - 1);
    }

    HURCHALLA_FORCE_INLINE C getCanonicalValue(V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_INVARIANT2(n_ < Rdiv2);
        S sx = x.get();
        T tx = static_cast<T>(sx);
        T tc = static_cast<T>(tx + n_);
#if defined(__clang__)
        // Since isValid() required sx >= -((S)n) and sx < n, the sum for tc
        // above will have carried (wrapped-around) if sx<0.  Thus testing
        // (tc < n) is equivalent to a test for (sx < 0).  And consequently,
        // testing (tc >= n) is equivalent to a test for (sx >= 0).
        // It is in theory the better test, because it lets the compiler use
        // the flags already set by (tx+n), if the compiler is smart.  Only
        // clang appears to be smart enough to do this though; for x64 and
        // arm32/64: gcc, msvc, icc produce better asm with the #else.
// TODO: possibly offer predefined macro for using _asm for this cmov
        HURCHALLA_CMOV(tc >= n_, tc, tx);  // if tc>=n_, set tc = tx.
#else
// TODO: possibly offer predefined macro for using _asm for this cmov
        HURCHALLA_CMOV(sx >= 0, tc, tx);   // if sx>=0, set tc = tx.
#endif
        HPBC_POSTCONDITION2(tc < n_);
        return C(tc);
    }

    HURCHALLA_FORCE_INLINE FV getFusingValue(V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_INVARIANT2(n_ < Rdiv2);
        C cx = getCanonicalValue(x);
        T a = cx.get();
        HPBC_ASSERT2(a < n_);  // a has range [0, n_)

        HPBC_INVARIANT2(n_ % 2 == 1);
        T half_n_floor = static_cast<T>(n_ >> 1);     // == (n_ - 1) / 2

        T tmp = static_cast<T>(a - n_);
        // static_cast<S>(tmp) has range [-n, 0)
        HPBC_ASSERT2(-static_cast<S>(n_) <= static_cast<S>(tmp)
                     && static_cast<S>(tmp) < 0);
        T result = a;
        // if ((n_-1)/2) < a, set result = tmp
        HURCHALLA_CMOV(half_n_floor < a, result, tmp);
        // Assume a <= half_n_floor; then trivially we would have
        //    0 <= a <= (n_ - 1) / 2
        //    Given our assumption, we would have  0 <= result <= (n_ - 1) / 2
        // Assume a > half_n_floor; then we would have
        //    static_cast<S>(a - n_) > static_cast<S>(half_n_floor - n_)
        //    static_cast<S>(tmp) > static_cast<S>(((n_ - 1) / 2) - n_)
        //    static_cast<S>(tmp) >= static_cast<S>(((n_ - 1) / 2) - n_ + 1)
        //    static_cast<S>(tmp) >= -static_cast<S>((n_ - 1) / 2)
        //    From an assertion above, we already know static_cast<S>(tmp) < 0
        //    and so together we would have
        //    -static_cast<S>((n_ - 1) / 2) <= static_cast<S>(tmp) < 0.
        //    Therefore given our assumption, we would have
        //    -static_cast<S>((n_ - 1) / 2) <= static_cast<S>(result) < 0.
        // Since one of the two assumptions must be true, we have a possible
        // range for result of
        // -static_cast<S>((n_-1)/2) <= static_cast<S>(result) <= (n_-1)/2
        HPBC_POSTCONDITION2(
                     -static_cast<S>(half_n_floor) <= static_cast<S>(result)
                     && static_cast<S>(result) <= static_cast<S>(half_n_floor));
        return FV(static_cast<S>(result));
    }

    using BC::fmsub;
    // Multiplies two montgomery values x and y, and then subtracts the fusing-
    // value fv from the product.  Returns the resulting montgomery value.
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmsub(V x, V y, FV fv, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        HPBC_INVARIANT2(n_ % 2 == 1);
        // note that the constructor also established the invariant n < R/2.
        HPBC_PRECONDITION2(-(static_cast<S>((n_ - 1)/2)) <= fv.get()
                           && fv.get() <= static_cast<S>((n_ - 1)/2));
        T u_lo;
        S u_hi = signed_multiply_to_hilo_product(u_lo, x.get(), y.get());

        // Performing the modular sub prior to the REDC will always give
        // equivalent results to performing the REDC and then the modular
        // subtraction.  See fmadd() in MontyCommonBase.h for a proof which
        // could be adapted to the subtraction in here.

        // By isValid(), x and y satisfy  -n <= x < n  and  -n <= y < n
        // Thus x*y <= n*n < n*R/2 == ((n-1)/2)*R + R/2
        // And since x*y == u_hi*R + u_lo, we know  u_hi <= ((n-1)/2).
        // Also,  x*y > -n*n > -n*R/2 == -((n+1)/2)*R + R/2.
        // And since x*y == u_hi*R + u_lo, we know  u_hi >= -((n+1)/2).
        // Putting this together, we know:  -((n+1)/2) <= u_hi <= ((n-1)/2).
        HPBC_ASSERT2(-(static_cast<S>((n_+1)/2)) <= u_hi
                     && u_hi <= static_cast<S>((n_-1)/2));
        // Since we have precondition that  -((n-1)/2) <= fv.get() <= (n-1)/2,
        // we know that
        // -((n+1)/2) - (n-1)/2 <= (u_hi - fv.get()) <= ((n-1)/2) - (-((n-1)/2))
        // Thus  -n <= (u_hi - fv.get()) <= n - 1 < n.
        // Also, since n < R/2, we have  -R/2 < (u_hi - fv.get()) < R/2,  which
        // means (u_hi - fv.get()) should never overflow type S.
        u_hi = static_cast<S>(u_hi - fv.get());
        HPBC_ASSERT2(-static_cast<S>(n_) <= u_hi && u_hi < static_cast<S>(n_));
        V vu_hi(u_hi);
        HPBC_ASSERT2(isValid(vu_hi));
        C cu_hi = getCanonicalValue(vu_hi);
        T tu_hi = cu_hi.get();
        HPBC_ASSERT2(tu_hi < n_);

        V result = montyREDC(tu_hi, u_lo, PTAG());

        HPBC_POSTCONDITION2(isValid(result));
        return result;
    }

    using BC::fmadd;
    template <class PTAG>   // Performance TAG (see optimization_tag_structs.h)
    HURCHALLA_FORCE_INLINE V fmadd(V x, V y, FV fv, PTAG) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        HPBC_INVARIANT2(n_ % 2 == 1);
        // note that the constructor also established the invariant n < R/2.
        HPBC_PRECONDITION2(-(static_cast<S>((n_ - 1)/2)) <= fv.get()
                           && fv.get() <= static_cast<S>((n_ - 1)/2));
        T u_lo;
        S u_hi = signed_multiply_to_hilo_product(u_lo, x.get(), y.get());

        // Performing the modular add prior to the REDC will always give
        // equivalent results to performing the REDC and then the modular
        // addition.  See fmadd() in MontyCommonBase.h for a proof which could
        // be adapted to the modular addition in here.

        // See fmsub() above for why the following assert is true
        HPBC_ASSERT2(-(static_cast<S>((n_+1)/2)) <= u_hi
                     && u_hi <= static_cast<S>((n_-1)/2));
        // Since we have precondition that  -((n-1)/2) <= fv.get() <= (n-1)/2,
        // we know that
        // -((n+1)/2) + -((n-1)/2) <= (u_hi + fv.get()) <= ((n-1)/2) + (n-1)/2
        // Thus  -n <= (u_hi + fv.get()) <= n - 1 < n.
        // Also, since n < R/2, we have  -R/2 < (u_hi + fv.get()) < R/2,  which
        // means (u_hi + fv.get()) should never overflow type S.
        u_hi = static_cast<S>(u_hi + fv.get());
        HPBC_ASSERT2(-static_cast<S>(n_) <= u_hi && u_hi < static_cast<S>(n_));

        V vu_hi(u_hi);
        HPBC_ASSERT2(isValid(vu_hi));
        C cu_hi = getCanonicalValue(vu_hi);
        T tu_hi = cu_hi.get();
        HPBC_ASSERT2(tu_hi < n_);

        V result = montyREDC(tu_hi, u_lo, PTAG());

        HPBC_POSTCONDITION2(isValid(result));
        return result;
    }

    using BC::add;
    HURCHALLA_FORCE_INLINE V add(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));        // we know  -n <= x.get() < n
        HPBC_PRECONDITION2(isValid(y));
        C cy = getCanonicalValue(y);
        T tx = static_cast<T>(x.get());
        T tmpx = static_cast<T>(tx - n_);
        // Note if x.get() is negative, then tx >= n, since we know n < R/2.
        // This is assuming wrap-around conversion to T.  (ie two's complement?)
        // Likewise if tx >= n, then x.get() can't be >=0 because x.get() >= 0
        // would mean conversion to T would leave the value unchanged after
        // conversion: thus it would result in tx < n, which is obviously
        // impossible if tx >= n.  This means
        // that if tx >= n, then x.get() < 0.
        // And by contrapositive of the first item, tx < n implies x.get() >= 0

        HURCHALLA_CMOV(tx >= n_, tmpx, tx);    // if tx>=n_, set tmpx=tx
        // Assume tx >= n.  Then x.get() < 0, and the cmov will have set
        //    tmpx = tx = static_cast<T>(x.get()).
        //    And hence,  static_cast<S>(tmpx) == x.get() < 0.  We also know
        //    from isValid(x) that  -static_cast<S>(n) <= x.get().  Thus,
        //    -static_cast<S>(n) <= static_cast<S>(tmpx) < 0.
        // Assume tx < n.  Then x.get() >= 0, and the cmov will have left
        //    tmpx unchanged, so we know
        //    tmpx == tx - n == static_cast<T>(x.get()) - n,  and so
        //    static_cast<S>(tmpx) == x.get() - static_cast<S>(n).
        //    Since x.get() >= 0,  static_cast<S>(tmpx) >= -static_cast<S>(n).
        //    We know from isValid(x) that  x.get() < static_cast<S>(n),  and so
        //    static_cast<S>(tmpx) == x.get() - static_cast<S>(n) < n - n == 0.
        //    Putting this all together, we have
        //    -static_cast<S>(n) <= static_cast<S>(tmpx) < 0.
        // For both assumptions, we know we get
        // -static_cast<S>(n) <= static_cast<S>(tmpx) < 0.
        HPBC_ASSERT2(-static_cast<S>(n_) <= static_cast<S>(tmpx)
                     && static_cast<S>(tmpx) < 0);
        HPBC_ASSERT2(0 <= cy.get() && cy.get() < n_); // because cy is canonical
        S result = static_cast<S>(
                               static_cast<S>(tmpx) + static_cast<S>(cy.get()));
        // we can see result will satisfy  -n <= result < n - 1
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V add(V x, C cy) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(0 <= cy.get() && cy.get() < n_); // cy is canonical
        T tx = static_cast<T>(x.get());
        T tmpx = static_cast<T>(tx - n_);
        HURCHALLA_CMOV(tx >= n_, tmpx, tx);    // if tx>=n_, set tmpx=tx
        // For proof that the following assert is true, see  add(V, V) above.
        HPBC_ASSERT2(-static_cast<S>(n_) <= static_cast<S>(tmpx)
                     && static_cast<S>(tmpx) < 0);
        S result = static_cast<S>(
                               static_cast<S>(tmpx) + static_cast<S>(cy.get()));
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }

    using BC::subtract;
    HURCHALLA_FORCE_INLINE V subtract(V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        C cx = getCanonicalValue(x);
        C cy = getCanonicalValue(y);
        S result = static_cast<S>(
                           static_cast<S>(cx.get()) - static_cast<S>(cy.get()));
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V subtract(V x, C cy) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(cy.get() < n_);
        C cx = getCanonicalValue(x);
        S result = static_cast<S>(
                           static_cast<S>(cx.get()) - static_cast<S>(cy.get()));
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }
    HURCHALLA_FORCE_INLINE V subtract(C cx, V y) const
    {
        HPBC_PRECONDITION2(cx.get() < n_);
        HPBC_PRECONDITION2(isValid(y));
        C cy = getCanonicalValue(y);
        S result = static_cast<S>(
                           static_cast<S>(cx.get()) - static_cast<S>(cy.get()));
        HPBC_POSTCONDITION2(isValid(V(result)));
        return V(result);
    }

    template <typename J, typename K>
    HURCHALLA_FORCE_INLINE V unordered_subtract(J x, K y) const
    {
        return subtract(x, y);
    }

private:
    // functions called by the 'curiously recurring template pattern' base (BC).
    friend BC;

    // note: PTAG gets ignored here
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(bool& resultIsZero, T u_hi, T u_lo, PTAG) const
    {
        HPBC_PRECONDITION2(u_hi < n_);  // verifies that (u_hi*R + u_lo) < n*R
        bool isNegative;  // ignored
        T result = REDC_incomplete(isNegative, u_hi, u_lo, n_, BC::inv_n_);
        resultIsZero = (result == 0);
        V v = V(static_cast<S>(result));
        HPBC_POSTCONDITION2(isValid(v));
        return v;
    }
    template <class PTAG> HURCHALLA_FORCE_INLINE
    V montyREDC(T u_hi, T u_lo, PTAG) const
    {
        bool resultIsZero;  // ignored
        return montyREDC(resultIsZero, u_hi, u_lo, PTAG());
    }

    // Let u be an arbitrary double-word value that is congruent (mod n_) to the
    // product of x and y, and that satisfies 0 <= u < n_*R.  Return the high
    // word of u, and write the low word of u to u_lo.
    HURCHALLA_FORCE_INLINE T multiplyToHiLo(T& u_lo, V x, V y) const
    {
        HPBC_PRECONDITION2(isValid(x));
        HPBC_PRECONDITION2(isValid(y));
        S product_hi = signed_multiply_to_hilo_product(u_lo, x.get(), y.get());
        // By isValid(), x and y both have range [-n, n).  Thus x*y has range
        // (-n*n, n*n], and due to class invariant n<R/2, x*y has range
        // (-n*R/2, n*R/2).  Since  product_hi*R + u_lo == x*y,
        // -n*R/2 < product_hi*R + u_lo < n*R/2.  Thus
        // product_hi*R < n*R/2 - u_lo <= n*R/2 < n*R/2 + R/2 == (n+1)*R/2
        // Since product_hi*R < ((n+1)/2)*R  and  n is odd, 
        // product_hi < (n + 1)/2 == (n-1)/2 + 1.  Thus product_hi <= (n-1)/2.
        // Likewise, since  product_hi*R + u_lo > -n*R/2,
        // product_hi*R > -n*R/2 - u_lo > -n*R/2 - R > -n*R/2 - 3*R/2,  and so
        // product_hi*R > -((n+3)/2)*R.  Therefore
        // product_hi > -((n+3)/2) == -((n+1)/2) - 1.  And thus
        // product_hi >= -(n+1)/2.  Putting it all together,
        // -(n+1)/2 <= product_hi <= (n-1)/2.  Also, due to invariant  n < R/2,
        // -R/2 < -(n+1)/2 <= product_hi <= (n-1)/2 < R/2.
        T tmp = static_cast<T>(product_hi);
        T u_hi = static_cast<T>(tmp + n_);
        // We can easily deduce from the results in the above comments that
        // (n-1)/2 <= u_hi <= n + (n-1)/2 < R.

// TODO: possibly offer predefined macro for using _asm for this cmov.  Note
// that this function is always or nearly always used for a montgomery multiply,
// and so we don't need to lower latency for our result because the first two
// mults of REDC will run in parallel.  But _asm might reduce uops by 1.
        // If (u_hi < n) is true, it indicates product_hi < 0.  Likewise,
        // (u_hi >= n) true would indicate product_hi >= 0.
        HURCHALLA_CMOV(u_hi >= n_, u_hi, tmp); // if u_hi>=n_, set u_hi = tmp
        HPBC_POSTCONDITION2(0 <= u_hi && u_hi < n_);
        return u_hi;
    }
    HURCHALLA_FORCE_INLINE T squareToHiLo(T& u_lo, V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        S tmp_hi = signed_multiply_to_hilo_product(u_lo, x.get(), x.get());
        // The same logic as given in multiplyToHiLo shows that
        // -(n+1)/2 <= tmp_hi <= (n-1)/2.  But additionally, since the square
        // of an integer is always >= 0, we therefore know
        // 0 <= tmp_hi <= (n-1)/2.
        HPBC_ASSERT2(tmp_hi >= 0);
        T u_hi = static_cast<T>(tmp_hi);
        HPBC_POSTCONDITION2(0 <= u_hi && u_hi < n_);
        return u_hi;
    }
    HURCHALLA_FORCE_INLINE bool isValid(V x) const
    {
        static_assert(!(ut_numeric_limits<T>::is_signed), "");
        T Rdiv2 = static_cast<T>(static_cast<T>(1) <<
                                       (ut_numeric_limits<T>::digits - 1));
        HPBC_INVARIANT2(n_ < Rdiv2);
        return (-static_cast<S>(n_) <= x.get() && x.get() < static_cast<S>(n_));
    }
    HURCHALLA_FORCE_INLINE bool isCanonical(V x) const
    {
        static_assert(std::is_same<S, decltype(x.get())>::value, "");
        return (0 <= x.get() && x.get() < static_cast<S>(n_));
    }
    // get a Natural number (i.e. number >= 0) congruent to x (mod n)
    HURCHALLA_FORCE_INLINE T getNaturalEquivalence(V x) const
    {
        HPBC_PRECONDITION2(isValid(x));
        T result = static_cast<T>(static_cast<T>(x.get()) + n_);
        // Since the precondition of isValid(x) required
        // -(static_cast<S>(n)) <= x.get() < n,  the sum for 'result' will have
        // carried (wrapped-around) if  x.get() < 0.
        HPBC_POSTCONDITION2(0 <= result && result < static_cast<T>(2*n_));
        return result;
    }
};


}} // end namespace

#endif


/*
    //   For types T larger than the CPU integer register size, square() is
    // likely to perform worse or the same as multiply(x, x).  The reason is
    // that this version of square() uses a signed multiply instead of an
    // unsigned multiply, and signed multiply is less efficient than unsigned
    // when it's a function call rather than a single native CPU instruction.
    // Types T that are larger than the register size can't be multiplied with a
    // single native CPU instruction.
*/
