// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_MONTY_COMMON_H_INCLUDED


#include "hurchalla/montgomery_arithmetic/detail/platform_specific/RedcLargeR.h"
#include "hurchalla/montgomery_arithmetic/detail/RedcSmallR.h"
#include "hurchalla/montgomery_arithmetic/detail/platform_specific/MontHelper.h"
#include "hurchalla/montgomery_arithmetic/detail/unsigned_multiply_to_hilo_product.h"
#include "hurchalla/montgomery_arithmetic/detail/safely_promote_unsigned.h"
#include "hurchalla/montgomery_arithmetic/detail/sized_uint.h"
#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {

// The public functions are at the bottom of this file.

namespace detail_monty_common {
// -----------------
// Private Functions
// -----------------


// default (primary) template
template <typename T, class Enable = void>
struct MontFunctionsCommon
{
  static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
  static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
  static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");

  template <class MTAG, class PTAG>
  static HURCHALLA_FORCE_INLINE
  T mul(T x, T y, T n, T neg_inv_n, MTAG, PTAG)
  {
    T u_lo;
    T u_hi = unsigned_multiply_to_hilo_product(&u_lo, x, y);
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  u = x*y < n*R.
    // Having u_hi < n is sufficient and necessary to satisfy our requirement of
    // x*y == u < n*R.  See REDC_non_minimized() in RedcLargeR.h for proof.
    HPBC_PRECONDITION2(u_hi < n);

    T result = RedcLargeR<T>::REDC(u_hi, u_lo, n, neg_inv_n, MTAG(), PTAG());
    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
  }

  template <class MTAG, class PTAG>
  static HURCHALLA_FORCE_INLINE
  T fmsub(T x, T y, T z, T n, T neg_inv_n, MTAG, PTAG)
  {
    HPBC_PRECONDITION2(z < n);  // z must be canonical (0 <= z < n)
    T u_lo;
    T u_hi = unsigned_multiply_to_hilo_product(&u_lo, x, y);
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  u = x*y < n*R.
    // Having u_hi < n is sufficient and necessary to satisfy our requirement of
    // x*y == u < n*R.  See REDC_non_minimized() in RedcLargeR.h for proof.
    HPBC_PRECONDITION2(u_hi < n);

    // TODO proof of correctness, showing that performing the modular sub prior
    // to the REDC will always give the same results as performing the REDC and
    // then the modular subtraction.
    // The following calculations should execute in parallel with the first two
    // multiplies in REDC_non_minimized(), since those mutiplies do not depend
    // on these calculations.  (Instruction level parallelism)
    T diff = MontHelper<T>::modsub_canonical_subtrahend(u_hi, z, n);
    // modsub_canonical_subtrahend()'s postcondition guarantees diff is a valid
    // montgomery value for any MTAG
    T result = RedcLargeR<T>::REDC(diff, u_lo, n, neg_inv_n, MTAG(), PTAG());

    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
  }

  template <class MTAG, class PTAG>
  static HURCHALLA_FORCE_INLINE
  T fmadd(T x, T y, T z, T n, T neg_inv_n, MTAG, PTAG)
  {
    HPBC_PRECONDITION2(z < n);  // z must be canonical (0 <= z < n)
    T u_lo;
    T u_hi = unsigned_multiply_to_hilo_product(&u_lo, x, y);
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  u = x*y < n*R.
    // Having u_hi < n is sufficient and necessary to satisfy our requirement of
    // x*y == u < n*R.  See REDC_non_minimized() in RedcLargeR.h for proof.
    HPBC_PRECONDITION2(u_hi < n);

    // TODO proof of correctness, showing that performing the modular add prior
    // to the REDC will always give the same results as performing the REDC and
    // then the modular addition.
    // The following calculations should execute in parallel with the first two
    // multiplies in REDC_non_minimized(), since those mutiplies do not depend
    // on these calculations.  (Instruction level parallelism)
    T sum = MontHelper<T>::modadd_canonical_second_addend(u_hi, z, n);
    // modadd_canonical_second_addend()'s postcondition guarantees sum is a
    // valid montgomery value for any MTAG.
    T result = RedcLargeR<T>::REDC(sum, u_lo, n, neg_inv_n, MTAG(), PTAG());

    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
  }

  template <class MTAG>
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T neg_inv_n, MTAG)
  {
    T result = RedcLargeR<T>::convert_out(x, n, neg_inv_n, MTAG());
    HPBC_POSTCONDITION2(result < n);
    return result;
  }
};



#ifndef HURCHALLA_TARGET_BIT_WIDTH
#error "HURCHALLA_TARGET_BIT_WIDTH must be defined"
#endif

// specialization for types T that are smaller than HURCHALLA_TARGET_BIT_WIDTH
template <typename T>
struct MontFunctionsCommon< T, typename std::enable_if<
            modular_arithmetic::ma_numeric_limits<T>::is_specialized &&
            HURCHALLA_TARGET_BIT_WIDTH >=
                            2 * modular_arithmetic::ma_numeric_limits<T>::digits
            >::type >
{
  using T2 = typename sized_uint<
                    2 * modular_arithmetic::ma_numeric_limits<T>::digits>::type;
  static_assert(modular_arithmetic::ma_numeric_limits<T>::is_integer, "");
  static_assert(!(modular_arithmetic::ma_numeric_limits<T>::is_signed), "");
  static_assert(modular_arithmetic::ma_numeric_limits<T>::is_modulo, "");
  static_assert(modular_arithmetic::ma_numeric_limits<T2>::is_integer, "");
  static_assert(!(modular_arithmetic::ma_numeric_limits<T2>::is_signed), "");
  static_assert(modular_arithmetic::ma_numeric_limits<T2>::is_modulo, "");

  template <class MTAG, class PTAG>
  static HURCHALLA_FORCE_INLINE
  T mul(T x, T y, T n, T neg_inv_n, MTAG, PTAG)
  {
    T2 u = static_cast<T2>(static_cast<T2>(x) * static_cast<T2>(y));
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  u = x*y < n*R, or equivalently:
    HPBC_PRECONDITION2(u < (static_cast<T2>(n)
                            << std::numeric_limits<T>::digits));

    T u_lo = static_cast<T>(u);
    T result = RedcSmallR<T>::REDC(u, u_lo, n, neg_inv_n, MTAG(), PTAG());
    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
  }

  template <class MTAG, class PTAG>
  static HURCHALLA_FORCE_INLINE
  T fmsub(T x, T y, T z, T n, T neg_inv_n, MTAG, PTAG)
  {
    HPBC_PRECONDITION2(z < n);  // z must be canonical (0 <= z < n)
    T2 u = static_cast<T2>(static_cast<T2>(x) * static_cast<T2>(y));
    static constexpr int bit_width_T = std::numeric_limits<T>::digits;
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  u = x*y < n*R, or equivalently:
    HPBC_PRECONDITION2(u < (static_cast<T2>(n) << bit_width_T));

    // TODO proof of correctness, showing that performing the modular sub prior
    // to the REDC will always give the same results as performing the REDC and
    // then the modular subtraction.
    // The following calculations should execute in parallel with the first two
    // multiplies in REDC_non_minimized2(), since those mutiplies do not depend
    // on these calculations.  (Instruction level parallelism)
    T2 zR = static_cast<T2>(static_cast<T2>(z) << bit_width_T);
    T2 nR = static_cast<T2>(static_cast<T2>(n) << bit_width_T);
    T2 u2 = MontHelper<T2>::modsub_canonical_subtrahend(u, zR, nR);
    // the low bits should be unchanged between u2 and u.
    HPBC_ASSERT2(static_cast<T>(u2) == static_cast<T>(u));

    T u_lo = static_cast<T>(u);
    T result = RedcSmallR<T>::REDC(u2, u_lo, n, neg_inv_n, MTAG(), PTAG());

    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
  }

  template <class MTAG, class PTAG>
  static HURCHALLA_FORCE_INLINE
  T fmadd(T x, T y, T z, T n, T neg_inv_n, MTAG, PTAG)
  {
    HPBC_PRECONDITION2(z < n);  // z must be canonical (0 <= z < n)
    T2 u = static_cast<T2>(static_cast<T2>(x) * static_cast<T2>(y));
    static constexpr int bit_width_T = std::numeric_limits<T>::digits;
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  u = x*y < n*R, or equivalently:
    HPBC_PRECONDITION2(u < (static_cast<T2>(n) << bit_width_T));

    // TODO proof of correctness, showing that performing the modular add prior
    // to the REDC will always give the same results as performing the REDC and
    // then the modular addition.
    // The following calculations should execute in parallel with the first two
    // multiplies in REDC_non_minimized2(), since those mutiplies do not depend
    // on these calculations.  (Instruction level parallelism)
    T2 zR = static_cast<T2>(static_cast<T2>(z) << bit_width_T);
    T2 nR = static_cast<T2>(static_cast<T2>(n) << bit_width_T);
    T2 u2 = MontHelper<T2>::modadd_canonical_second_addend(u, zR, nR);
    // the low bits should be unchanged between u2 and u.
    HPBC_ASSERT2(static_cast<T>(u2) == static_cast<T>(u));

    T u_lo = static_cast<T>(u);
    T result = RedcSmallR<T>::REDC(u2, u_lo, n, neg_inv_n, MTAG(), PTAG());

    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
  }

  template <class MTAG>
  static HURCHALLA_FORCE_INLINE
  T convert_out(T x, T n, T neg_inv_n, MTAG)
  {
    T result = RedcSmallR<T>::convert_out(x, n, neg_inv_n, MTAG());
    HPBC_POSTCONDITION2(result < n);
    return result;
  }
};


} // end namespace detail_monty_common


// ----------------
// Public Functions
// ----------------


// Multiplies two mongomery values x and y.
// Returns the product as a montgomery value.
template <typename T, class MTAG, class PTAG>
HURCHALLA_FORCE_INLINE
T montmul(T x, T y, T n, T neg_inv_n, MTAG, PTAG)
{
    // Precondition: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  x*y < n*R.  (The constant R
    // represents the value R = 2^(ma::ma_numeric_limits<T>::digits))
    T result = detail_monty_common::MontFunctionsCommon<T>::
                                        mul(x, y, n, neg_inv_n, MTAG(), PTAG());
    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
}

// Multiplies two mongomery values x and y, and then subtracts montgomery
// value z from the product.  Returns the resulting montgomery value.
template <typename T, class MTAG, class PTAG>
HURCHALLA_FORCE_INLINE
T montfmsub(T x, T y, T z, T n, T neg_inv_n, MTAG, PTAG)
{
    // Precondition #1: z must be canonical (i.e. 0 <= z < n).
    HPBC_PRECONDITION2(z < n);
    // Precondition #2: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  x*y < n*R.  (The constant R
    // represents the value R = 2^(ma::ma_numeric_limits<T>::digits))

    T result = detail_monty_common::MontFunctionsCommon<T>::
                                   fmsub(x, y, z, n, neg_inv_n, MTAG(), PTAG());
    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
}

// Multiplies two mongomery values x and y, and then adds montgomery
// value z to the product.  Returns the resulting montgomery value.
template <typename T, class MTAG, class PTAG>
HURCHALLA_FORCE_INLINE
T montfmadd(T x, T y, T z, T n, T neg_inv_n, MTAG, PTAG)
{
    // Precondition #1: z must be canonical (i.e. 0 <= z < n).
    HPBC_PRECONDITION2(z < n);
    // Precondition #2: Assuming theoretical unlimited precision standard
    // multiplication, this function requires  x*y < n*R.  (The constant R
    // represents the value R = 2^(ma::ma_numeric_limits<T>::digits))

    T result = detail_monty_common::MontFunctionsCommon<T>::
                                   fmadd(x, y, z, n, neg_inv_n, MTAG(), PTAG());
    // Postcondition: result is a valid montgomery value for the montgomery type
    // associated with MTAG
    return result;
}

template <typename T, class MTAG>
HURCHALLA_FORCE_INLINE
T montout(T x, T n, T neg_inv_n, MTAG)
{
    T result = detail_monty_common::MontFunctionsCommon<T>::
                                           convert_out(x, n, neg_inv_n, MTAG());
    HPBC_POSTCONDITION2(result < n);
    return result;
}


}} // end namespace

#endif
