// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_UNSIGNED_MULT_TO_HILO_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_UNSIGNED_MULT_TO_HILO_H_INCLUDED


#include "hurchalla/modular_arithmetic/detail/ma_numeric_limits.h"
#include "hurchalla/modular_arithmetic/detail/platform_specific/compiler_macros.h"
#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#if defined(_MSC_VER)
#  include <intrin.h>
#  pragma warning(push)
#  pragma warning(disable : 4127)
#endif

namespace hurchalla { namespace montgomery_arithmetic {


// Generic (non-platform specific) implementation of the contract for
//   T unsigned_multiply_to_hilo_product(T* pLowProduct, T u, T v).
// Return Value:
//   Returns the high portion of the product.
// Notes:
//   I adapted this code from https://stackoverflow.com/a/58381061
//   On ARM32 with clang it compiles nicely, using the UMAAL instruction.
template <typename T>
HURCHALLA_FORCE_INLINE
T slow_unsigned_multiply_to_hilo_product(T* pLowProduct, T u, T v)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T>::is_signed), "");

    // for example, if T==uint64_t, shift ought to == 32
    static const unsigned int shift = ma::ma_numeric_limits<T>::digits / 2;
    // for example, if T==uint64_t, lowmask ought to == 0xFFFFFFFF
    static const T lowmask = (static_cast<T>(1) << shift) - static_cast<T>(1);

    T u0 = u & lowmask;
    T v0 = v & lowmask;
    T u1 = u >> shift;
    T v1 = v >> shift;

    // Calculate all the cross products.
    T lo_lo = u0 * v0;
    T hi_lo = u1 * v0;
    T lo_hi = u0 * v1;
    T hi_hi = u1 * v1;

    // The next statement will not overflow.  Proof: let S=2^(shift). We can see
    // that both (lo_lo >> shift) and (hi_lo & lowmask) must be less than S.
    // Therefore the max possible value of cross= (S-1) + (S-1) + (S-1)*(S-1) ==
    // S-1 + S-1 + S*S - 2*S + 1 == S*S - 1, which is the max value that can be
    // represented in type T.  Thus the calculation will never overflow.
    T cross = (lo_lo >> shift) + (hi_lo & lowmask) + lo_hi;
    // The next statement will not overflow, for the same reason as above.
    T high = (hi_lo >> shift) + (cross >> shift) + hi_hi;

    *pLowProduct = (cross << shift) | (lo_lo & lowmask);
    return high;
}


template <typename T>
#ifdef HURCHALLA_COMPILE_ERROR_ON_SLOW_MATH
  // cause a compile error instead of falling back to the slow template function
  T impl_unsigned_multiply_to_hilo_product(T* pLowProduct, T u, T v) = delete;
#else
  HURCHALLA_FORCE_INLINE
  T impl_unsigned_multiply_to_hilo_product(T* pLowProduct, T u,T v)
  {
      return slow_unsigned_multiply_to_hilo_product(pLowProduct, u, v);
  }
#endif




// Intended for use by the functions below
template <typename T, typename T2>
HURCHALLA_FORCE_INLINE T umult_to_hilo_product(T* pLowProduct, T u, T v)
{
    namespace ma = hurchalla::modular_arithmetic;
    static_assert(ma::ma_numeric_limits<T>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T>::is_signed), "");
    static_assert(ma::ma_numeric_limits<T2>::is_integer, "");
    static_assert(!(ma::ma_numeric_limits<T2>::is_signed), "");
    static_assert(ma::ma_numeric_limits<T2>::digits >=
                  2*ma::ma_numeric_limits<T>::digits, "");
    T2 product = static_cast<T2>(static_cast<T2>(u) * static_cast<T2>(v));
    *pLowProduct = static_cast<T>(product);
    return  static_cast<T>(product >> ma::ma_numeric_limits<T>::digits);
}



// -------- PLATFORM SPECIFIC nontemplate overloads ----------

// Note: these fast nontemplate function overloads get first priority for being
// called (see http://www.gotw.ca/publications/mill17.htm ), when both one of
// these nontemplate functions and the generic template function match the
// caller's provided argument type(s).



// Note that when using these simple functions, the generated asm from
// clang/icc/gcc is generally quite good.
// GCC for ARM seems to make the worst generated asm, but it's not so bad as to
// make inline asm seem worthwhile.
HURCHALLA_FORCE_INLINE std::uint8_t impl_unsigned_multiply_to_hilo_product(
                      std::uint8_t* pLowProduct, std::uint8_t u, std::uint8_t v)
{
    // Note we could have used 'T2 = unsigned int' since 'int' is >= 16bit.
    using T2 = std::uint16_t;
    return umult_to_hilo_product<decltype(u),T2>(pLowProduct, u, v);
}
HURCHALLA_FORCE_INLINE std::uint16_t impl_unsigned_multiply_to_hilo_product(
                   std::uint16_t* pLowProduct, std::uint16_t u, std::uint16_t v)
{
    using T2 = std::uint32_t;
    return umult_to_hilo_product<decltype(u),T2>(pLowProduct, u, v);
}
// --------------------------------------------------------------------------
HURCHALLA_FORCE_INLINE std::uint32_t impl_unsigned_multiply_to_hilo_product(
                   std::uint32_t* pLowProduct, std::uint32_t u, std::uint32_t v)
{
    using T2 = std::uint64_t;
    return umult_to_hilo_product<decltype(u),T2>(pLowProduct, u, v);
}
// Note for MSVC: 'uint32_t' functions using intrinsics don't improve the asm
// generated compared to the simple function implementation, and so intrinsic
// versions are not present here.  For reference, the intrinsics would have
// been __emulu (for X86) and _arm_umull (for ARM32).
// --------------------------------------------------------------------------


// The following fast 'uint64_t' functions use intrinsics (MSVC), or compiler
// extensions (__uint128_t on GNU compatible compilers).
// Assembly versions for x86 or ARM aren't needed - clang/gcc/icc generate
// assembly that is good enough via __uint128_t, and MSVC does well using
// the intrinsics.

// MSVC + x64
#if defined(_MSC_VER) && defined(HURCHALLA_TARGET_ISA_X86_64)
HURCHALLA_FORCE_INLINE std::uint64_t impl_unsigned_multiply_to_hilo_product(
                   std::uint64_t* pLowProduct, std::uint64_t u, std::uint64_t v)
{
    using std::uint64_t;
    uint64_t highProduct;
    *pLowProduct = _umul128(u, v, &highProduct);
    if (HPBC_POSTCONDITION3_MACRO_IS_ACTIVE) {
        uint64_t tmpHi, tmpLo;
        tmpHi = slow_unsigned_multiply_to_hilo_product(&tmpLo, u, v);
        HPBC_POSTCONDITION3(highProduct == tmpHi && *pLowProduct == tmpLo);
    }
    return highProduct;
}

// MSVC + ARM64
#elif defined(_MSC_VER) && defined(_M_ARM64)
HURCHALLA_FORCE_INLINE std::uint64_t impl_unsigned_multiply_to_hilo_product(
                   std::uint64_t* pLowProduct, std::uint64_t u, std::uint64_t v)
{
    using std::uint64_t;
    uint64_t highProduct = __umulh(u, v);
    *pLowProduct = u*v;
    if (HPBC_POSTCONDITION3_MACRO_IS_ACTIVE) {
        uint64_t tmpHi, tmpLo;
        tmpHi = slow_unsigned_multiply_to_hilo_product(&tmpLo, u, v);
        HPBC_POSTCONDITION3(highProduct == tmpHi && *pLowProduct == tmpLo);
    }
    return highProduct;
}

#elif (HURCHALLA_COMPILER_HAS_UINT128_T())
HURCHALLA_FORCE_INLINE std::uint64_t impl_unsigned_multiply_to_hilo_product(
                   std::uint64_t* pLowProduct, std::uint64_t u, std::uint64_t v)
{
    using T2 = __uint128_t;
    return umult_to_hilo_product<decltype(u),T2>(pLowProduct, u, v);
}
#endif


}} // end namespace


#if defined(_MSC_VER)
#  pragma warning(pop)
#endif

#endif
