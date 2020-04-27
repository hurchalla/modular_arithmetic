
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_UNSIGNED_MULT_TO_HILO_H__INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_IMPL_UNSIGNED_MULT_TO_HILO_H__INCLUDED


#include "hurchalla/programming_by_contract/programming_by_contract.h"
#include <cstdint>
#include <limits>

namespace hurchalla { namespace montgomery_arithmetic {


/* Generic (non-platform specific) implementation of the contract for
   T unsigned_multiply_to_hilo_product(T* pLowProduct, T u, T v).
Return Value:
   Returns the high portion of the product.
Notes:
   I adapted this code from https://stackoverflow.com/a/58381061
   On ARM32 with clang it compiles nicely, using the UMAAL instruction.
*/
#ifndef COMPILE_ERROR_ON_SLOW_MATH
template <typename T>
T impl_unsigned_multiply_to_hilo_product(T* pLowProduct, T u, T v)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");

    static const unsigned int shift = sizeof(T)*(8/2);
    // for example, if T==uint64_t, shift ought to == 32
    static const T lowmask = static_cast<T>(static_cast<T>(1) << shift) - 1;
    // for example, if T==uint64_t, lowmask ought to == 0xFFFFFFFF

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
#else
// cause a compile error if instantiating the (slow) template function
template <typename T>
T impl_unsigned_multiply_to_hilo_product(T* pLowProduct, T u, T v) = delete;
#endif // #ifndef COMPILE_ERROR_ON_SLOW_MATH




// Intended for use by the functions below
template <typename T, typename T2>
T umult_to_hilo_product(T* pLowProduct, T u, T v)
{
    static_assert(std::numeric_limits<T>::is_integer, "");
    static_assert(!(std::numeric_limits<T>::is_signed), "");
    static_assert(std::numeric_limits<T2>::is_integer, "");
    static_assert(!(std::numeric_limits<T2>::is_signed), "");
    static_assert(std::numeric_limits<T2>::digits >=
                  2*std::numeric_limits<T>::digits, "");
    T2 product = (T2)u * (T2)v;
    *pLowProduct = (T)(product);
    return (T)(product >> std::numeric_limits<T>::digits);
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
inline uint8_t impl_unsigned_multiply_to_hilo_product(uint8_t* pLowProduct,
                                                       uint8_t u, uint8_t v)
{
    // Note we could have used 'T2 = unsigned int' since 'int' is >= 16bit.
    using T2 = uint16_t;
    return umult_to_hilo_product<decltype(u),T2>(pLowProduct, u, v);
}
inline uint16_t impl_unsigned_multiply_to_hilo_product(uint16_t* pLowProduct,
                                                       uint16_t u, uint16_t v)
{
    using T2 = uint32_t;
    return umult_to_hilo_product<decltype(u),T2>(pLowProduct, u, v);
}
inline uint32_t impl_unsigned_multiply_to_hilo_product(uint32_t* pLowProduct,
                                                       uint32_t u, uint32_t v)
{
    using T2 = uint64_t;
    return umult_to_hilo_product<decltype(u),T2>(pLowProduct, u, v);
}


// These fast 'uint64_t' functions use intrinsics (MSVC), or compiler
// extensions (__uint128_t on GNU compatible compilers).
// Assembly versions for x86 or ARM aren't needed - clang/gcc/icc generate
// assembly that is good enough via __uint128_t, and MSVC does well using
// the intrinsics.

// MSVC + x64
#if defined(_MSC_VER) && defined(TARGET_ISA_X86_64)
#  include <intrin.h>
inline uint64_t impl_unsigned_multiply_to_hilo_product(uint64_t* pLowProduct,
                                                       uint64_t u, uint64_t v)
{
    uint64_t highProduct;
    *pLowProduct = _umul128(u, v, &highProduct);
    if (POSTCONDITION3_MACRO_IS_ACTIVE) {
        uint64_t tmpHi, tmpLo;
        tmpHi = impl_unsigned_multiply_to_hilo_product<uint64_t>(&tmpLo, u, v);
        postcondition3(highProduct == tmpHi && *pLowProduct == tmpLo);
    }
    return highProduct;
}

// MSVC + ARM64
#elif defined(_MSC_VER) && defined(TARGET_ISA_ARM_64)
#  include <intrin.h>
inline uint64_t impl_unsigned_multiply_to_hilo_product(uint64_t* pLowProduct,
                                                       uint64_t u, uint64_t v)
{
    uint64_t highProduct = __umulh(u, v);
    *pLowProduct = u*v;
    if (POSTCONDITION3_MACRO_IS_ACTIVE) {
        uint64_t tmpHi, tmpLo;
        tmpHi = impl_unsigned_multiply_to_hilo_product<uint64_t>(&tmpLo, u, v);
        postcondition3(highProduct == tmpHi && *pLowProduct == tmpLo);
    }
    return highProduct;
}

// GNU compatible compiler (for __uint128_t support)
#elif TARGET_BIT_WIDTH >= 64 && ( defined(__SIZEOF_INT128__) || \
                  (__clang_major__ >= 3) || (__INTEL_COMPILER >= 1300) || \
                  ( !defined(__INTEL_COMPILER) && !defined(__clang_major__) && \
                    (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1)) ) )
// The macro __SIZEOF_INT128__ indicates if __int128 is supported, though some
// compilers supported __int128 prior to providing that macro.  Clang had it at
// least since 3.0 and icc got it in 13.0.  Gcc has had since at least 4.1.  See
// https://stackoverflow.com/questions/16088282/is-there-a-128-bit-integer-in-gcc
// Note also that clang and icc define the GNUC and GNUC_MINOR macros.  Clang
// v3.0 defines them as 4,2.  Icc v13 defines them as 4,7.
inline uint64_t impl_unsigned_multiply_to_hilo_product(uint64_t* pLowProduct,
                                                       uint64_t u, uint64_t v)
{
    using T2 = __uint128_t;
    return umult_to_hilo_product<decltype(u),T2>(pLowProduct, u, v);
}
#endif




// Note for MSVC: 'uint32_t' functions using intrinsics don't improve the asm
// generated compared to the simple function implementation, and so intrinsic
// versions are not present here.  For reference, the intrinsics would have
// been __emulu (for X86) and _arm_umull (for ARM32).


// Note for MSVC for X86-32 bit: the asm generated for the 'uint64_t' version
// via the generic template function is terrible.  In contrast, the code in
// comments below produces decent asm for MSVC x86-32.  I expect this code to
// be fine if you wish to use it.  However, it remains in comments because
// x86-32 is obsolete (making it an unnecessary complication).
/*
// MSVC + x32
#if defined(_MSC_VER) && defined(TARGET_ISA_X86_32)
inline uint64_t impl_unsigned_multiply_to_hilo_product(uint64_t* pLowProduct,
                                                       uint64_t u, uint64_t v)
{
    // this code is a modification of the generic template function, performing
    // the arithmetic with uint32_t variables rather than uint64_t variables.

    uint32_t u0 = (uint32_t)u;
    uint32_t v0 = (uint32_t)v;
    uint32_t u1 = (uint32_t)(u >> 32);
    uint32_t v1 = (uint32_t)(v >> 32);

    uint32_t lolo_lo;
    uint32_t lolo_hi = impl_unsigned_multiply_to_hilo_product(&lolo_lo, u0, v0);
    uint32_t hilo_lo;
    uint32_t hilo_hi = impl_unsigned_multiply_to_hilo_product(&hilo_lo, u1, v0);
    uint32_t lohi_lo;
    uint32_t lohi_hi = impl_unsigned_multiply_to_hilo_product(&lohi_lo, u0, v1);
    uint32_t hihi_lo;
    uint32_t hihi_hi = impl_unsigned_multiply_to_hilo_product(&hihi_lo, u1, v1);

    uint64_t cross = (uint64_t)lolo_hi + (uint64_t)hilo_lo +
                     (((uint64_t)lohi_hi << 32) | (uint64_t)lohi_lo);
    uint64_t high = (uint64_t)hilo_hi + (cross >> 32) +
                    (((uint64_t)hihi_hi << 32) | (uint64_t)hihi_lo);

    *pLowProduct = (cross << 32) | (uint64_t)(lolo_lo);
    return high;
}
#endif
*/


}} // end namespace

#endif
