// --- This file is distributed under the MIT Open Source License, as detailed
// by the file "LICENSE.TXT" in the root of this repository ---

#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_COMPILER_MACROS_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_COMPILER_MACROS_H_INCLUDED


#ifdef _MSC_VER
#  define HURCHALLA_FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
   // Note that Intel compiler (icc) also accepts __forceinline
#  define HURCHALLA_FORCE_INLINE inline __attribute__((always_inline))
#else
#  define HURCHALLA_FORCE_INLINE inline
#endif


#if defined(__clang__) || defined(__INTEL_COMPILER)
#  define HURCHALLA_REQUEST_UNROLL_LOOP _Pragma("unroll")
#elif defined(__GNUC__) && __GNUC__ >= 8
#  define HURCHALLA_REQUEST_UNROLL_LOOP _Pragma("GCC unroll 127")
#else
#  define HURCHALLA_REQUEST_UNROLL_LOOP
#endif


#if defined(__x86_64__) || defined(_M_X64)
#  ifndef HURCHALLA_TARGET_ISA_X86_64
#    define HURCHALLA_TARGET_ISA_X86_64 1
#  endif
#elif defined(__i386) || defined(_M_IX86)
#  ifndef HURCHALLA_TARGET_ISA_X86_32
#    define HURCHALLA_TARGET_ISA_X86_32 1
#  endif
#endif

#if defined(__aarch64__) || defined(_M_ARM64)
#  ifndef HURCHALLA_TARGET_ISA_ARM_64
#    define HURCHALLA_TARGET_ISA_ARM_64 1
#  endif
#elif defined(__arm__) || defined(_M_ARM)
#  ifndef HURCHALLA_TARGET_ISA_ARM_32
#    define HURCHALLA_TARGET_ISA_ARM_32 1
#  endif
#endif


#ifndef HURCHALLA_TARGET_BIT_WIDTH
#  if defined(__x86_64__) || defined(_M_X64)
#    define HURCHALLA_TARGET_BIT_WIDTH 64
#  elif defined(__i386) || defined(_M_IX86)
#    define HURCHALLA_TARGET_BIT_WIDTH 32
#  elif defined(__aarch64__) || defined(_M_ARM64)
#    define HURCHALLA_TARGET_BIT_WIDTH 64
#  elif defined(__arm__) || defined(_M_ARM)
#    define HURCHALLA_TARGET_BIT_WIDTH 32
#  elif defined(__powerpc__) || defined(__ppc__) || defined(__PPC__)
#    if defined(__powerpc64__) || defined(__ppc64__) || defined(__PPC64__) || \
               defined(_ARCH_PPC64) || defined(__64BIT__) || defined(_LP64) || \
               defined(__LP64__)
#      define HURCHALLA_TARGET_BIT_WIDTH 64
#    else
#      define HURCHALLA_TARGET_BIT_WIDTH 32
#    endif
#  elif defined(__ia64) || defined(__itanium__) || defined(_M_IA64) || \
                                                               defined(__ia64__)
#    define HURCHALLA_TARGET_BIT_WIDTH 64
#  else
     // fallback if we couldn't find the target ALU's native bit depth after
     // looking at predefined compiler macros
#    include <cstdint>
#    if SIZE_MAX == UINT64_MAX
#      define HURCHALLA_TARGET_BIT_WIDTH 64
#    elif SIZE_MAX == UINT32_MAX
#      define HURCHALLA_TARGET_BIT_WIDTH 32
#    elif SIZE_MAX == UINT16_MAX
#      define HURCHALLA_TARGET_BIT_WIDTH 16
#    else
#      error "HURCHALLA_TARGET_BIT_WIDTH was undefined, and couldn't be deduced"
#    endif
#  endif

   // sanity check the setting we deduced for HURCHALLA_TARGET_BIT_WIDTH
#  include <cstddef>
#  include <climits>
   static_assert(sizeof(std::size_t) * CHAR_BIT == HURCHALLA_TARGET_BIT_WIDTH,
      "[This may be a false positive, but] This error suggests that the "
      "preprocessor logic in this header file incorrectly set "
      "HURCHALLA_TARGET_BIT_WIDTH.  You can predefine "
      "HURCHALLA_TARGET_BIT_WIDTH yourself, or change this file");
#else
   // sanity check the predefined value set for HURCHALLA_TARGET_BIT_WIDTH
#  include <type_traits>
   static_assert(std::is_integral<decltype(HURCHALLA_TARGET_BIT_WIDTH)>::value
                 && HURCHALLA_TARGET_BIT_WIDTH > 0,
      "When you predefine the macro HURCHALLA_TARGET_BIT_WIDTH, it must be a "
      "positive integer");
#endif


// In theory, the macro __SIZEOF_INT128__ indicates if __int128 is supported,
// but some compilers supported __int128 prior to providing that macro.  Clang
// had it at least since 3.0 and icc got it in 13.0.  Gcc has had since at least
// 4.1.  See
// https://stackoverflow.com/questions/16088282/is-there-a-128-bit-integer-in-gcc
// Note also that clang and icc define the GNUC and GNUC_MINOR macros.  Clang
// v3.0 defines them as 4,2.  Icc v13 defines them as 4,7.
//
// The macro  HURCHALLA_COMPILER_HAS_UINT128_T()  lets us know if __uint128_t is
// supported.
#if (HURCHALLA_TARGET_BIT_WIDTH < 64)
#  define HURCHALLA_COMPILER_HAS_UINT128_T() 0
#elif defined(__SIZEOF_INT128__) || \
             (defined(__clang_major__) && (__clang_major__ >= 3)) ||  \
             (defined(__INTEL_COMPILER) && (__INTEL_COMPILER >= 1300))
#  define HURCHALLA_COMPILER_HAS_UINT128_T() 1
#elif defined(__GNUC__) && defined(__GNUC_MINOR__) && \
             !defined(__INTEL_COMPILER) && !defined(__clang__) && \
             (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1))
#  define HURCHALLA_COMPILER_HAS_UINT128_T() 1
#else
#  define HURCHALLA_COMPILER_HAS_UINT128_T() 0
#endif


#if defined(__has_builtin)
#  define HURCHALLA_COMPILER_HAS_BUILTIN(builtin) __has_builtin(builtin)
#else
#  define HURCHALLA_COMPILER_HAS_BUILTIN(builtin) 0
#endif
#if HURCHALLA_COMPILER_HAS_BUILTIN(__builtin_expect) || \
    ( defined(__GNUC__) && !defined(__INTEL_COMPILER) && \
            !defined(__clang__) && __GNUC__ >= 3 ) || \
    ( defined(__INTEL_COMPILER) && __INTEL_COMPILER >= 1300 )
#  define HURCHALLA_LIKELY_IF(cond)   if (__builtin_expect(!!(cond), 1))
#  define HURCHALLA_UNLIKELY_IF(cond) if (__builtin_expect(!!(cond), 0))
#elif defined(__has_cpp_attribute)
#  if __has_cpp_attribute(likely)
#    define HURCHALLA_LIKELY_IF(cond)   if (!!(cond)) [[likely]]
#  else
#    define HURCHALLA_LIKELY_IF(cond)   if (!!(cond))
#  endif
#  if __has_cpp_attribute(unlikely)
#    define HURCHALLA_UNLIKELY_IF(cond) if (!!(cond)) [[unlikely]]
#  else
#    define HURCHALLA_UNLIKELY_IF(cond) if (!!(cond))
#  endif
#else
#  define HURCHALLA_LIKELY_IF(cond)   if (!!(cond))
#  define HURCHALLA_UNLIKELY_IF(cond) if (!!(cond))
#endif


#ifdef HURCHALLA_ALLOW_INLINE_ASM_ALL
#  undef HURCHALLA_DISALLOW_INLINE_ASM_MODMUL
#  undef HURCHALLA_ALLOW_INLINE_ASM_MODMUL
#  define HURCHALLA_ALLOW_INLINE_ASM_MODMUL 1
#else
   // By default, enable the inline asm modmult unless explicitly disallowed.
   // The x86 asm version is many times faster than the non-asm version, and for
   // asm it is relatively simple (though all inline asm is extremely difficult
   // to verify).
#  undef HURCHALLA_ALLOW_INLINE_ASM_MODMUL
#  ifndef HURCHALLA_DISALLOW_INLINE_ASM_MODMUL
#    define HURCHALLA_ALLOW_INLINE_ASM_MODMUL 1
#  endif
#endif


#endif
