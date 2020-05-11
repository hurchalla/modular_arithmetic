
#ifndef HURCHALLA_MONTGOMERY_ARITHMETIC_COMPILER_MACROS_H_INCLUDED
#define HURCHALLA_MONTGOMERY_ARITHMETIC_COMPILER_MACROS_H_INCLUDED


#ifdef _MSC_VER
#  define FORCE_INLINE __forceinline
#elif defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#  define FORCE_INLINE inline __attribute__((always_inline))
#else
#  define FORCE_INLINE
#endif

#if defined(__clang__) || defined(__INTEL_COMPILER)
#  define REQUEST_UNROLL_LOOP _Pragma("unroll")
#elif defined(__GNUC__) && __GNUC__ >= 8
#  define REQUEST_UNROLL_LOOP _Pragma("GCC unroll 127")
#else
#  define REQUEST_UNROLL_LOOP
#endif

#ifndef TARGET_BIT_WIDTH
#error "TARGET_BIT_WIDTH must be defined"
#endif
// In theory, the macro __SIZEOF_INT128__ indicates if __int128 is supported,
// but some compilers supported __int128 prior to providing that macro.  Clang
// had it at least since 3.0 and icc got it in 13.0.  Gcc has had since at least
// 4.1.  See
// https://stackoverflow.com/questions/16088282/is-there-a-128-bit-integer-in-gcc
// Note also that clang and icc define the GNUC and GNUC_MINOR macros.  Clang
// v3.0 defines them as 4,2.  Icc v13 defines them as 4,7.
//
// The macro  COMPILER_HAS_UINT128_T  lets us know if __uint128_t is supported.
#define COMPILER_HAS_UINT128_T (TARGET_BIT_WIDTH >= 64 && \
                ( defined(__SIZEOF_INT128__) || (__clang_major__ >= 3) || \
                  (__INTEL_COMPILER >= 1300) || \
                  ( !defined(__INTEL_COMPILER) && !defined(__clang_major__) && \
                    (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1)) )))


#endif
