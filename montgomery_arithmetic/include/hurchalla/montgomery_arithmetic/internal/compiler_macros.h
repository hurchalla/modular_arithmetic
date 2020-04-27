
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


#endif
