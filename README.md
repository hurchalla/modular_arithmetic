# Modular Arithmetic
Modular Arithmetic library for C++

## Design goals

A correct and flexible library for modular arithmetic of native integer types, achieving best possible performance.  For integer types that are double the native bit width (e.g. 128 bit), performance should still be reasonably good though not optimal.  Larger than 128 bit types are permissible; however a library like GMP is likely to provide much better performance at such sizes.

## Status

Released.  All planned functionality and unit tests are finished and working correctly.  I have not yet updated the CMakeLists.tst to add compile definitions options that enable inline asm (inline asm is available on x86-64 only).  Until then, if you use this library and want the best possible performance, you may wish to consider passing in a macro definition flag for HURCHALLA_ALLOW_INLINE_ASM_ALL.  It is not enabled by default because inline asm is extremely difficult to verify for correctness.  While I believe I'm very skilled at writing high quality inline asm, I advise you to be skeptical; unit tests of inline asm are far less helpful than you might think- https://gcc.gnu.org/wiki/DontUseInlineAsm

## Author

* **Jeffrey Hurchalla**

## License

This project is licensed under the MIT License - see the [LICENSE.TXT](LICENSE.TXT) file for details

## TODO

Determine whether users of this library can safely use uint128_t with gcc versions between 5.1 and 11.  See https://gcc.gnu.org/bugzilla/show_bug.cgi?id=98474
Set up continuous integration (Travis or GitHub actions).
Document the project.
Solve the long compile time and high memory use (during compile) for the files test_MontgomeryForm.cpp and test_montgomery_pow.cpp.
Compare performance of impl_modular_multiplication_prereduced_inputs(uint64_t, uint64_t, uint64_t) with internal __uint128_t, to the template function version.
Add the following options/target_compile_definitions to CMakeLists.txt, and document here:
HURCHALLA_COMPILE_ERROR_ON_SLOW_MATH
HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE
HURCHALLA_TARGET_BIT_WIDTH
HURCHALLA_TARGET_ISA_X86_32
HURCHALLA_TARGET_ISA_X86_64
HURCHALLA_TARGET_ISA_ARM_32
HURCHALLA_TARGET_ISA_ARM_64

HURCHALLA_DISALLOW_INLINE_ASM_MODMUL
HURCHALLA_ALLOW_INLINE_ASM_ALL

HURCHALLA_ALLOW_INLINE_ASM_ABSDIFF
HURCHALLA_ALLOW_INLINE_ASM_MODADD
HURCHALLA_ALLOW_INLINE_ASM_MODSUB
HURCHALLA_ALLOW_INLINE_ASM_MONTHELPER
HURCHALLA_ALLOW_INLINE_ASM_REDC
experimental:
HURCHALLA_ALLOW_INLINE_ASM_MONTADD_SQRT_RANGE
HURCHALLA_ALLOW_INLINE_ASM_MONTSUB_SQRT_RANGE

Document that the INLINE_ASM macros above may or may not improve performance.  You need to benchmark with different ASM macros defined/not defined, and generally you would want to start with simply comparing performance with HURCHALLA_ALLOW_INLINE_ASM_ALL defined or not defined.  On x86_64 intel, in brief testing defining HURCHALLA_ALLOW_INLINE_ASM_ALL, I found gcc got a 0-8% improvement, and clang suffered a 0-20% loss of performance.

For MSVC optimization: /Gy (function-level linking) and /Gw (global data optimization) compiler switches
https://docs.microsoft.com/en-us/archive/msdn-magazine/2015/february/compilers-what-every-programmer-should-know-about-compiler-optimizations
