# Modular Arithmetic
Modular Arithmetic library for C++

## Design goals

A flexible and rock solid library for modular arithmetic of native types, achieving best known performance.  Performance should be reasonably good but not optimal, for double the native bit width, e.g. for 128 bit types.  Larger than 128 bit types are permissible; however a library like GMP is likely to provide much better performance at such sizes.

## Status

In beta.  All planned functionality and unit tests are finished and working correctly.  I have not yet updated the CMakeLists.tst to add compile definitions options that enable inline asm (inline asm is available on x86-64 only).  Until then, if you use this library and want the best possible performance, you may wish to consider passing in a macro definition flag for HURCHALLA_ALLOW_INLINE_ASM_ALL.  It is not enabled by default because inline asm is extremely difficult to verify for correctness.  While I believe I'm very skilled at writing high quality inline asm, I advise you to be skeptical; unit tests of inline asm are far less helpful than you might think- https://gcc.gnu.org/wiki/DontUseInlineAsm

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

For MSVC optimization: /Gy (function-level linking) and /Gw (global data optimization) compiler switches
https://docs.microsoft.com/en-us/archive/msdn-magazine/2015/february/compilers-what-every-programmer-should-know-about-compiler-optimizations

in some README (maybe this one) talk about Montgomery special multiply, takes T and V, returns T.  but probably don't implement it in MontgomeryForm, since it makes the MF API bigger (and thus less simple) without adding much value.  It can be done via low level API instead.
Also talk about montgomery multiply that takes a W and T (where the type W value is a montmul of a T with R_cubed_mod_n, or a montmul of a V with R_squared_mod_n), and returns V.  Again this is interesting but probably not all that useful in practice, and this first way to get a W would require an extra computation of R_cubed_mod_n in the MF constructor.  Regardless this kind of special montgomery multiply adds complexity to the MF API without adding much value I can see.  It can be done instead via the low level API.

Get rid of MontySixthRange and HalfRange and just use two versions of famul() within each of FullRange and QuarterRange, possibly by making famul templated or tagged to correspond to full/half/quarter/sixth.
