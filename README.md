# Modular Arithmetic
Modular Arithmetic library for C++

## Design goals

A flexible and rock solid library for modular arithmetic of up to 128 bit types, that exceeds the performance of any known solutions.  Larger than 128 bit types are usable also but not optimized.

## Status

In beta.  All planned functionality and unit tests are finished and working correctly.  I have not yet updated the CMakeLists.tst to add compile definitions options that enable inline asm (inline asm is available on x86-64 only).  Until then, if you use this library and want the best possible performance, you may wish to consider passing in a macro definition flag for HURCHALLA_ALLOW_INLINE_ASM_ALL.  It is not enabled by default because inline asm is extremely difficult to verify for correctness.  While I believe I'm very skilled at writing high quality inline asm, I advise you to be skeptical; unit tests of inline asm are far less helpful than you might think- https://gcc.gnu.org/wiki/DontUseInlineAsm

## Author

* **Jeffrey Hurchalla**

## License

This project is licensed under the MIT License - see the [LICENSE.TXT](LICENSE.TXT) file for details

## TODO

Set up continuous integration (Travis or GitHub actions)
Document the project.
Compare performance of impl_modular_multiplication_prereduced_inputs(uint64_t, uint64_t, uint64_t) with internal __uint128_t, to the template function version.
Add the following options/target_compile_definitions to CMakeLists.txt, and document here:
HURCHALLA_COMPILE_ERROR_ON_SLOW_MATH
HURCHALLA_TARGET_ISA_HAS_NO_DIVIDE
HURCHALLA_TARGET_BIT_WIDTH
HURCHALLA_TARGET_ISA_X86_32
HURCHALLA_TARGET_ISA_X86_64
HURCHALLA_TARGET_ISA_ARM_32
HURCHALLA_TARGET_ISA_ARM_64

HURCHALLA_ALLOW_INLINE_ASM_MODADD
HURCHALLA_ALLOW_INLINE_ASM_MODSUB
HURCHALLA_ALLOW_INLINE_ASM_MONTMUL
HURCHALLA_DISALLOW_INLINE_ASM_MODMUL
HURCHALLA_ALLOW_INLINE_ASM_ALL

For MSVC optimization: /Gy (function-level linking) and /Gw (global data optimization) compiler switches
https://docs.microsoft.com/en-us/archive/msdn-magazine/2015/february/compilers-what-every-programmer-should-know-about-compiler-optimizations
