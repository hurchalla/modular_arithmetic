# Modular Arithmetic
Modular Arithmetic library for C++

## Status

In development, currently alpha version.  All planned functionality is complete and works correctly with simple tests.  I need to implement full tests to be in beta.

## Author

* **Jeffrey Hurchalla**

## License

This project is licensed under the MIT License - see the [LICENSE.TXT](LICENSE.TXT) file for details

## TODO

Add Google tests.
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
