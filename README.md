# Modular Arithmetic
Modular Arithmetic library for C++

## Design goals

A correct and flexible library with best possible performance for modular arithmetic of native integer types.  For integer types that are double the native bit width (e.g. 128 bit), performance should still be reasonably good though not optimal.  Larger than 128 bit types are permissible; however a library like GMP is likely to provide much better performance at such sizes.

## Status

Released.  All planned functionality and unit tests are finished and working correctly.  I have not yet updated the CMakeLists.tst to add compile definitions options that enable inline asm (inline asm is available on x86-64 only).  Until then, if you use this library and want the best possible performance, you may wish to consider passing in a macro definition flag for HURCHALLA_ALLOW_INLINE_ASM_REDC.  It is not enabled by default because inline asm is extremely difficult to verify for correctness.  While I believe I'm very skilled at writing high quality inline asm, I advise you to be skeptical; unit tests of inline asm are far less helpful than you might think- https://gcc.gnu.org/wiki/DontUseInlineAsm

## Author

* **Jeffrey Hurchalla**

## License

This project is licensed under the MPL 2.0 License - see the [LICENSE.TXT](LICENSE.TXT) file for details

<br/>

## How to use the library

### With CMake

If you're using CMake for your project and you wish to add this modular arithmetic library to it, then clone this git repository onto your system.  In your project's CMakeLists.txt file, add the following two lines with appropriate changes to their italic portions to match your project and paths ( an easy replacement for *your_binary_dir* is ${CMAKE_CURRENT_BINARY_DIR} ):  
add_subdirectory(*path_of_the_cloned_modular_arithmetic_repository* &nbsp; *your_binary_dir*/modular_arithmetic)  
target_link_libraries(*your_project_target_name* &nbsp; hurchalla_modular_arithmetic)  

For best performance you *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.  You can do this by calling CMake with -DCMAKE_BUILD_TYPE=Release.  

It may help to see a simple [example project with CMake](example_with_cmake).

### Without CMake

If you're not using CMake for your project, you'll need to install/copy these modular arithmetic headers and dependencies to some directory in order to use them.  To do this, first clone this git repository onto your system.  You'll need CMake on your system (at least temporarily), so install CMake if you don't have it.  Then from your shell run the following commands:  

>cd *path_of_the_cloned_modular_arithmetic_repository*  
>mkdir tmp  
>cd tmp  
>cmake -S.. -B.  
>cmake --install . --prefix *the_folder_you_want_to_install_to*  
If you prefer, for the last command you could instead use CMake's default install location (on linux this is /usr/local) by omitting the --prefix and subsequent folder.  

This will copy all the header files needed for the modular arithmetic library to an "include" subfolder in the installation folder of your choosing.
When compiling your project, you'll of course need to ensure that you have that include subfolder as part of your include path.  

For good performance you *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.  You can generally do this by adding the option flag -DNDEBUG to your compile command.  

It may help to see a simple [example](example_without_cmake).

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

HURCHALLA_ALLOW_INLINE_ASM_ABSDIFF
HURCHALLA_ALLOW_INLINE_ASM_MODADD
HURCHALLA_ALLOW_INLINE_ASM_MODSUB
HURCHALLA_ALLOW_INLINE_ASM_QUARTERRANGE_GET_CANONICAL
HURCHALLA_ALLOW_INLINE_ASM_HALFRANGE_GET_CANONICAL
HURCHALLA_ALLOW_INLINE_ASM_REDC

for testing use only:
HURCHALLA_ALLOW_INLINE_ASM_ALL
HURCHALLA_AVOID_CSELECT

Document that the INLINE_ASM macros above may or may not improve performance.  You need to benchmark with different ASM macros defined/not defined, and generally you would want to start with simply comparing performance with HURCHALLA_ALLOW_INLINE_ASM_REDC defined or not defined.  On x86_64 intel, in brief testing defining HURCHALLA_ALLOW_INLINE_ASM_REDC, I found gcc saw ~5% improvement and clang was essentially unaffected.

For MSVC optimization: /Gy (function-level linking) and /Gw (global data optimization) compiler switches
https://docs.microsoft.com/en-us/archive/msdn-magazine/2015/february/compilers-what-every-programmer-should-know-about-compiler-optimizations
