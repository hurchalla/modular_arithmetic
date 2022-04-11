# The "Clockwork" Modular Arithmetic Library

![Alt text](images/clockxtrasmall_border2.jpg?raw=true "Clock Gears, photo by Krzysztof Golik, licensed CC BY-SA 4.0")

Clockwork is a high performance, easy to use Modular Arithmetic (header-only) library for C++ for up to 128 bit integer types, with extensive support for Montgomery arithmetic.  If you want or need Montgomery arithmetic, or general modular arithmetic functions, Clockwork is almost certainly the fastest and easiest library you can use.

## Design goals

The goal for Clockwork was to create a correct and flexible library with the best possible performance for modular arithmetic of native integer types.  For integer types that are double the native bit width (e.g. 128 bit), performance is still good but not quite as well optimized.  Larger than 128 bit types are permissible; however a library like GMP is likely to a better choice for such sizes.

## Status

Released.  All planned functionality and unit tests are finished and working correctly.

## Author

* **Jeffrey Hurchalla**

## License

This project is licensed under the MPL 2.0 License - see the [LICENSE.TXT](LICENSE.TXT) file for details

<br/>

## How to use the library

### With CMake

If you're using CMake for your project and you wish to add the Clockwork modular arithmetic library to it, then clone this git repository onto your system.  In your project's CMakeLists.txt file, add the following two lines with appropriate changes to their italic portions to match your project and paths ( an easy replacement for *your_binary_dir* is ${CMAKE_CURRENT_BINARY_DIR} ):  
add_subdirectory(*path_of_the_cloned_modular_arithmetic_repository* &nbsp; *your_binary_dir*/modular_arithmetic)  
target_link_libraries(*your_project_target_name* &nbsp; hurchalla_modular_arithmetic)  

For best performance you *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.  You can do this by calling CMake with -DCMAKE_BUILD_TYPE=Release.  

It may help to see a simple [example project with CMake](examples/example_with_cmake).

### Without CMake

If you're not using CMake for your project, you'll need to install/copy Clockwork's modular arithmetic headers and dependencies to some directory in order to use them.  To do this, first clone this git repository onto your system.  You'll need CMake on your system (at least temporarily), so install CMake if you don't have it.  Then from your shell run the following commands:  

>cd *path_of_the_cloned_modular_arithmetic_repository*  
>mkdir tmp  
>cd tmp  
>cmake -S.. -B.  
>cmake --install . --prefix *the_folder_you_want_to_install_to*  
If you prefer, for the last command you could instead use CMake's default install location (on linux this is /usr/local) by omitting the --prefix and subsequent folder.  

This will copy all the header files needed for this modular arithmetic library to an "include" subfolder in the installation folder of your choosing.
When compiling your project, you'll of course need to ensure that you have that include subfolder as part of your include path.  

For good performance you *must* ensure that the standard macro NDEBUG (see &lt;cassert&gt;) is defined when compiling.  You can generally do this by adding the option flag -DNDEBUG to your compile command.  

It may help to see a simple [example](examples/example_without_cmake).

## The API

Clockwork modular arithmetic is a header-only library, and the API is exposed by very short and simple header files (all headers not under any *detail* folder).  There are two main folder groupings: montgomery_arithmetic, and modular_arithmetic (i.e. standard non-montgomery).  A quick summary of the header files and functions is provided below; in all cases T is a template parameter of integral type.  Please view the header files for their documentation.  Probably the single most useful file is MontgomeryForm.h, discussed below.

From the modular_arithmetic group, the files *absolute_value_difference.h*, *modular_addition.h*, *modular_subtraction.h*, *modular_multiplication.h*, *modular_multiplicative_inverse.h*, and *modular_pow.h* provide the following functions, using standard (non-Montgomery) modular arithmetic:

*hurchalla::absolute_value_difference(T a, T b)*.  Returns the absolute value of (a-b), performed as if a and b are infinite precision signed ints.  
*hurchalla::modular_subtraction_prereduced_inputs(T a, T b, T modulus)*.  Let a conceptual "%%" operator represent a modulo operator that always returns a non-negative remainder. This function returns (a-b) %% modulus, performed as if a and b are infinite precision signed ints.  
*hurchalla::modular_addition_prereduced_inputs(T a, T b, T modulus)*.  Returns (a+b)%modulus, performed as if a and b have infinite precision and thus as if (a+b) is never subject to integer overflow.  
*hurchalla::modular_multiplication_prereduced_inputs(T a, T b, T modulus)*.   Returns (a\*b)%modulus, performed as if a and b have infinite precision.
*hurchalla::modular_multiplicative_inverse(T a, T modulus)*.  Returns the multiplicative inverse of a if it exists, and otherwise returns 0.  
*hurchalla::modular_pow(T base, T exponent, T modulus)*.  Returns the modular exponentiation of base^exponent (mod modulus).  

From the montgomery_arithmetic group, the file *MontgomeryForm.h* provides the easy to use (and zero cost abstraction) class *hurchalla::MontgomeryForm*, which has member functions for effortlessly performing operations in the Montgomery domain.  These operations include converting to/from Montgomery domain, add, sub, multiply, square, fused-multiply-add/sub, pow, gcd, and more.  For improved performance in some situations, the file montgomery_form_aliases.h provides simple aliases for faster (with limitations on allowed modulus) instantiations of the class MontgomeryForm.

For a simple demonstration of using MontgomeryForm, you can see one of the [examples](examples/example_without_cmake).

If you prefer not to use the high level interface of MontgomeryForm, and instead wish to directly call low level Montgomery arithmetic functions (such as REDC), the API header files within montgomery_arithmetic/low_level_api support all essential low level functions.

## Performance Notes

If you're interested in experimenting, predefining certain macros when compiling might improve performance - see [macros_for_performance.md](macros_for_performance.md).

## TODO

For the unit tests, solve the long compile time and high memory use (during compile) of the files test_MontgomeryForm.cpp and test_montgomery_pow.cpp.
