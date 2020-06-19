
macro(EnableMaxWarnings target)


# Set compiler warnings.  We want to be able to compile cleanly with as many
# warnings enabled as possible, so that we know users will not get compiler
# warnings from our headers.

if(MSVC)
    target_compile_options(${target} PRIVATE /W4 /WX /Wall)
else()
    target_compile_options(${target} PRIVATE -Werror -Wall -Wextra
          -pedantic-errors -Wshadow -Wcast-qual -Wmissing-include-dirs
          -Wnon-virtual-dtor -Wconversion -Wsign-conversion -Wfloat-equal
          -Winvalid-pch -Wuninitialized -Wunused -Wunused-parameter -Wformat=2
          -Wdisabled-optimization -Wvla -Wundef -Wzero-as-null-pointer-constant
          -Wmissing-declarations)

    # clang or gcc
    if ((CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR
               (CMAKE_CXX_COMPILER_ID STREQUAL "GNU"))
        target_compile_options(${target} PRIVATE
                -Weffc++ -Wpedantic -Wredundant-decls -Wstrict-overflow=4
                -Wstack-protector -Wpacked -Wnull-dereference -Wregister
                -Wold-style-cast -Wdouble-promotion -Wformat-nonliteral
                -Wmissing-noreturn -Wctor-dtor-privacy -Winline)
    endif()

    # clang, or gcc10(or higher)
    if((CMAKE_CXX_COMPILER_ID MATCHES "Clang") OR
               ((CMAKE_CXX_COMPILER_ID STREQUAL "GNU") AND
               (CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.1)))
        target_compile_options(${target} PRIVATE
                -Wextra-semi)
    endif()


    # additional compiler specific warnings

    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
      # CLANG -----------
        target_compile_options(${target} PRIVATE  -ferror-limit=3
                -Wcast-align -Wmismatched-tags -Wabstract-vbase-init
                -Warray-bounds-pointer-arithmetic -Wassign-enum
                -Watomic-properties -Wauto-import -Wbinary-literal
                -Wc++14-compat-pedantic -Wc++14-extensions
                -Wc++17-compat-pedantic -Wc++17-extensions -Wclass-varargs
                -Wcomma -Wconditional-uninitialized -Wconsumed
                -Wcovered-switch-default -Wcuda-compat -Wdeprecated
                -Wduplicate-enum -Wformat-non-iso -Wformat-pedantic -Wgcc-compat
                -Wgnu -Wheader-hygiene -Widiomatic-parentheses
                -Wimplicitly-unsigned-literal
                -Winconsistent-missing-destructor-override -Wloop-analysis
                -Wmicrosoft -Wmissing-variable-declarations -Wnewline-eof
                -Wnon-gcc -Wnonportable-system-include-path
                -Wnullable-to-nonnull-conversion -Wopenmp-clauses -Wover-aligned
                -Wpadded -Wpedantic-core-features -Wpointer-arith
                -Wpointer-to-int-cast -Wpragma-pack -Wpragmas
                -Wprofile-instr-missing -Wredundant-parens -Wreserved-id-macro
                -Wreserved-user-defined-literal -Wretained-language-linkage
                -Wshadow-all -Wshift-sign-overflow -Wsigned-enum-bitfield
                -Wsource-uses-openmp -Wspir-compat -Wstatic-in-inline
                -Wstrict-selector-match -Wswitch-enum
                -Wtautological-constant-in-range-compare -Wthread-safety
                -Wthread-safety-negative -Wundefined-func-template
                -Wundefined-reinterpret-cast -Wunguarded-availability
                -Wunnamed-type-template-args -Wunreachable-code-aggressive
                -Wunsupported-dll-base-class-template
                -Wunused-exception-parameter -Wunused-member-function
                -Wunused-template -Wused-but-marked-unused -Wvector-conversion
                -Wwritable-strings)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.0)
            target_compile_options(${target} PRIVATE
                    -Walloca -Watomic-implicit-seq-cst -Wc++20-compat-pedantic
                    -Wc++20-extensions -Wctad-maybe-unsupported
                    -Wextra-semi-stmt -Wformat-type-confusion
                    -Wimplicit-int-float-conversion -Wmisexpect
                    -Wpoison-system-directories -Wnon-modular-include-in-module
                    -Wquoted-include-in-framework-header -Wsuspicious-memaccess)
        endif()

    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
      # GCC ----------------
        target_compile_options(${target} PRIVATE  -fmax-errors=3
                -Wlogical-op -Wnoexcept -Wduplicated-cond -Wduplicated-branches
                -Wstringop-overflow=3 -Wsuggest-final-methods -Wplacement-new=2
                -Wsuggest-final-types -Wsuggest-override -Walloc-zero
                -Wnormalized -Wformat-truncation=2 -Wformat-overflow=2
                -Wformat-signedness -Wstrict-null-sentinel -Wopenmp-simd
                -Wunsafe-loop-optimizations -Wvector-operation-performance)
        if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10.1)
            target_compile_options(${target} PRIVATE
                               -Wcast-align=strict -Wcomma-subscript -Wvolatile)
        else()
            target_compile_options(${target} PRIVATE
                               -Wcast-align)
        endif()

    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
      # INTEL C++ ----------------
        target_compile_options(${target} PRIVATE  -fmax-errors=3
                -pedantic -Wno-type-limits)
    endif()
endif()


endmacro()
