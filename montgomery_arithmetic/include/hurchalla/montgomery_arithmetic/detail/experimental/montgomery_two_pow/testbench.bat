@echo off
setlocal enabledelayedexpansion

REM ============================================
REM !!README!!
REM You MUST run this bat script from the Visual
REM Studio x64 Native Tools Command Prompt.  (I
REM used Windows Start Menu/Visual Studio 2022/
REM Visual Studio x64 Native Tools Command
REM Prompt for VS 2022).
REM ============================================

REM ============================================
REM CONFIGURATION
REM ============================================

set repo_directory=C:\Users\jeff\Documents\repos

set cppcompiler=cl
set cpp_standard=/std:c++17

REM Common include paths
set includes=/I"%repo_directory%\modular_arithmetic\modular_arithmetic\include" ^
            /I"%repo_directory%\modular_arithmetic\montgomery_arithmetic\include" ^
            /I"%repo_directory%\util\include"




REM ============================================
REM FUNCTION: run_test
REM %1 = output file
REM %2 = MONT type
REM %3 = TEST TYPE (SCALAR or ARRAY)
REM ============================================

goto :main

:run_test
set outfile=%1
set mont_type=%2
set test_type=%3

echo.
echo ============================================
echo Running: %mont_type%  %test_type%
echo ============================================
echo.

REM Build compile command
%cppcompiler% /O2 %cpp_standard% /EHsc ^
    /DDEF_MONT_TYPE=%mont_type% ^
    /DDEF_UINT_TYPE=uint64_t ^
    /DTEST_%test_type% ^
    %includes% ^
    -c testbench_montgomery_two_pow.cpp

if %ERRORLEVEL% NEQ 0 (
    echo Compilation failed!
    exit /b 1
)

REM Link
%cppcompiler% /O2 %cpp_standard% ^
    testbench_montgomery_two_pow.obj ^
    /Fe:testbench_montgomery_two_pow.exe

if %ERRORLEVEL% NEQ 0 (
    echo Link failed!
    exit /b 1
)

REM Run and capture output (NO POWERSHELL)
REM testbench_montgomery_two_pow.exe 191 8 22 > "%outfile%" 2>&1
REM  use the following to run the process at high priority:
start /b /high /wait testbench_montgomery_two_pow.exe 191 8 22 > "%outfile%" 2>&1


if %ERRORLEVEL% NEQ 0 (
    echo Execution failed!
    exit /b 1
)

echo Done: %outfile%

goto :eof


REM ============================================
REM MAIN TESTS
REM ============================================

:main

call :run_test 64_quarter_msvc_noasm_scalar.txt MontgomeryQuarter SCALAR
call :run_test 64_quarter_msvc_noasm_array.txt  MontgomeryQuarter ARRAY

call :run_test 64_half_msvc_noasm_scalar.txt    MontgomeryHalf SCALAR
call :run_test 64_half_msvc_noasm_array.txt     MontgomeryHalf ARRAY

call :run_test 64_full_msvc_noasm_scalar.txt    MontgomeryFull SCALAR
call :run_test 64_full_msvc_noasm_array.txt     MontgomeryFull ARRAY

echo.
echo All benchmarks complete.
exit /b 0
