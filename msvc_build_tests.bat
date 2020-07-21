
@echo off
set build_dir=build\msvc

cmake -S. -B.\%build_dir% -DTEST_HURCHALLA_LIBS=ON
if %errorlevel% neq 0 exit /b %errorlevel%
cmake --build .\%build_dir% --config Release
if %errorlevel% neq 0 exit /b %errorlevel%
cmake --build .\%build_dir% --config Debug
if %errorlevel% neq 0 exit /b %errorlevel%

REM due to death tests, gtest makes the tests for programming_by_contract output
REM a ton of lines that start with "Running main()".  We'll filter those out.
%build_dir%\Release\test_programming_by_contract.exe > tmp_test_results.txt
set result=%errorlevel%
type tmp_test_results.txt | find /v "Running main()"
del tmp_test_results.txt
if %result% neq 0 exit /b %result%

%build_dir%\Release\test_ndebug_programming_by_contract.exe>tmp_test_results.txt
set result=%errorlevel%
type tmp_test_results.txt | find /v "Running main()"
del tmp_test_results.txt
if %result% neq 0 exit /b %result%

%build_dir%\Release\test_hurchalla_modular_arithmetic.exe
if %result% neq 0 exit /b %result%


REM Once again we'll filter out the useless lines starting with "Running main()"
%build_dir%\Debug\test_programming_by_contract.exe > tmp_test_results.txt
set result=%errorlevel%
type tmp_test_results.txt | find /v "Running main()"
del tmp_test_results.txt
if %result% neq 0 exit /b %result%

%build_dir%\Debug\test_ndebug_programming_by_contract.exe > tmp_test_results.txt
set result=%errorlevel%
type tmp_test_results.txt | find /v "Running main()"
del tmp_test_results.txt
if %result% neq 0 exit /b %result%

%build_dir%\Debug\test_hurchalla_modular_arithmetic.exe
if %result% neq 0 exit /b %result%
