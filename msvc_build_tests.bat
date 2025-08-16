
@echo off
REM This Source Code Form is subject to the terms of the Mozilla Public
REM License, v. 2.0. If a copy of the MPL was not distributed with this
REM file, You can obtain one at https://mozilla.org/MPL/2.0/.

set build_dir=build\msvc

REM Example of how to use an earlier version of MSVC than the default:
REM cmake --help   (will show the available Generators you can use)
REM cmake -S. -B.\%build_dir% -DTEST_HURCHALLA_MODULAR_ARITHMETIC=ON -G "Visual Studio 15"
REM the above line appears to build x86-32.  To get x64:
REM cmake -S. -B.\%build_dir% -DTEST_HURCHALLA_MODULAR_ARITHMETIC=ON -G "Visual Studio 15 2017 Win64"

REM for Visual Studio 2019 and above, set the architecture with -A, for example:
REM -G "Visual Studio 16 2019" -A Win32
REM -G "Visual Studio 16 2019" -A x64
REM -G "Visual Studio 16 2019" -A ARM
REM -G "Visual Studio 16 2019" -A ARM64

cmake -S. -B.\%build_dir% -DTEST_HURCHALLA_MODULAR_ARITHMETIC=ON -DHURCHALLA_TEST_MODULAR_ARITHMETIC_HEAVYWEIGHT=ON -G "Visual Studio 17 2022" -A x64
if %errorlevel% neq 0 exit /b %errorlevel%
cmake --build .\%build_dir% --config Release
if %errorlevel% neq 0 exit /b %errorlevel%
cmake --build .\%build_dir% --config Debug
if %errorlevel% neq 0 exit /b %errorlevel%


%build_dir%\Release\test_hurchalla_modular_arithmetic.exe
if %errorlevel% neq 0 exit /b %errorlevel%

%build_dir%\Debug\test_hurchalla_modular_arithmetic.exe
if %errorlevel% neq 0 exit /b %errorlevel%
