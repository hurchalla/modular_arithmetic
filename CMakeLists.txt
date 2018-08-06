cmake_minimum_required(VERSION 3.10)
project(modular_arithmetic VERSION 1.0.0 LANGUAGES C CXX)

if(TARGET modular_arithmetic)
    return()
endif()


if(NOT DEFINED CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
if(NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

#
set(public_header_files
    include/hurchalla/modular_arithmetic/ma.h
    )

add_library(modular_arithmetic STATIC src/ma.cpp ${public_header_files})

target_include_directories(modular_arithmetic PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_include_directories(modular_arithmetic SYSTEM PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)
#use this instead?
#target_include_directories(modular_arithmetic PUBLIC $<INSTALL_INTERFACE:include>)


if(NOT DEFINED default_value_pbc_wrap_stdlib_assert)
    set(default_value_pbc_wrap_stdlib_assert ON)
endif()
add_subdirectory(../programming_by_contract ${CMAKE_CURRENT_BINARY_DIR}/programming_by_contract)

#add_subdirectory(../common ${CMAKE_CURRENT_BINARY_DIR}/common)


target_link_libraries(modular_arithmetic PRIVATE programming_by_contract)
#target_link_libraries(modular_arithmetic PUBLIC common)
#  if you have a public dependency - something your implementation depends upon, and that is also
#  visible through your public headers, use PUBLIC above instead of PRIVATE



# ***Test***

# include(CTest) defines in the option BUILD_TESTING, which allows to control whether to build all tests of the project or not.
# CMake will generate tests only if the enable_testing() command has been invoked. The CTest module invokes the command automatically when the BUILD_TESTING option is ON.

option(PBC_BUILD_TEST "build test" OFF)
option(BUILD_TESTING "build test" OFF) # The ctest variable for building tests
if(${PBC_BUILD_TEST} OR ((CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR) AND ${BUILD_TESTING}))
    enable_testing()
#    include(CTest)
    add_subdirectory(test/)
endif()
