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


set(public_header_files
    include/hurchalla/modular_arithmetic/ma.h
    )

add_library(modular_arithmetic STATIC src/ma.cpp ${public_header_files})

target_include_directories(modular_arithmetic
            PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_include_directories(modular_arithmetic SYSTEM
            PUBLIC $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)
#use this instead?
#target_include_directories(modular_arithmetic
#            PUBLIC $<INSTALL_INTERFACE:include>)


if(NOT DEFINED default_value_pbc_wrap_stdlib_assert)
    set(default_value_pbc_wrap_stdlib_assert ON)
endif()
add_subdirectory(../programming_by_contract
            ${CMAKE_CURRENT_BINARY_DIR}/programming_by_contract)

#add_subdirectory(../common ${CMAKE_CURRENT_BINARY_DIR}/common)


target_link_libraries(modular_arithmetic PRIVATE programming_by_contract)
#target_link_libraries(modular_arithmetic PUBLIC common)
# if you have a public dependency - something your implementation depends upon,
# and that is also visible through your public headers, use PUBLIC above instead
# of PRIVATE




# ***Tests***

# call enable_testing() just once, in the top level CMakeLists.txt
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(BUILD_TESTING_ONLY_CURRENT_PROJECT "Build only tests for the current (top level) project." OFF)
    option(BUILD_TESTING "Build the testing tree." OFF)
    #include(CTest)
    enable_testing()
    if(${BUILD_TESTING} OR ${BUILD_TESTING_ONLY_CURRENT_PROJECT})
        add_subdirectory(test/)
    endif()
#elseif(${BUILD_TESTING} AND (NOT ${BUILD_TESTING_ONLY_CURRENT_PROJECT}) AND ${CMAKE_TESTING_ENABLED})
elseif(${BUILD_TESTING} AND (NOT ${BUILD_TESTING_ONLY_CURRENT_PROJECT}))
    add_subdirectory(test/)
endif()

