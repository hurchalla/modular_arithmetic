# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET modular_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.12)
project(modular_arithmetic VERSION 1.0.0 LANGUAGES C CXX)


# if this is the top level CMakeLists.txt, let IDEs group projects into folders.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
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


# Check for required C++ features only if the C++ compiler ID is non-empty.
#target_compile_features(modular_arithmetic INTERFACE 
#          $<$<NOT:$<CXX_COMPILER_ID:>>:
#               cxx_static_assert
#           >)




# By default, don't require clients of this library to know anything about
# the programming_by_contract library. If the default option value for
# HURCHALLA_INTERNAL_PBC_WRAP_STDLIB_ASSERT were OFF instead of ON below,
# clients of this library would be obligated by default to provide a custom
# assert handler.
option(HURCHALLA_INTERNAL_PBC_WRAP_STDLIB_ASSERT "Option HURCHALLA_INTERNAL_PBC_WRAP_STDLIB_ASSERT: Set to ON by default.  If you set HURCHALLA_INTERNAL_PBC_WRAP_STDLIB_ASSERT=OFF, then you will be able to (and required to) provide a custom assert handler, and you will have full access to the programming by contract functionality.  You can set this option via command line 'cmake -DHURCHALLA_INTERNAL_PBC_WRAP_STDLIB_ASSERT=OFF ...')" ON)
if(HURCHALLA_INTERNAL_PBC_WRAP_STDLIB_ASSERT)
    target_compile_definitions(modular_arithmetic PRIVATE
              PBC_WRAP_STDLIB_ASSERT
              )
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

# if this is the top level CMakeLists.txt
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(BUILD_TESTING_ONLY_CURRENT_PROJECT "Build only tests for the current (top level) project." OFF)
    option(BUILD_TESTING "Build the testing tree." OFF)
    if(${BUILD_TESTING_ONLY_CURRENT_PROJECT} AND (NOT ${BUILD_TESTING}))
        set(BUILD_TESTING ON CACHE BOOL "" FORCE)
    endif()
    if(${BUILD_TESTING} OR ${BUILD_TESTING_ONLY_CURRENT_PROJECT})
        enable_testing()
        #include(CTest)
        add_subdirectory(test)
    endif()
elseif(${BUILD_TESTING} AND (NOT ${BUILD_TESTING_ONLY_CURRENT_PROJECT}))
    # If a parent project (to this project) calls include(CTest), the included
    # CTest.cmake will define a BUILT_TESTING option and will call
    # enable_testing if BUILD_TESTING is ON.  My own (possibly parent) projects
    # may set BUILD_TESTING ON - if they do, they will call include(CTest) or
    # call enable_testing().
    # So, if we are inside this elseif clause, we know enable_testing() has
    # already been called.  We will sanity-check this now.
    if (NOT ${CMAKE_TESTING_ENABLED})
        message(FATAL_ERROR "Fatal error: enable_testing() was not called")
    endif()
    add_subdirectory(test)
endif()

