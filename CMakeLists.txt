# --- This file is distributed under the MIT Open Source License, as detailed
# in the file "LICENSE.TXT" in the root of this repository ---

if(TARGET hurchalla_modular_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.13)

project(hurchalla_modular_arithmetic VERSION 1.0.0 LANGUAGES CXX)

# if this is the top level CMakeLists.txt, let IDEs group projects into folders
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()



# TODO: this section seems slightly messy for detecting/setting up testing
# --------------------
option(TEST_HURCHALLA_LIBS
       "Build the tests for all Hurchalla library projects."
       OFF)

# if this is the top level CMakeLists.txt, add testing options, and enable
# testing when testing options have been set to ON.
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    option(TEST_HURCHALLA_MODULAR_ARITHMETIC
        "Build the tests for the Hurchalla modular arithmetic library project."
        OFF)
    if(TEST_HURCHALLA_MODULAR_ARITHMETIC OR TEST_HURCHALLA_LIBS)
        enable_testing()
        # include(CTest)
    endif()
elseif(TEST_HURCHALLA_LIBS)
    # If TEST_HURCHALLA_LIBS is set to ON, enable_testing() should have been
    # called either directly or indirectly by the top level project. (Note that
    # if a project calls include(CTest), the included CTest.cmake defines a
    # BUILT_TESTING option and calls enable_testing if BUILD_TESTING is ON.)
    if (NOT CMAKE_TESTING_ENABLED)
        message(FATAL_ERROR "Fatal error: TEST_HURCHALLA_LIBS was set, but enable_testing() was never called")
    endif()
endif()



# TODO: I'm not sure whether these lines provide any benefit:
# --------------------
if(NOT DEFINED CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
if(NOT DEFINED CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()
if(NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
endif()



add_library(hurchalla_modular_arithmetic INTERFACE)


add_subdirectory(modular_arithmetic
                 ${CMAKE_CURRENT_BINARY_DIR}/modular_arithmetic)
add_subdirectory(montgomery_arithmetic
                 ${CMAKE_CURRENT_BINARY_DIR}/montgomery_arithmetic)

target_link_libraries(hurchalla_modular_arithmetic
                      INTERFACE hurchalla_basic_modular_arithmetic)
target_link_libraries(hurchalla_modular_arithmetic
                      INTERFACE hurchalla_montgomery_arithmetic)



# TODO:  The following is a guess at how installing target include directories
# should be done.  It may be wrong.  It assumes that the build step from the
# subdirectories montgomery_arithmetic and modular_arithmetic (which have
# build phase target_include_directories commands) provides the information to
# cmake which cmake then uses in the install phase target_include_directories
# below.  I haven't tested this or anything related to cmake install yet...
# ---------------------
target_include_directories(hurchalla_modular_arithmetic
                           INTERFACE $<INSTALL_INTERFACE:include>)
# TODO: use this instead?
# ----------------
#target_include_directories(hurchalla_modular_arithmetic SYSTEM
#                      INTERFACE $<INSTALL_INTERFACE:$<INSTALL_PREFIX>/include>)



# ***Tests***

# if this is the top level CMakeLists.txt
if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    if(TEST_HURCHALLA_MODULAR_ARITHMETIC OR TEST_HURCHALLA_LIBS)
        add_subdirectory(test)
    endif()
elseif(TEST_HURCHALLA_LIBS)
    add_subdirectory(test)
endif()
