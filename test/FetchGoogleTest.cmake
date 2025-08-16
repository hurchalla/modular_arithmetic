# Copyright (c) 2020-2022 Jeffrey Hurchalla.
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at https://mozilla.org/MPL/2.0/.


# Use FetchContent to get Googletest
# https://cmake.org/cmake/help/git-master/module/FetchContent.html

# Also inspired by https://github.com/Crascit/DownloadProject/issues/29
# and https://github.com/CLIUtils/cmake/blob/master/AddGoogleTest.cmake


if (NOT TARGET gtest_main)
    # Prevent GoogleTest from overriding our compiler/linker options
    # when building with Visual Studio
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF)


    include(FetchContent)

    if(FORCE_TEST_HURCHALLA_CPP11_STANDARD)
        # googletest v1.12 is the final release that supports C++11
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG        release-1.12.1
        )
        set(CMAKE_POLICY_VERSION_MINIMUM 3.10)
        #
        # when using GIT_TAG release-1.12.1, setting CMAKE_POLICY_VERSION_MINIMUM
        # to 3.10 avoids deprecation warnings (gtest v1.12's cmakelists.txt
        # evidently set cmake_minimum_required to somethimg under 3.10).
    else()
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG        main
        )
    endif()


    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

#    OLD METHOD (deprecated):
#
#    FetchContent_GetProperties(googletest)
#    if(NOT googletest_POPULATED)
#        FetchContent_Populate(googletest)
#    endif()


#    set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL ON)
#
#    # Needed because googletest project sets minimum CMake version < 3.0
#    set(CMAKE_POLICY_DEFAULT_CMP0048 OLD)
#
#    # Ignore undef warnings, issue with source file:
#    # googletest/include/gtest/internal/gtest-port.h:309:5: error: "_MSC_VER" is not defined
#    if(NOT MSVC)
#        add_compile_options( -Wno-undef )
#    endif()
#    set(CMAKE_INSTALL_DEFAULT_COMPONENT_NAME googletest)

#    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR})



#    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)



# alternative
#    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS 1 CACHE BOOL "")
#    add_subdirectory(${googletest_SOURCE_DIR} ${googletest_BINARY_DIR} EXCLUDE_FROM_ALL)
#    unset(CMAKE_SUPPRESS_DEVELOPER_WARNINGS)



#    add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} 
#            --force-new-ctest-process --output-on-failure)
#    set_target_properties(check PROPERTIES FOLDER "Scripts")

    mark_as_advanced(
    gmock_build_tests
    gtest_build_samples
    gtest_build_tests
    gtest_disable_pthreads
    gtest_force_shared_crt
    gtest_hide_internal_symbols
    BUILD_GMOCK
    BUILD_GTEST
    )

    set_target_properties(gtest gtest_main gmock gmock_main
        PROPERTIES FOLDER "Googletest")

#    if(MSVC AND MSVC_VERSION GREATER_EQUAL 1900)
#        target_compile_definitions(gtest PUBLIC _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
#        target_compile_definitions(gtest_main PUBLIC _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
#        target_compile_definitions(gmock PUBLIC _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
#        target_compile_definitions(gmock_main PUBLIC _SILENCE_TR1_NAMESPACE_DEPRECATION_WARNING)
#    endif()

endif()
