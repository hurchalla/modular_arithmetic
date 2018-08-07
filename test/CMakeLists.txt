
if(TARGET test_modular_arithmetic)
    return()
endif()

cmake_minimum_required(VERSION 3.12)
#project(test_modular_arithmetic VERSION 1.0.0 LANGUAGES CXX)



# The following ideas come from
#   https://github.com/google/googletest/blob/master/googletest/README.md
# They are further explained at https://crascit.com/2015/07/25/cmake-gtest/ and
# https://github.com/Crascit/DownloadProject/blob/master/DownloadProject.cmake
if (NOT TARGET gtest_main)
   # Download and unpack googletest at configure time
   configure_file(CMakeLists.txt.in
                  "${CMAKE_BINARY_DIR}/googletest/download/CMakeLists.txt")
   
   execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
     RESULT_VARIABLE result
     WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest/download")
   if(result)
     message(FATAL_ERROR "CMake step for googletest failed: ${result}")
   endif()
   
   execute_process(COMMAND ${CMAKE_COMMAND} --build .
     RESULT_VARIABLE result
     WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest/download")
   if(result)
     message(FATAL_ERROR "Build step for googletest failed: ${result}")
   endif()
   
   # Prevent overriding the parent project's compiler/linker settings on Windows
   set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
   
   # Add googletest to our build. This defines gtest and gtest_main targets.
   add_subdirectory("${CMAKE_BINARY_DIR}/googletest/src"
                    "${CMAKE_BINARY_DIR}/googletest/build"
                    EXCLUDE_FROM_ALL)
endif()


#set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/tests)
#set(CTEST_BINARY_DIRECTORY ${PROJECT_BINARY_DIR}/tests)


# needed for gtest_discover_tests()
include(GoogleTest)


add_executable(test_modular_arithmetic  test_ma.cpp)
target_link_libraries(test_modular_arithmetic
                      modular_arithmetic
                      gtest_main)
#add_test(test_modular_arithmetic  test_modular_arithmetic)
gtest_discover_tests(test_modular_arithmetic)

add_executable(test_modular_arithmetic2  test_ma2.cpp)
target_link_libraries(test_modular_arithmetic2
                      modular_arithmetic
                      gtest_main)
#add_test(test_modular_arithmetic2  test_modular_arithmetic2)
gtest_discover_tests(test_modular_arithmetic2)


#make test
#ctest --output-on-failure
#--schedule-random


#If a Makefile has a target named check, then make check will "build" that target. It's typically a phony target, meaning that it is a make-scripted command rather than a file named "check" that gets created.
#The gnu project advises that all gnu software should include a make check target that runs post-build tests in the build directory.


#undocumented CMake options:
# -H sets the source directory for the project.
# The -B option (also undocumented) sets the binary directory.

