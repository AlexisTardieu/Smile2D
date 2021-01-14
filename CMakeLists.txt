cmake_minimum_required(VERSION 2.6)

project(Smile2D)
include(FindPkgConfig)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Werror -O3 -g -std=c++11 -Wno-cast-function-type -Wno-unused -Wno-missing-field-initializers -Wno-deprecated-copy")

# message(STATUS "COMPILE OPTIONS : ${CMAKE_CXX_FLAGS}")

set(RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/build)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/build)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/libs)

set(CMAKE_BUILD_TYPE Release)

set(LIB_DIR ${CMAKE_SOURCE_DIR}/libs)

include(${LIB_DIR}/lib.cmake)

# Add Smile2D
add_subdirectory (${CMAKE_SOURCE_DIR}/src)
