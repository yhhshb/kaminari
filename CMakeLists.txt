cmake_minimum_required(VERSION 3.11)
project(KAMINARI)

include(CMakePrintHelpers)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
endif ()

MESSAGE(STATUS "CMAKE_BUILD_TYPE: " ${CMAKE_BUILD_TYPE})

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if (UNIX AND (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL "x86_64"))
  MESSAGE(STATUS "Compiling with flags: -march=native -mbmi2 -msse4.2")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native") # If working on x86_64
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2 -msse4.2") # for hardware popcount and pdep
endif()

if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-parameter -Wno-unused-function -Wno-unknown-attributes")
  endif()
  if (KAMINARI_USE_SANITIZERS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")
  endif()
endif()

## fetch latest argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)

include_directories(${PROJECT_SOURCE_DIR})
# include_directories(${PROJECT_SOURCE_DIR}/bundled/biolib)

## GGCAT
add_custom_target(
   ggcat_cpp_api
   COMMAND make
   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bundled/ggcat/crates/capi/ggcat-cpp-api
)
set(GGCAT_CPP_BINDINGS ${CMAKE_SOURCE_DIR}/bundled/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cpp_bindings.a)
set(GGCAT_CXX_INTEROP ${CMAKE_SOURCE_DIR}/bundled/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_cxx_interop.a)
set(GGCAT_API ${CMAKE_SOURCE_DIR}/bundled/ggcat/crates/capi/ggcat-cpp-api/lib/libggcat_api.a)

## kaminari 

set(BIOLIB_SOURCES
  bundled/biolib/src/io.cpp
  bundled/biolib/src/elias_fano.cpp
  bundled/biolib/src/kmer_view.cpp
  bundled/biolib/src/logtools.cpp
)

set(LPHASH_SOURCES
  bundled/lphash/lib/src/constants.cpp
  bundled/lphash/lib/src/minimizer.cpp
  bundled/lphash/lib/src/mphf.cpp
  bundled/lphash/lib/src/partitioned_mphf.cpp
  bundled/lphash/lib/src/quartet_wtree.cpp
  bundled/lphash/lib/src
)

set(KAMINARI_SOURCES
  src/build.cpp
  src/GGCAT.cpp
  src/hybrid.cpp
)

add_library(biolib_static STATIC
  ${BIOLIB_SOURCES}
)

# Create a static lib for LPHash
add_library(lphash_static STATIC
  ${LPHASH_SOURCES}
)

# add_subdirectory(bundled/biolib)
add_executable(kaminari src/kaminari.cpp ${KAMINARI_SOURCES})
add_dependencies(kaminari ggcat_cpp_api lphash_static biolib_static)
target_link_libraries(kaminari
  biolib_static
  lphash_static
  argparse
  z
  ${GGCAT}
  ${GGCAT_API}
  ${GGCAT_CPP_BINDINGS}
  ${GGCAT_CXX_INTEROP}
  ${CMAKE_DL_LIBS}
)

#----------------------------- Tests ----------------------------