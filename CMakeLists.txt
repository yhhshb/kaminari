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
    GIT_REPOSITORY https://github.com/yhhshb/argparse.git
)
FetchContent_MakeAvailable(argparse)

include_directories(${PROJECT_SOURCE_DIR})

## kaminari 

set(BIOLIB_SOURCES
  bundled/biolib/src/io.cpp
  bundled/biolib/src/elias_fano.cpp
  bundled/biolib/src/kmer_view.cpp
  bundled/biolib/src/logtools.cpp
  bundled/biolib/bundled/MurmurHash3.cpp
)

set(KAMINARI_SOURCES
  src/utils.cpp
  src/build.cpp
  src/query.cpp
  src/hybrid.cpp
)

add_library(biolib_static STATIC
  ${BIOLIB_SOURCES}
)

# add_subdirectory(bundled/biolib)
add_executable(kaminari src/kaminari.cpp ${KAMINARI_SOURCES})

# Include FQFeeder
add_subdirectory(bundled/FQFeeder)

# Include FQFeeder
add_subdirectory(bundled/unordered_dense)

# Add pthash as a subdirectory
add_subdirectory(bundled/pthash)

add_subdirectory(bundled/Minimizers)

# pthash includes
target_include_directories(kaminari PRIVATE
    ${PROJECT_SOURCE_DIR}/bundled/pthash
    ${PROJECT_SOURCE_DIR}/bundled/pthash/include
    ${PROJECT_SOURCE_DIR}/bundled/pthash/external/essentials/include
    ${PROJECT_SOURCE_DIR}/bundled/pthash/external/xxHash
)

add_dependencies(kaminari biolib_static) # lphash_static

target_link_libraries(kaminari
  biolib_static
  argparse
  fqfeeder
  PTHASH
  unordered_dense
  z
)

if (LINUX)
  target_link_libraries(kaminari
    rt
  )
endif()

#----------------------------- Tests ----------------------------


# Add test.cpp as a new executable called test
add_executable(test src/test.cpp)

# Link necessary libraries to the test executable
target_link_libraries(test
  biolib_static
  argparse
  fqfeeder
  PTHASH
  unordered_dense
  z
)

# If on Linux, link the rt library for the test executable as well
if (LINUX)
  target_link_libraries(test rt)
endif()