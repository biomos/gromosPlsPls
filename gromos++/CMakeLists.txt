cmake_minimum_required(VERSION 3.14)

# Read version from VERSION file
file(READ "${CMAKE_CURRENT_SOURCE_DIR}/VERSION" _version_file)
string(REGEX MATCH "GROMOS_MAJOR_VERSION=([0-9]+)" _ "${_version_file}")
set(GROMOS_MAJOR_VERSION ${CMAKE_MATCH_1})
string(REGEX MATCH "GROMOS_MINOR_VERSION=([0-9]+)" _ "${_version_file}")
set(GROMOS_MINOR_VERSION ${CMAKE_MATCH_1})
string(REGEX MATCH "GROMOS_MICRO_VERSION=([0-9]+)" _ "${_version_file}")
set(GROMOS_MICRO_VERSION ${CMAKE_MATCH_1})
set(GROMOS_VERSION "${GROMOS_MAJOR_VERSION}.${GROMOS_MINOR_VERSION}.${GROMOS_MICRO_VERSION}")

project(gromosxx VERSION ${GROMOS_VERSION} LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Required for shared library builds
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo" "MinSizeRel")
endif()

option(ENABLE_OPENMP "Enable OpenMP support" OFF)
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# Compiler-specific flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(COMPILER_GCC 1)
  set(CMAKE_CXX_FLAGS_RELEASE "-g -O3 -funroll-loops -fno-gcse -DNDEBUG -fomit-frame-pointer"
      CACHE STRING "GNU release flags" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall" CACHE STRING "GNU debug flags" FORCE)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
  set(COMPILER_ICC 1)
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ip -restrict -align -DNDEBUG"
      CACHE STRING "Intel release flags" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG "-g" CACHE STRING "Intel debug flags" FORCE)
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS CXX)
  add_compile_definitions(OMP)
endif()

find_package(Threads REQUIRED)
find_library(MATH_LIBRARY m)
find_package(ZLIB REQUIRED)
find_package(GSL REQUIRED)
set(HAVE_GSL 1)

find_package(Clipper)
if(CLIPPER_FOUND)
  set(HAVE_CLIPPER 1)
  set(FFTW_PREFIX "fftw3_")
else()
  set(FFTW_PREFIX "fftw_")
endif()

find_package(FFTW3)

find_package(MDPP)
if(MDPP_FOUND)
  set(HAVE_MDPP 1)
  set(GROMOSXX 1)
endif()

find_package(GMX)
if(GMX_FOUND)
  set(HAVE_GMX 1)
endif()

find_path(VMD_INCLUDE_DIR molfile_plugin.h)
if(VMD_INCLUDE_DIR)
  set(HAVE_VMD 1)
endif()

find_library(READLINE_LIBRARY readline)
find_path(READLINE_INCLUDE_DIR readline/readline.h)
if(READLINE_LIBRARY AND READLINE_INCLUDE_DIR)
  set(HAVE_READLINE 1)
endif()

include(CheckIncludeFileCXX)
check_include_file_cxx("wordexp.h" HAVE_WORDEXP_H)
if(HAVE_WORDEXP_H)
  include(CheckSymbolExists)
  check_symbol_exists(wordexp "wordexp.h" HAVE_WORDEXP)
  check_symbol_exists(wordfree "wordexp.h" HAVE_WORDFREE)
endif()

include(CheckSymbolExists)
set(CMAKE_REQUIRED_LIBRARIES "${MATH_LIBRARY}")
check_symbol_exists(hypot   "math.h" HAVE_HYPOT)
check_symbol_exists(expm1   "math.h" HAVE_EXPM1)
check_symbol_exists(acosh   "math.h" HAVE_ACOSH)
check_symbol_exists(asinh   "math.h" HAVE_ASINH)
check_symbol_exists(atanh   "math.h" HAVE_ATANH)
check_symbol_exists(isinf   "math.h" HAVE_ISINF)
check_symbol_exists(finite  "math.h" HAVE_FINITE)
check_symbol_exists(isfinite "math.h" HAVE_ISFINITE)
check_symbol_exists(isnan   "math.h" HAVE_ISNAN)
check_symbol_exists(log1p   "math.h" HAVE_LOG1P)
unset(CMAKE_REQUIRED_LIBRARIES)

string(TIMESTAMP GROMOS_DATE "%Y-%m-%dT%H:%M:%S")

# Generate config.h at the source root so that relative includes like
# #include "../../config.h" (used in src/pb/ and src/utils/) resolve correctly.
# This mirrors what automake does for in-source builds.
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/config.h.in"
  "${CMAKE_CURRENT_SOURCE_DIR}/config.h"
)

include(GNUInstallDirs)

# Make config.h findable via <config.h> as well (used in src/args/Arguments.cc)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}")

add_subdirectory(src)
add_subdirectory(programs)
add_subdirectory(contrib)
add_subdirectory(data)

enable_testing()

find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
  configure_file(doc/doxygen.conf.in doxygen.conf @ONLY)
  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxygen.conf
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM
  )
endif()
