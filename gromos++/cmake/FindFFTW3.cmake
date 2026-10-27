# FindFFTW3.cmake - Find the FFTW3 library
#
# Variables set:
#   FFTW3_FOUND        - True if found
#   FFTW3_INCLUDE_DIRS - Include directories
#   FFTW3_LIBRARIES    - Library to link against
#
# Imported target:
#   FFTW3::fftw3
#
# Hint: FFTW3_ROOT or FFTW_ROOT

find_path(FFTW3_INCLUDE_DIR
  NAMES fftw3.h
  HINTS ${FFTW3_ROOT} ${FFTW_ROOT}
  PATH_SUFFIXES include
)

find_library(FFTW3_LIBRARY
  NAMES fftw3
  HINTS ${FFTW3_ROOT} ${FFTW_ROOT}
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3
  REQUIRED_VARS FFTW3_LIBRARY FFTW3_INCLUDE_DIR
)

if(FFTW3_FOUND)
  set(FFTW3_INCLUDE_DIRS ${FFTW3_INCLUDE_DIR})
  set(FFTW3_LIBRARIES    ${FFTW3_LIBRARY})
  if(NOT TARGET FFTW3::fftw3)
    add_library(FFTW3::fftw3 UNKNOWN IMPORTED)
    set_target_properties(FFTW3::fftw3 PROPERTIES
      IMPORTED_LOCATION             "${FFTW3_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${FFTW3_INCLUDE_DIRS}"
    )
  endif()
endif()

mark_as_advanced(FFTW3_INCLUDE_DIR FFTW3_LIBRARY)
