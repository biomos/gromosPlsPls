# FindClipper.cmake - Find the Clipper x-ray crystallography library
#
# Variables set:
#   CLIPPER_FOUND        - True if found
#   CLIPPER_INCLUDE_DIRS - Include directories
#   CLIPPER_LIBRARIES    - Libraries to link against
#
# Imported target:
#   Clipper::clipper
#
# Hints: CLIPPER_ROOT, CCP4_ROOT

find_path(CLIPPER_INCLUDE_DIR
  NAMES clipper/clipper.h
  HINTS ${CLIPPER_ROOT} ${CCP4_ROOT}
  PATH_SUFFIXES include
)

find_library(CLIPPER_CORE_LIBRARY    NAMES clipper-core    HINTS ${CLIPPER_ROOT} ${CCP4_ROOT} PATH_SUFFIXES lib)
find_library(CLIPPER_CONTRIB_LIBRARY NAMES clipper-contrib HINTS ${CLIPPER_ROOT} ${CCP4_ROOT} PATH_SUFFIXES lib)
find_library(CLIPPER_CCP4_LIBRARY    NAMES clipper-ccp4    HINTS ${CLIPPER_ROOT} ${CCP4_ROOT} PATH_SUFFIXES lib)
find_library(CCP4C_LIBRARY           NAMES ccp4c           HINTS ${CCP4_ROOT}                 PATH_SUFFIXES lib)
find_library(RFFTW_LIBRARY           NAMES rfftw)
find_library(FFTW_LIBRARY            NAMES fftw)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Clipper
  REQUIRED_VARS CLIPPER_CORE_LIBRARY CLIPPER_INCLUDE_DIR
)

if(CLIPPER_FOUND)
  set(CLIPPER_INCLUDE_DIRS ${CLIPPER_INCLUDE_DIR})
  set(CLIPPER_LIBRARIES
    ${CLIPPER_CCP4_LIBRARY}
    ${CCP4C_LIBRARY}
    ${CLIPPER_CONTRIB_LIBRARY}
    ${CLIPPER_CORE_LIBRARY}
    ${RFFTW_LIBRARY}
    ${FFTW_LIBRARY}
    Threads::Threads
  )
  if(NOT TARGET Clipper::clipper)
    add_library(Clipper::clipper INTERFACE IMPORTED)
    set_target_properties(Clipper::clipper PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${CLIPPER_INCLUDE_DIRS}"
      INTERFACE_LINK_LIBRARIES      "${CLIPPER_LIBRARIES}"
    )
  endif()
endif()

mark_as_advanced(
  CLIPPER_INCLUDE_DIR CLIPPER_CORE_LIBRARY CLIPPER_CONTRIB_LIBRARY
  CLIPPER_CCP4_LIBRARY CCP4C_LIBRARY RFFTW_LIBRARY FFTW_LIBRARY
)
