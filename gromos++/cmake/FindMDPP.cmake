# FindMDPP.cmake - Find the GROMOS MD++ library
#
# Variables set:
#   MDPP_FOUND       - True if found
#   MDPP_LIBRARIES   - Library to link against
#
# Imported target:
#   MDPP::mdpp
#
# Hint: MDPP_ROOT

find_path(MDPP_INCLUDE_DIR
  NAMES md++/math/gmath.h
  HINTS ${MDPP_ROOT}
  PATH_SUFFIXES include
)

find_library(MDPP_LIBRARY
  NAMES mdpp
  HINTS ${MDPP_ROOT}
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MDPP REQUIRED_VARS MDPP_LIBRARY)

if(MDPP_FOUND)
  set(MDPP_LIBRARIES ${MDPP_LIBRARY})
  if(NOT TARGET MDPP::mdpp)
    add_library(MDPP::mdpp UNKNOWN IMPORTED)
    set_target_properties(MDPP::mdpp PROPERTIES
      IMPORTED_LOCATION "${MDPP_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${MDPP_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(MDPP_INCLUDE_DIR MDPP_LIBRARY)
