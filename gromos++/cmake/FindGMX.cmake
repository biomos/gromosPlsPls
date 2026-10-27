# FindGMX.cmake - Find the Gromacs library
#
# Variables set:
#   GMX_FOUND       - True if found
#   GMX_LIBRARIES   - Library to link against
#
# Imported target:
#   GMX::gmx
#
# Hint: GROMACS_ROOT or GMX_ROOT

find_path(GMX_INCLUDE_DIR
  NAMES gromacs/version.h gromacs/utility/fatalerror.h
  HINTS ${GROMACS_ROOT} ${GMX_ROOT}
  PATH_SUFFIXES include
)

find_library(GMX_LIBRARY
  NAMES gmx
  HINTS ${GROMACS_ROOT} ${GMX_ROOT}
  PATH_SUFFIXES lib
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMX REQUIRED_VARS GMX_LIBRARY)

if(GMX_FOUND)
  set(GMX_LIBRARIES ${GMX_LIBRARY})
  if(NOT TARGET GMX::gmx)
    add_library(GMX::gmx UNKNOWN IMPORTED)
    set_target_properties(GMX::gmx PROPERTIES
      IMPORTED_LOCATION "${GMX_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${GMX_INCLUDE_DIR}"
    )
  endif()
endif()

mark_as_advanced(GMX_INCLUDE_DIR GMX_LIBRARY)
