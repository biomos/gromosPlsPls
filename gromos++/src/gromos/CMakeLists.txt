# Exception.cc only includes Exception.h; compiled into gromos as gromos_exc OBJECT lib
add_library(gromos_exc OBJECT Exception.cc)
target_include_directories(gromos_exc PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
)
