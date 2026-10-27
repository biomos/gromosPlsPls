add_library(bound OBJECT
  Boundary.cc
  Vacuum.cc
  TruncOct.cc
  RectBox.cc
  Triclinic.cc
)
target_include_directories(bound PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
)
