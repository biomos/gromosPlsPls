add_library(args OBJECT
  Arguments.cc
  BoundaryParser.cc
  ReferenceParser.cc
  GatherParser.cc
  OutformatParser.cc
)
target_include_directories(args PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
)
