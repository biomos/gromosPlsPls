add_library(gmath OBJECT
  Vec.cc
  Matrix.cc
  Distribution.cc
  WDistribution.cc
  Expression.cc
  Correlation.cc
  Physics.cc
  Mesh.cc
)
target_include_directories(gmath PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${GSL_INCLUDE_DIRS}
)
