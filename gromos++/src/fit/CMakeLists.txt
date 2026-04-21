add_library(fit OBJECT
  PositionUtils.cc
  Reference.cc
  RotationalFit.cc
  TranslationalFit.cc
  FastRotationalFit.cc
  AtomDistances.cc
)
target_include_directories(fit PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${GSL_INCLUDE_DIRS}
)
