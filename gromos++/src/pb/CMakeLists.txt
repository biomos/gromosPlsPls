add_library(pb OBJECT
  FDPoissonBoltzmann.cc
  FDPoissonBoltzmann_ICCG_NPBC.cc
  FDPoissonBoltzmann_ICCG_PBC.cc
  PB_Parameters.cc
  FFTBoundaryCondition.cc
  FFTChargeDipole.cc
  FFTChargeDipole_LS.cc
  FFTChargeDipole_RF.cc
  FFTChargeShapingFunction.cc
  FFTChargeShapingFunction_Hat.cc
  FFTChargeShapingFunction_Parabola.cc
  FFTDipoleDipole.cc
  FFTDipoleDipole_LS.cc
  FFTDipoleDipole_RF.cc
  FFTGridType.cc
  FFTInsideOutside.cc
  FFTInteractionTypeCodes.cc
  FFTPoisson.cc
  FFTPoissonIterator.cc
  FFTVacuumField.cc
  FFTVacuumField_LS.cc
  FFTVacuumField_RF.cc
  Ewald_edir.cc
)
target_include_directories(pb PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_SOURCE_DIR}/src
  ${GSL_INCLUDE_DIRS}
)
if(FFTW3_FOUND)
  target_include_directories(pb PRIVATE ${FFTW3_INCLUDE_DIRS})
endif()
