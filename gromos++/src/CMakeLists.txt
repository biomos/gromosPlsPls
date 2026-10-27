add_subdirectory(gromos)
add_subdirectory(gmath)
add_subdirectory(gcore)
add_subdirectory(gio)
add_subdirectory(bound)
add_subdirectory(utils)
add_subdirectory(fit)
add_subdirectory(args)
add_subdirectory(pb)

# Combined gromos library built from all component OBJECT libraries
add_library(gromos
  $<TARGET_OBJECTS:gromos_exc>
  $<TARGET_OBJECTS:gmath>
  $<TARGET_OBJECTS:gcore>
  $<TARGET_OBJECTS:gio>
  $<TARGET_OBJECTS:bound>
  $<TARGET_OBJECTS:utils>
  $<TARGET_OBJECTS:fit>
  $<TARGET_OBJECTS:args>
  $<TARGET_OBJECTS:pb>
)
target_include_directories(gromos PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(gromos PUBLIC
  GSL::gsl
  GSL::gslcblas
  ZLIB::ZLIB
  $<$<BOOL:${MATH_LIBRARY}>:${MATH_LIBRARY}>
  $<$<BOOL:${ENABLE_OPENMP}>:OpenMP::OpenMP_CXX>
)
if(READLINE_LIBRARY AND READLINE_INCLUDE_DIR)
  target_link_libraries(gromos PUBLIC ${READLINE_LIBRARY})
  target_include_directories(gromos PUBLIC ${READLINE_INCLUDE_DIR})
endif()
if(FFTW3_FOUND)
  target_link_libraries(gromos PUBLIC FFTW3::fftw3)
endif()
set_target_properties(gromos PROPERTIES
  VERSION   "${GROMOS_MAJOR_VERSION}.${GROMOS_MICRO_VERSION}.0"
  SOVERSION "${GROMOS_MAJOR_VERSION}"
  OUTPUT_NAME gromos
)

# Dynamically loadable plugin library (same content as gromos, built as MODULE)
add_library(gromosplugin MODULE
  $<TARGET_OBJECTS:gromos_exc>
  $<TARGET_OBJECTS:gmath>
  $<TARGET_OBJECTS:gcore>
  $<TARGET_OBJECTS:gio>
  $<TARGET_OBJECTS:bound>
  $<TARGET_OBJECTS:utils>
  $<TARGET_OBJECTS:fit>
  $<TARGET_OBJECTS:args>
  $<TARGET_OBJECTS:pb>
)
target_link_libraries(gromosplugin PRIVATE
  GSL::gsl GSL::gslcblas ZLIB::ZLIB
  $<$<BOOL:${MATH_LIBRARY}>:${MATH_LIBRARY}>
)

install(TARGETS gromos gromosplugin
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install gromos/Exception.h to flat includedir (matches automake behaviour)
install(FILES gromos/Exception.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Install module headers to includedir/gromos++
install(FILES
  gmath/Vec.h gmath/Vec.cc gmath/Matrix.h gmath/Distribution.h gmath/WDistribution.h
  gmath/Stat.h gmath/Stat.cc gmath/StatDisk.h gmath/StatDisk.cc
  gmath/Expression.h gmath/Correlation.h gmath/Physics.h gmath/Mesh.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gromos++
)
install(FILES
  gcore/Angle.h gcore/AngleType.h gcore/AtomPair.h gcore/AtomTopology.h
  gcore/MassType.h gcore/Bond.h gcore/BondType.h gcore/Box.h gcore/Remd.h
  gcore/Constraint.h gcore/CrossDihedral.h gcore/Dihedral.h gcore/DihedralType.h
  gcore/Exclusion.h gcore/GromosForceField.h gcore/Improper.h gcore/ImproperType.h
  gcore/LJType.h gcore/LJException.h gcore/LJExceptionType.h gcore/CGType.h
  gcore/MoleculeTopology.h gcore/Molecule.h gcore/SolventTopology.h gcore/Solvent.h
  gcore/BuildingBlock.h gcore/BbSolute.h gcore/BbLink.h gcore/System.h
  gcore/LinearTopology.h gcore/PtTopology.h gcore/VirtualAtomType.h
  gcore/VirtualAtoms.h gcore/Weight.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gromos++
)
install(FILES
  gio/Ginstream.h gio/gzstream.h gio/InTopology.h gio/OutTopology.h
  gio/OutGromacs.h gio/InG96.h gio/OutCoordinates.h gio/OutG96.h gio/OutG96S.h
  gio/OutPdb.h gio/Outvmdam.h gio/OutCif.h gio/InParameter.h
  gio/InBuildingBlock.h gio/InLinkingBlock.h gio/InPtTopology.h
  gio/OutPtTopology.h gio/StringTokenizer.h gio/InCIF.h
  gio/InIACElementNameMapping.h gio/InBFactorOccupancy.h
  gio/OutBuildingBlock.h gio/InPDB.h gio/InAmberTopology.h gio/InChargeGroups.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gromos++
)
install(FILES
  bound/Boundary.h bound/Vacuum.h bound/TruncOct.h bound/RectBox.h bound/Triclinic.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gromos++
)
install(FILES
  utils/parse.h utils/Rmsd.h utils/Neighbours.h utils/VirtualAtom.h utils/Noe.h
  utils/AtomSpecifier.h utils/AtomicRadii.h utils/CremerPople.h
  utils/VectorSpecifier.h utils/PropertyContainer.h utils/Property.h
  utils/Energy.h utils/TrajArray.h utils/CheckTopo.h
  utils/Hbond_calc.h utils/Hbond_calc_2c.h utils/Hbond_calc_3c.h
  utils/Hbond_calc_bridges.h utils/Hbond.h utils/Dssp.h
  utils/EnergyTraj.h utils/RestrTraj.h utils/DipTraj.h
  utils/SimplePairlist.h utils/SoluteWeightedDistance.h
  utils/FfExpert.h utils/FfExpertGraph.h utils/make_top.h
  utils/Temperature.h utils/Value.h
  utils/ExpressionParser.h utils/ExpressionParser.cc
  utils/groTime.h utils/RdcFuncs.h utils/JvalFuncs.h utils/CommandLine.h
  utils/debug.h utils/AminoAcid.h utils/RDF.h utils/NeutronScattering.h
  utils/CubeSystem.hcc utils/Disicl.h utils/Gch.h utils/IntegerInputParser.h
  utils/StringOps.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gromos++
)
install(FILES
  fit/PositionUtils.h fit/Reference.h fit/RotationalFit.h
  fit/TranslationalFit.h fit/FastRotationalFit.h fit/AtomDistances.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gromos++
)
install(FILES
  args/Arguments.h args/BoundaryParser.h args/ReferenceParser.h
  args/GatherParser.h args/OutformatParser.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gromos++
)
install(FILES
  pb/FDPoissonBoltzmann.h pb/FDPoissonBoltzmann_ICCG_NPBC.h
  pb/FDPoissonBoltzmann_ICCG_PBC.h pb/PB_Parameters.h
  pb/FFTBoundaryCondition.h pb/FFTChargeDipole.h pb/FFTChargeDipole_LS.h
  pb/FFTChargeDipole_RF.h pb/FFTChargeShapingFunction.h
  pb/FFTChargeShapingFunction_Hat.h pb/FFTChargeShapingFunction_Parabola.h
  pb/FFTDipoleDipole.h pb/FFTDipoleDipole_LS.h pb/FFTDipoleDipole_RF.h
  pb/FFTGridType.h pb/FFTInsideOutside.h pb/FFTInteractionTypeCodes.h
  pb/FFTPoisson.h pb/FFTPoissonIterator.h pb/FFTVacuumField.h
  pb/FFTVacuumField_LS.h pb/FFTVacuumField_RF.h pb/Ewald_edir.h
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gromos++
)

# Tests - defined after gromos so they can link against it
function(add_gromos_test subdir name)
  add_executable(${subdir}_${name}
    ${CMAKE_CURRENT_SOURCE_DIR}/${subdir}/${name}.t.cc
  )
  target_link_libraries(${subdir}_${name} PRIVATE gromos)
  add_test(NAME ${subdir}_${name} COMMAND ${subdir}_${name})
endfunction()

# gmath tests
add_gromos_test(gmath Vec)
add_gromos_test(gmath Matrix)
add_gromos_test(gmath Distribution)
add_gromos_test(gmath Expression)
add_gromos_test(gmath Correlation)
add_gromos_test(gmath Mesh)

# gcore tests
add_gromos_test(gcore AtomTopology)
add_gromos_test(gcore Exclusion)
add_gromos_test(gcore MoleculeTopology)
add_gromos_test(gcore Molecule)
add_gromos_test(gcore SolventTopology)
add_gromos_test(gcore Solvent)
add_gromos_test(gcore BbSolute)

# gio tests
add_gromos_test(gio InG96)
add_gromos_test(gio InTopology)
add_gromos_test(gio InParameter)
add_gromos_test(gio InBuildingBlock)
add_gromos_test(gio OutGromacs)
add_gromos_test(gio OutG96S)
add_gromos_test(gio OutPdb)
add_gromos_test(gio Outvmdam)

# bound tests
add_gromos_test(bound TruncOct)
add_gromos_test(bound RectBox)

# utils tests
add_gromos_test(utils ExpressionParser)
add_gromos_test(utils Rmsd)
add_gromos_test(utils Neighbours)
add_gromos_test(utils VirtualAtom)
add_gromos_test(utils Noe)
add_gromos_test(utils AtomSpecifier)
add_gromos_test(utils PropertyContainer)
add_gromos_test(utils Energy)
add_gromos_test(utils CheckTopo)
add_gromos_test(utils SimplePairlist)
add_gromos_test(utils FfExpert)

# fit tests
add_gromos_test(fit PositionUtils)
add_gromos_test(fit Reference)
add_gromos_test(fit RotationalFit)
add_gromos_test(fit TranslationalFit)

# args test
add_gromos_test(args Arguments)

# gromos exception test
add_gromos_test(gromos Exception)
