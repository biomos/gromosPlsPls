include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/programs)

# VMD plugin convenience library (noinst equivalent)
if(VMD_INCLUDE_DIR)
  add_library(vmdplugin STATIC vmdplugin.cc)
  target_include_directories(vmdplugin PRIVATE ${VMD_INCLUDE_DIR})
endif()

# Standard contrib programs — one .cc per executable, all link gromos
set(CONTRIB_PROGRAMS
  aggregates
  rottemp
  ener_rplc
  espmap
  lb_top
  m_ener
  cg_potential
  countbb
  dipole_pot
  stacking
  model
  gca_ener
  renumber
  distmat
  iwdcf
  bb2tex
  crossdih_top
  saltbridge
  bilayer_gel
  make_rdc_spec
  fit_rdc
  add_atom
  top_to_build
  helix_axis
  diff_top
  grid_dens
  pdb2seq
  nscatt
  structal
  cgLJpot
  hvap
  depsi2native
  rep_rewrite
  local_comp
  prep_hybrid
  sasa_new
  ener_test
  force
  mut_top
  mod_top
  rmsd_cmp
  eds_mult_all
  stat_ana
  sphericity
  expand_ua
  fit_ener
  fit_ener_traj
)

foreach(prog ${CONTRIB_PROGRAMS})
  add_executable(${prog} ${prog}.cc)
  target_link_libraries(${prog} PRIVATE gromos)
endforeach()

# MD++-dependent program — always built; MDPP linked only when available
add_executable(shake shake.cc)
target_link_libraries(shake PRIVATE gromos)
if(MDPP_FOUND)
  target_link_libraries(shake PRIVATE MDPP::mdpp)
  target_compile_definitions(shake PRIVATE HAVE_MDPP)
endif()

# Gromacs-dependent program — always built; GMX linked only when available
add_executable(gromacs2gromos gromacs2gromos.cc)
target_link_libraries(gromacs2gromos PRIVATE gromos)
if(GMX_FOUND)
  target_link_libraries(gromacs2gromos PRIVATE GMX::gmx)
  target_compile_definitions(gromacs2gromos PRIVATE HAVE_GMX)
endif()

install(TARGETS ${CONTRIB_PROGRAMS} shake gromacs2gromos
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
