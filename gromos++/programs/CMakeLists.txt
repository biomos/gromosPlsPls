# All programs include src/ for gromos++ headers
include_directories(${CMAKE_SOURCE_DIR}/src)

# Standard programs — one .cc per executable, all link gromos
set(GROMOS_PROGRAMS
  rmsd
  pdb2g96
  ener
  dg_ener
  noe
  prep_noe
  post_noe
  prep_bb
  rgyr
  tser
  ditrans
  diffus
  tstrip
  temperature
  gathtraj
  edyn
  rmsf
  dipole
  rmsdmat
  cluster
  postcluster
  check_box
  make_top
  make_sasa_top
  con_top
  com_top
  red_top
  pt_top
  link_top
  copy_box
  build_box
  bin_box
  ran_box
  ran_solvation
  pert_top
  explode
  hbond
  frameout
  nhoparam
  rdf
  mdf
  iondens
  mk_script
  check_top
  ene_ana
  int_ener
  visco
  epsilon
  eps_field
  gch
  sim_box
  ion
  gca
  dssp
  atominfo
  shake_analysis
  filter
  follow
  jval
  jepot
  swd
  sasa
  sasa_hasel
  tcf
  m_widom
  cog
  rot_rel
  inbox
  unify_box
  make_pt_top
  solute_entropy
  pairlist
  prep_eds
  eds_update_1
  eds_update_2
  dfmult
  reweight
  ext_ti_ana
  ext_ti_merge
  bar
  bilayer_dist
  bilayer_oparam
  rep_ana
  matrix_overlap
  close_pair
  predict_noe
  prep_xray_le
  xrayts
  duplicate
  epath
  disicl
  trs_ana
  dGslv_pbsolv
  dfgrid
  contactnum
  cos_dipole
  cos_epsilon
  addvirt_top
  amber2gromos
)

foreach(prog ${GROMOS_PROGRAMS})
  add_executable(${prog} ${prog}.cc)
  target_link_libraries(${prog} PRIVATE gromos)
endforeach()

# mk_script.h is also installed as a public header
install(FILES mk_script.h DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

# Clipper-dependent x-ray programs
set(CLIPPER_PROGRAMS
  cry
  structure_factor
  prep_xray
  r_factor
  xray_map
  r_real_factor
  cry_rms
)
if(CLIPPER_FOUND)
  foreach(prog ${CLIPPER_PROGRAMS})
    add_executable(${prog} ${prog}.cc)
    target_link_libraries(${prog} PRIVATE gromos Clipper::clipper)
  endforeach()
endif()

# Install all programs
install(TARGETS ${GROMOS_PROGRAMS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
if(CLIPPER_FOUND)
  install(TARGETS ${CLIPPER_PROGRAMS} RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
