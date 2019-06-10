add_executable(test_crystal_plasticity_single_crystal test_crystal_plasticity_single_crystal.cc)
target_link_libraries(test_crystal_plasticity_single_crystal PRIVATE ConstitutiveModel)

add_test(
  NAME
     test_crystal_plasticity_single_crystal
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/../run.sh test_crystal_plasticity_single_crystal crystal_plasticity_single_crystal single_crystal_results)


