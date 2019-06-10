add_executable(test_crystal_plasticity_single_crystal test_crystal_plasticity_single_crystal.cc)
target_link_libraries(test_crystal_plasticity_single_crystal PRIVATE ConstitutiveModel)


add_test(
  NAME test_crystal_plasticity_single_crystal
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/../runTest.sh
    ${CMAKE_CURRENT_BINARY_DIR}/test_crystal_plasticity_single_crystal
    ${CMAKE_CURRENT_SOURCE_DIR}/single_crystal_results.ref
    ${CMAKE_CURRENT_BINARY_DIR}/single_crystal_results.txt)


