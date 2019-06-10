add_executable(test_crystal_plasticity_with_thermal test_crystal_plasticity_with_thermal.cc)
target_link_libraries(test_crystal_plasticity_with_thermal PRIVATE ConstitutiveModel)

add_test(
  NAME
     test_crystal_plasticity_with_thermal
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/test_crystal_plasticity_with_thermal.sh)
