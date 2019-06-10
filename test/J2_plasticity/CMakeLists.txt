add_executable(test_J2_plasticity test_J2_plasticity.cc)
target_link_libraries(test_J2_plasticity PRIVATE ConstitutiveModel)

add_test(
  NAME
     test_J2_plasticity
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/../run.sh test_J2_plasticity J2_plasticity stress)
