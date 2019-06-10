add_executable(test_elasticity test_elasticity.cc)
target_link_libraries(test_elasticity PRIVATE ConstitutiveModel)

add_test(
  NAME
    test_elasticity
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/../run.sh test_elasticity elasticity stress)
