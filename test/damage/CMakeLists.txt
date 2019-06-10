add_executable(test_damage_model test_damage_model.cc F_of_t.cc)
target_link_libraries(test_damage_model PRIVATE ConstitutiveModel)

add_test(
  NAME test_damage_model
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/../runTest.sh
    "${CMAKE_CURRENT_BINARY_DIR}/test_damage_model ${CMAKE_CURRENT_SOURCE_DIR}/params_compression.in"
    ${CMAKE_CURRENT_SOURCE_DIR}/stress.sigma.ref
    ${CMAKE_CURRENT_BINARY_DIR}/stress.sigma.txt)
