add_executable(test_damage_model test_damage_model.cc F_of_t.cc)
target_link_libraries(test_damage_model PRIVATE ConstitutiveModel)

add_test(
  NAME
    test_damage_model
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/../run.sh "test_damage_model ${CMAKE_CURRENT_SOURCE_DIR}/params_compression.in" damage stress.sigma)
