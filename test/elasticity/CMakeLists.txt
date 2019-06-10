add_executable(test_elasticity test_elasticity.cc)
target_link_libraries(test_elasticity PRIVATE ConstitutiveModel)

add_test(
  NAME test_elasticity
  COMMAND
    ${CMAKE_CURRENT_SOURCE_DIR}/../runTest.sh
    ${CMAKE_CURRENT_BINARY_DIR}/test_elasticity
    ${CMAKE_CURRENT_SOURCE_DIR}/stress.ref
    ${CMAKE_CURRENT_BINARY_DIR}/stress.txt)
