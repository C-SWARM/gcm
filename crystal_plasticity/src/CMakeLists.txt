add_library(crystal_plasticity_obj
  OBJECT
    construct_linearization_parts.cc
    flowlaw.cc
    hardening.cc
    solve_system.cc)

target_include_directories(crystal_plasticity_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include  
  PRIVATE
    ${CMAKE_SOURCE_DIR}/local_include
    ${CMAKE_SOURCE_DIR}/constitutive_model_handle/include
    ${CMAKE_SOURCE_DIR}/elasticity/include
    ${CMAKE_SOURCE_DIR}/material/include
    ${CMAKE_SOURCE_DIR}/utils/include
    )
