add_library(poro_visco_plasticity_obj
  OBJECT
    poro_visco_plasticity.cc)

target_include_directories(poro_visco_plasticity_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
  PRIVATE
    ${CMAKE_SOURCE_DIR}/local_include
    ${CMAKE_SOURCE_DIR}/utils/include
    ${CMAKE_SOURCE_DIR}/constitutive_model_handle/include
    ${CMAKE_SOURCE_DIR}/elasticity/include
    ${CMAKE_SOURCE_DIR}/material/include
    )
