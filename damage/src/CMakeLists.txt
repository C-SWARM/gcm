add_library(damage_obj
  OBJECT
    continuum_damage_model.cc)

target_include_directories(damage_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
  PRIVATE
    ${CMAKE_SOURCE_DIR}/local_include
    ${CMAKE_SOURCE_DIR}/utils/include
    ${CMAKE_SOURCE_DIR}/constitutive_model_handle/include
    ${CMAKE_SOURCE_DIR}/elasticity/include
    ${CMAKE_SOURCE_DIR}/material/include
    )
