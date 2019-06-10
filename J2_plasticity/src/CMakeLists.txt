add_library(J2_plasticity_obj
  OBJECT
    J2_plasticity.cc)

target_include_directories(J2_plasticity_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
  PRIVATE
    ${CMAKE_SOURCE_DIR}/local_include
    ${CMAKE_SOURCE_DIR}/constitutive_model_handle/include
    ${CMAKE_SOURCE_DIR}/elasticity/include
    ${CMAKE_SOURCE_DIR}/material/include
    ${CMAKE_SOURCE_DIR}/utils/include)

target_link_libraries(J2_plasticity_obj PRIVATE ttl)
