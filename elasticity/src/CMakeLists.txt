add_library(elasticity_obj
  OBJECT
    hyperelasticity.cc
    strain_energy_density_function.cc)

target_include_directories(elasticity_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
  PRIVATE
    ${CMAKE_SOURCE_DIR}/local_include
    ${CMAKE_SOURCE_DIR}/utils/include
    ${CMAKE_SOURCE_DIR}/material/include
    )
