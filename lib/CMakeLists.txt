add_library(ConstitutiveModel
  $<TARGET_OBJECTS:J2_plasticity_obj>
  $<TARGET_OBJECTS:constitutive_model_handle_obj>
  $<TARGET_OBJECTS:crystal_plasticity_obj>
  $<TARGET_OBJECTS:damage_obj>
  $<TARGET_OBJECTS:elasticity_obj>
  $<TARGET_OBJECTS:material_obj>
  $<TARGET_OBJECTS:utils_obj>
  $<TARGET_OBJECTS:poro_visco_plasticity_obj>)

target_link_libraries(ConstitutiveModel
  INTERFACE
    J2_plasticity_obj
    constitutive_model_handle_obj
    crystal_plasticity_obj
    damage_obj
    elasticity_obj
    material_obj
    utils_obj
    poro_visco_plasticity_obj
    data_structure_obj)

target_include_directories(ConstitutiveModel
  INTERFACE
    ${CMAKE_SOURCE_DIR}/local_include
    ${CMAKE_SOURCE_DIR}/utils/include)
