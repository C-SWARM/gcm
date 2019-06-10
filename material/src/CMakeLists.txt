add_library(material_obj
  OBJECT
    material_properties.cc
    slip_system.cc)

target_include_directories(material_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
  PRIVATE
    ${CMAKE_SOURCE_DIR}/local_include
    ${CMAKE_SOURCE_DIR}/utils/include
   )
