add_library(data_structure_obj
  OBJECT
    data_structure.cc)

target_include_directories(data_structure_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    )
