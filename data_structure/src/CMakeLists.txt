add_library(data_structure_obj INTERFACE)

target_include_directories(data_structure_obj
  INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_link_libraries(data_structure_obj INTERFACE lapack)
