add_library(utils_obj
  OBJECT
    input_file_read_handles.cc
    math_help.cc)

target_include_directories(utils_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include)

target_link_libraries(utils_obj PUBLIC lapack)
