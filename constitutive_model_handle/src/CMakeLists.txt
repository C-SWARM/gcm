add_library(constitutive_model_handle_obj
  OBJECT
    GcmSolverInfo.cc
    constitutive_model_handle.cc)

target_include_directories(constitutive_model_handle_obj
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
  PRIVATE
    ${CMAKE_SOURCE_DIR}/local_include
    ${CMAKE_SOURCE_DIR}/utils/include
    ${CMAKE_SOURCE_DIR}/elasticity/include
    ${CMAKE_SOURCE_DIR}/material/include
    )
