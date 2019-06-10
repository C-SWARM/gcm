cmake_minimum_required(VERSION 3.13)

project(GCM LANGUAGES CXX)

option(GCM_USE_MKL "Use Intel MKL for BLAS/LAPACK" ON)

# tell cmake we want to use c++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# enable cmake unit testing.
include(CTest)

# Deal with LAPACK/BLAS ... we will remove this dependency some day.
if (GCM_USE_MKL)
  set (BLA_VENDOR Intel10_64lp_seq)
endif (GCM_USE_MKL)
find_package(LAPACK)

add_library(lapack INTERFACE)
target_include_directories(lapack
  INTERFACE
    ${CMAKE_CURRENT_BINARY_DIR}/include)
target_link_libraries(lapack INTERFACE ${LAPACK_LIBRARIES})

# find ttl
include(FetchContent)

FetchContent_Declare(
  ttl
  GIT_REPOSITORY https://github.com/C-SWARM/ttl.git
  GIT_TAG develop
  UPDATE_DISCONNECTED 1
  FETCHCONTENT_QUIET OFF)

set(TTL_BUILD_TESTS OFF CACHE BOOL "build ttl tests" FORCE)

FetchContent_GetProperties(ttl)
if(NOT ttl_POPULATED)
  FetchContent_Populate(ttl)
  add_subdirectory(${ttl_SOURCE_DIR} ${ttl_BINARY_DIR})
endif()

add_subdirectory(include)
add_subdirectory(J2_plasticity)
add_subdirectory(constitutive_model_handle)
add_subdirectory(crystal_plasticity)
add_subdirectory(damage)
add_subdirectory(data_structure)
add_subdirectory(elasticity)
add_subdirectory(material)
add_subdirectory(utils)
add_subdirectory(poro_viscoplasticity)

add_subdirectory(lib)
add_subdirectory(test)
