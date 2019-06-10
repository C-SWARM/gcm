cmake_minimum_required(VERSION 3.13)

project(GCM LANGUAGES CXX)

# tell cmake we want to use c++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#
include(CTest)

set(BLA_VENDOR Intel10_64lp_seq)
find_package(BLAS REQUIRED)

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
