cmake_minimum_required(VERSION 3.22)

project(Matrix)

add_definitions(-DWITH_AVX2)

add_definitions(-DWITH_AVX2)

add_definitions(-DINT)

set(CMAKE_C_STANDARD 11)

add_compile_options(-mfma -fopenmp -mavx2)

add_executable(matrix matrix.cpp matrix.hpp)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message("OpenMP found.")
    target_link_libraries(matrix PUBLIC OpenMP::OpenMP_CXX)
endif()