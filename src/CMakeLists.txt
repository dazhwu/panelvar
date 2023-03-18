cmake_minimum_required(VERSION 3.1)
##cd build
##cmake ..
##make install

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#add_definitions(-DEIGEN_USE_LAPACKE)
project(pvar_module VERSION 0.1.0)


set(CMAKE_CXX_FLAGS "-Wall -O3 -m64 -fopenmp  -mavx -mfma  -march=native ")
#-m64 -fopenmp -march=native  -mavx2 -DNDEBUG


# Include dir
include_directories(/usr/local/include)
include_directories(/usr/include/eigen3)
include_directories(/usr/include/eigen3/unsupported)

#include_directories(/usr/include/python3.9)
#include_directories(/usr/local/lib/python3.9/dist-packages/pybind11/include)
include_directories(/usr/include/python3.10)
include_directories(/usr/local/lib/python3.10/dist-packages/pybind11/include)

set(LIB_OUT "lib")

# Source files

aux_source_directory(. DIR_SRCS)



add_library(pvar_module SHARED ${DIR_SRCS})


set_target_properties(pvar_module PROPERTIES PREFIX "")

set_target_properties(pvar_module PROPERTIES OUTPUT_NAME "pvar_module")

install(TARGETS pvar_module DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../panelvar/)
