cmake_minimum_required(VERSION 3.5.0)
project(potpourri3d_bindings)

# Recurse in to pybind
set(PYBIND11_NEWPYTHON ON)
add_subdirectory(deps/pybind11)

# set location of eigen for geometry-central
set(GC_EIGEN_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/deps/eigen" CACHE PATH "my path")

# geometry-central
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_subdirectory(deps/geometry-central)

pybind11_add_module(potpourri3d_bindings
  src/cpp/core.cpp
  src/cpp/io.cpp
  src/cpp/heat_helpers.cpp
  src/cpp/mesh.cpp
  src/cpp/point_cloud.cpp
)

include_directories(potpourri3d_bindings ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp)

target_link_libraries(potpourri3d_bindings PRIVATE geometry-central)

install(TARGETS potpourri3d_bindings LIBRARY DESTINATION .)