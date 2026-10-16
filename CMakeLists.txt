cmake_minimum_required(VERSION 3.20)
project(demo_land_model LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  # Debug by default; configure with -DCMAKE_BUILD_TYPE=Release for production runs.
  set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Core model library: everything under src/ except the driver.
add_library(landmodel
  src/kinds_mod.f90
  src/constants_mod.f90
  src/types_mod.f90
  src/tridiagonal_mod.f90
  src/soil_properties_mod.f90
  src/soil_water_mod.f90
  src/soil_heat_mod.f90
  src/canopy_mod.f90
  src/forcing_mod.f90
  src/land_mod.f90
  src/io_mod.f90
)

target_compile_options(landmodel PUBLIC
  -std=f2008 -ffree-line-length-none -Wall -Wextra
  $<$<CONFIG:Debug>:-O0 -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow>
  $<$<CONFIG:Release>:-O3>
)

add_executable(land_driver app/main.f90)
target_link_libraries(land_driver PRIVATE landmodel)

# Tests
enable_testing()
add_library(test_utilities tests/test_utilities_mod.f90)
target_link_libraries(test_utilities PUBLIC landmodel)

foreach(test_name tridiagonal soil_water soil_heat canopy land)
  add_executable(test_${test_name} tests/test_${test_name}.f90)
  target_link_libraries(test_${test_name} PRIVATE test_utilities)
  add_test(NAME ${test_name} COMMAND test_${test_name})
endforeach()
