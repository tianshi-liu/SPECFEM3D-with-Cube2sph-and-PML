#  CMakeLists.txt

cmake_minimum_required(VERSION 3.10.0)

################## compiler options ##################
# compiler
set(CMAKE_C_COMPILER          ${CC})
set(CMAKE_Fortran_COMPILER     ${MPIFC})


# build function test
set(BUILD_TEST                  ${BUILD_TEST})

################## dependencies ##################
# preferred installation prefix of dependencies

# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
set(${var} ${val})
endif()
endmacro()

# project
project(cube2sph LANGUAGES C Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/obj)

################# cmake setup #################
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_C_COMPILER_ID STREQUAL "GNU")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -march=native")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=gnu -g -march=native ")
elseif (CMAKE_C_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -xHost -traceback")
else()
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -march=native")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=gnu -g -march=native ")
endif()

set(CMAKE_VERBOSE_MAKEFILE off)

################# find packages #################
# defination
add_definitions(-DWITH_MPI)
add_definitions(-DFORCE_VECTORIZATION)

# mpi
find_package(MPI REQUIRED COMPONENTS C Fortran)
include_directories(${MPI_C_INCLUDE_DIRS})

# netcdf
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/include)
find_package(NetCDF COMPONENTS Fortran)
include_directories(${NetCDF_Fortran_INCLUDE_DIR})
message(${NetCDF_Fortran_INCLUDE_DIR})

################# include paths #################
include_directories(${PROJECT_BINARY_DIR}/obj  ${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR})

file(GLOB_RECURSE SOURCE CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE "src/*.c" "src/*.f90" "src/*.F90" )
set(BINARY_LIST "cube2sph_adepml" "cube2sph_force"  "cube2sph_station" "setup_adepml" 
                "write_force_solution_file"   "write_stations_file"
                "cube2sph_boundary_gmt"  "cube2sph_nopml"  "cube2sph_station_inv"  
                "write_cmt_solution_file"  "write_horizontal_slice_vtk"  "write_vertical_slice_vtk"
                "cube2latlon"
)
foreach(f  ${BINARY_LIST})
    list(FILTER SOURCE EXCLUDE REGEX ${f})
endforeach()

# add library
add_library(cube2sph ${SOURCE})

foreach(f ${BINARY_LIST})
    add_executable(${f} src/${f}.f90)
    target_link_libraries(${f} cube2sph ${NetCDF_Fortran_LIBRARIES})
endforeach()