#  CMakeLists.txt

cmake_minimum_required(VERSION 3.10.0)

################## compiler options ##################
# compiler
set(CMAKE_C_COMPILER          ${CC})
set(CMAKE_Fortran_COMPILER     ${FC})


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
project(JDSurfG LANGUAGES C Fortran)

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_SOURCE_DIR}/obj)

################# cmake setup #################
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -march=native -ffast-math")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=gnu -g -march=native -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused")
elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -xHost")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std08 -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence  -diag-disable 6477,8889 -implicitnone -gen-interfaces -warn all -traceback -O3 -check nobounds")
else()
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -march=native -ffast-math")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=gnu -g -march=native -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused")
endif()

set(CMAKE_VERBOSE_MAKEFILE off)

################# find packages #################
# defination
add_definitions(-DWITH_MPI)
add_definitions(-DFORCE_VECTORIZATION)

# mpi
find_package(MPI REQUIRED COMPONENTS C Fortran)
include_directories(${MPI_C_INCLUDE_DIRS})

################# include paths #################
include_directories(${PROJECT_SOURCE_DIR}/obj  ${PROJECT_SOURCE_DIR}/setup)

# subs
add_subdirectory(${PROJECT_SOURCE_DIR}/src/shared)
add_subdirectory(${PROJECT_SOURCE_DIR}/src/meshfem3D) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/cube2sph) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/generate_databases) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/cuda) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/specfem3D) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/auxiliaries) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/decompose_mesh) 
add_subdirectory(${PROJECT_SOURCE_DIR}/src/tomography) 