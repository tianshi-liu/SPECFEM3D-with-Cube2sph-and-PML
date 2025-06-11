#  CMakeLists.txt

cmake_minimum_required(VERSION 3.20.0)

################## compiler options ##################
# compiler
set(CMAKE_C_COMPILER          ${CC})
set(CMAKE_Fortran_COMPILER     ${MPIFC})

# build function test
set(BUILD_TEST                  ${BUILD_TEST})
set(ENABLE_GPU              ${ENABLE_GPU})

################## dependencies ##################
# preferred installation prefix of dependencies

# set default values
macro(setDefault var val)
if(NOT DEFINED ${var})
  set(${var} ${val})
endif()
endmacro()

option(ENABLE_CUDA "Enable CUDA" OFF)
option(ENABLE_CUDA_AWARE  "Enable CUDA AWARE MPI" OFF)

# project
if(ENABLE_CUDA)
  project(specfem3d-cube2sph LANGUAGES C Fortran CUDA)
else()
  project(specfem3d-cube2sph LANGUAGES C Fortran)
endif()

# build type 
set(CMAKE_BUILD_TYPE RELEASE)
#set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED True)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/obj)

################# cmake setup #################
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
if (CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=gnu -g -march=native -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused")
elseif (CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -xHost")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std08 -xHost -fpe0 -ftz -assume buffered_io -assume byterecl -align sequence  -diag-disable 6477,8889,10441 -implicitnone -gen-interfaces -warn all -traceback -O3 -check nobounds")
else()
  set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -g -march=native")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=gnu -g -march=native -fimplicit-none -frange-check -fmax-errors=10 -pedantic -pedantic-errors -Waliasing -Wampersand -Wcharacter-truncation -Wline-truncation -Wsurprising -Wno-tabs -Wunderflow -ffpe-trap=invalid,zero,overflow -Wunused")
endif()

set(CMAKE_VERBOSE_MAKEFILE off)

# defination
add_definitions(-DWITH_MPI)
add_definitions(-DFORCE_VECTORIZATION)


################# find packages #################
# mpi
find_package(MPI REQUIRED COMPONENTS C Fortran)
include_directories(${MPI_C_INCLUDE_DIRS})

# cuda
if(ENABLE_CUDA) 
  find_package(CUDAToolkit REQUIRED)
  add_definitions(-DUSE_CUDA)
  include_directories(${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS}")
  set(GPULIBS CUDA::cudart)
  if(ENABLE_CUDA_AWARE) 
    message(STATUS "Enable CUDA-Aware MPI")
    add_definitions(-DUSE_CUDA_AWARE_MPI)
  endif()
endif()

################# include paths #################
include_directories(${CMAKE_BINARY_DIR}/obj  ${PROJECT_SOURCE_DIR}/setup)

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

if(ENABLE_CUDA) 
  # set your architecture here
  set_target_properties(cuda PROPERTIES CUDA_ARCHITECTURES native)
endif()
