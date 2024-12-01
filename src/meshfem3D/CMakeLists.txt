file(GLOB SOURCE CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE "*.c" "*.f90" "*.F90" )
list(FILTER SOURCE EXCLUDE REGEX "adios.F90")
list(FILTER SOURCE EXCLUDE REGEX "meshfem3D.F90")

# library
add_library(mesh ${SOURCE})
add_dependencies(mesh shared)

add_executable(xmeshfem3D meshfem3D.F90)
target_link_libraries(xmeshfem3D mesh shared)
