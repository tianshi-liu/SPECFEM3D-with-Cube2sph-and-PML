file(GLOB SOURCE CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE "*.c" "*.f90" "*.F90" )
list(FILTER SOURCE EXCLUDE REGEX "adios.F90")
list(FILTER SOURCE EXCLUDE REGEX "ASDF.f90")
list(FILTER SOURCE EXCLUDE REGEX "specfem3D.F90")

# library
add_library(solver ${SOURCE})
add_dependencies(solver cuda database)

add_executable(xspecfem3D specfem3D.F90)
target_link_libraries(xspecfem3D solver cuda database cube2sph_rot mesh shared)
