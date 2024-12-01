file(GLOB SOURCE CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE "*.c" "*.f90" "*.F90" )
list(FILTER SOURCE EXCLUDE REGEX "adios.F90")
list(FILTER SOURCE EXCLUDE REGEX "^generate_databases.f90")

# library
add_library(database ${SOURCE})
add_dependencies(database cube2sph_rot)

add_executable(xgenerate_databases generate_databases.f90)
target_link_libraries(xgenerate_databases database cube2sph_rot mesh shared)
