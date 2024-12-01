file(GLOB SOURCE CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE "*.f90" )

add_library(cube2sph_rot ${SOURCE})
add_dependencies(cube2sph_rot shared)