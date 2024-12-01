file(GLOB SOURCE CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE "*.c" "*.f90" "*.F90" )
list(FILTER SOURCE EXCLUDE REGEX "adios_helper")
list(FILTER SOURCE EXCLUDE REGEX "adios_manager.F90")
list(FILTER SOURCE EXCLUDE REGEX "asdf_manager.f90")

add_library(shared ${SOURCE})