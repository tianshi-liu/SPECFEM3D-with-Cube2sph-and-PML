set(POST_DIR "postprocess_sensitivity_kernels")
file(GLOB SOURCE CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE "${POST_DIR}/postprocess_par.f90" )

# library
add_library(tomo ${SOURCE})
add_dependencies(tomo shared)

add_executable(xsmooth_sem_sph_pde ${POST_DIR}/smooth_sem_sph_pde.F90)
target_link_libraries(xsmooth_sem_sph_pde tomo solver shared)
