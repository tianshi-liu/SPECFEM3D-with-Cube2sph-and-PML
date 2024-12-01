set(SOURCE "convolve_source_timefunction.f90" 
        "detect_duplicates_stations_file.f90"
        "combine_vol_data_adios_stubs.f90"
        "combine_vol_data_impl.f90"
        "combine_vol_data_vtk_binary.F90"
        "vtk_writer.c")
add_library(aux ${SOURCE})
add_dependencies(aux solver shared)

add_executable(xcombine_surf_data combine_surf_data.f90)
add_executable(xcreate_slice_loc create_slice_loc.f90)
add_executable(xcreate_slice create_slice.f90)
add_executable(xcombine_vol_data_vtk combine_vol_data_vtk_binary.F90)
add_executable(xcreate_movie_shakemap_AVS_DX_GMT create_movie_shakemap_AVS_DX_GMT.f90)
target_link_libraries(xcombine_surf_data aux database mesh shared)
target_link_libraries(xcreate_slice_loc aux solver database mesh shared)
target_link_libraries(xcreate_slice solver aux database mesh shared)
target_link_libraries(xcombine_vol_data_vtk aux solver database mesh shared)
target_link_libraries(xcreate_movie_shakemap_AVS_DX_GMT solver database mesh shared)
