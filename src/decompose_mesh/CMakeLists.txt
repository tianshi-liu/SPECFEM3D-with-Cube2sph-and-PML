file(GLOB SOURCE CONFIGURE_DEPENDS LIST_DIRECTORIES FALSE "*.f90" "*.F90" )
list(FILTER SOURCE EXCLUDE REGEX "program_decompose")

set(SCOTCHDIR ${PROJECT_SOURCE_DIR}/external_libs/scotch)
file(GLOB SOURCE1 "${SCOTCHDIR}/src/libscotch/*.c" "${SCOTCHDIR}/src/libscotch/last_resort/*.c")
list(FILTER SOURCE1 EXCLUDE REGEX "dgraph")
# list(FILTER SOURCE1 EXCLUDE REGEX "library_error.c")
list(FILTER SOURCE1 EXCLUDE REGEX "comm.c")
list(FILTER SOURCE1 EXCLUDE REGEX "common_sort.c")
list(FILTER SOURCE1 EXCLUDE REGEX "dorder")
list(FILTER SOURCE1 EXCLUDE REGEX "dmapping")
list(FILTER SOURCE1 EXCLUDE REGEX "edge.c")

# library
add_library(decompose ${SOURCE})
add_dependencies(decompose mesh)
add_definitions(-DUSE_SCOTCH)

add_library(scotch ${SOURCE1})

INCLUDE_DIRECTORIES(${SCOTCHDIR}/include )
INCLUDE_DIRECTORIES(${SCOTCHDIR}/src/libscotch )
target_compile_options(scotch PRIVATE -DCOMMON_FILE_COMPRESS_GZ -DCOMMON_PTHREAD -DCOMMON_RANDOM_FIXED_SEED -DSCOTCH_RENAME -DSCOTCH_RENAME_PARSER -DSCOTCH_PTHREAD -DSCOTCH_COMMON_RENAME -Drestrict=__restrict -pthread)

add_executable(xdecompose_mesh program_decompose_mesh.f90)
target_link_libraries(xdecompose_mesh decompose scotch  mesh shared)
