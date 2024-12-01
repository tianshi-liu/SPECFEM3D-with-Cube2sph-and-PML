set(CMAKE_Fortran_COMPILER "/scinet/niagara/software/2019b/opt/intel-2019u4/openmpi/4.0.1-hpcx2.5/bin/mpif90")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_Fortran_COMPILER_VERSION "19.0.4.20190416")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "ELF")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/scinet/niagara/software/2019b/opt/intel-2019u4/openmpi/4.0.1-hpcx2.5/include;/scinet/niagara/software/2019b/opt/intel-2019u4/openmpi/4.0.1-hpcx2.5/lib;/scinet/niagara/software/2019b/opt/intel-2019u4-hdf5-1.8.21/netcdf/4.6.3/include;/scinet/niagara/software/2019b/opt/intel-2019u4/hdf5/1.8.21/include;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/ipp/include;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/mkl/include;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/pstl/include;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/tbb/include;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/daal/include;/gpfs/fs1/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/compiler/include/intel64;/gpfs/fs1/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/compiler/include/icc;/gpfs/fs1/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/compiler/include;/usr/local/include;/gpfs/fs1/scinet/niagara/software/2019b/core/lib/gcc/x86_64-pc-linux-gnu/7.4.0/include;/gpfs/fs1/scinet/niagara/software/2019b/core/lib/gcc/x86_64-pc-linux-gnu/7.4.0/include-fixed;/gpfs/fs1/scinet/niagara/software/2019b/core/include;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "mpi_usempif08;mpi_usempi_ignore_tkr;mpi_mpifh;mpi;ifport;ifcoremt;imf;svml;m;ipgo;irc;pthread;svml;c;gcc;gcc_s;irc_s;dl;c")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/scinet/niagara/software/2019b/opt/intel-2019u4/openmpi/4.0.1-hpcx2.5/lib;/scinet/niagara/software/2019b/opt/intel-2019u4-hdf5-1.8.21/netcdf/4.6.3/lib;/scinet/niagara/software/2019b/opt/intel-2019u4/hdf5/1.8.21/lib;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/ipp/lib/intel64;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/mkl/lib/intel64_lin;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/tbb/lib/intel64/gcc4.1;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/daal/lib/intel64_lin;/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/tbb/lib/intel64_lin/gcc4.4;/gpfs/fs1/scinet/intel/2019u4/compilers_and_libraries_2019.4.243/linux/compiler/lib/intel64_lin;/gpfs/fs1/scinet/niagara/software/2019b/core/lib/gcc/x86_64-pc-linux-gnu/7.4.0;/gpfs/fs1/scinet/niagara/software/2019b/core/lib/gcc;/gpfs/fs1/scinet/niagara/software/2019b/core/lib64;/lib64;/usr/lib64;/gpfs/fs1/scinet/niagara/software/2019b/core/lib;/lib;/usr/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
