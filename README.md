# Cube2sph

Accurate and flexible continental-scale seismic wave simulations based on SPECFEM3D

- mesh truncated at customized depth
- implementation of curvilinear PML based on auxiliary differential equations
- CUDA accelaration, on both forward/adjoint simulation and sensitivity kernels
- Teleseismic wavefield injection.

## Installation
This package and has two parts: a spectral-element solver and the cube2sph-toolkit.  

To compile the first part, you can use:
```bash 
mkdir -p build
cd build
cmake .. -DCC=gcc -DMPIFC=mpif90 
make -j8
```
If you want to build CUDA accelarated version, please use:
```bash 
cmake .. -DCC=gcc -DMPIFC=mpif90  -DENABLE_CUDA=ON
```
By default, the CUDA architecture will be chosen as `native`. If you want to use some user defined number, please open `CMakeLists.txt` and find `set_target_properties(cuda PROPERTIES CUDA_ARCHITECTURES native)`. Then set `native` to the number you want.

The CUDA-Aware MPI technique would facilitate communications. If you want to enable cuda aware mpi, please use:
```bash 
cmake .. -DCC=gcc -DMPIFC=mpif90 -DENABLE_CUDA=ON  \ 
        -DENABLE_CUDA_AWARE=ON
```


For the cube2sph toolkit, please make sure you've installed [netcdf-fortran](https://docs.unidata.ucar.edu/netcdf-fortran/current/) on you machine. Then you can go to `utils/cube2sph` and try the following:

```bash 
mkdir -p build
cd build
cmake .. -DCC=gcc -DMPIFC=mpif90 
make -j8
```

## Anisotropic Model Support
You should provide `c21` model with density in `tomography.xyz`.
1. Anisotropic models in `proc*_external_mesh.bin` are in cartesian coordinates, i.e. $C_{xxxx}, C_{xyzx}$ 
2. Anisotropic models in `proc*_c11-c66.bin`, `tomography.xyz` are in radial coordinates,i.e. $C_{rrrr}, C_{r\theta \phi r}$ 

## Steps (EXAMPLES/NED-Model)
*1.* Preparing parameter files (`DATA/Par_file_initmesh`, `DATA/meshfem3D_files/Mesh_Par_file`) and model files (e.g., interface files, tomographic files). Notes
- Make should your simulation region is bigger than the study region, which can be tuned in `step0_plot_region.sh`.
- If teleseismic simulation is required, please make sure the PML boundaries match the `WAVEFIELD_DISCONTINUITY_BOX_*` in `DATA/Par_file_initmesh`. Enable `COUPLE_WITH_INJECTION = .true.` and `INJECTINO_TYPE = 4`.
- In this step, the `MODEL` parameter should be `default`. If you want to set user defined model, use negative material number, and match the file 
with `nummaterial_velocity_file_tomo`.
- For realistic models and topographies, please refer to `create_model`.

*2.* Running the internal mesher (or an external mesher) to generate a cube-shaped mesh. In order to properly output the mesh files, the internal mesher needs to be run in serial (NPROCS=1). Then run `step1_run_mesher.sh`, set the target `NPROC`. 

*3.* Set the `lat lon rotate_angle` of the center of the study region (which can be printed from `step0_plot_region.sh`) in `step2_slurm_cube2sph.sh`, MPI is required in this step. Do not change anything in `Cube2sph_model_par` except the ellipticity.

*4.* Prepare `STATIONS`, `CMTSOLUTION` and `FORCESOLUTION`. All these files should be rotated to Cartesian coordinates by routines in `utils/cube2sph/bin/write*`. For teleseismic simulation, you need an empty `FORCESOLUTION`.

*5.* (optional) Prepare Injection wavefields, which can be handled by [AxiSEMLib](https://github.com/nqdu/AxiSEMLib/tree/main).

*6.* Launching the SPECFEM3D solver.

*7.* Performing rotation for the seismograms.