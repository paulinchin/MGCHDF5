# MGCHDF5

These are fortran subroutines and MATLAB scripts to work with HDF5 MAGIC3D output.

**Fortran subroutines**:
1. out3 (out3_mpi_h5.f90)
Output full 3D domain of all q to hdf5 format file
Output filename format: fort.qXXXX.h5

2. outslicehor3 (out3_mpi_h5slicehor.f90)
Output data only at altitudes of interest.
Set next parameters in out3_mpi_h5slicehor.f90:
numberofslices - set number of slices to output
numslice - set cell number to output
Output filename format: fort.qhXXXX.h5

3. outslicever3 (out3_mpi_h5slicever.f90)
Output data only at altitudes of interest.
Set next parameters in out3_mpi_h5slicever.f90:
arraysize
idarray
Both parameters can be found using "calcsliceh5.m" script
Output filename format: fort.qvXXXX.h5

4. restart (restart3_mpi_hdf.f90)
Subroutine to restart from needed frame. Note that output used to restart should be output from out3 (out3_mpi_h5.f90) routine

Compilation on Vega assumes changes in Makefile and modules used:
1. Executable MAGIC should be compiled with next modules which should be loaded on Vega:

module load szip/intel-compiler/2.1.1<br>
module load hdf5/intel-mpi/intel-compiler/1.8.19<br>
module load intel/mpi/64/2018/1.163<br>
module load intel/compiler/64/2018/18.0.1<br>

Be sure to remove other compilers before compiling (module rm ......).
2. Makefile change:
Updated Makefile can be found in /MAGIC3D folder
