# MGCHDF5

These are fortran subroutines and MATLAB scripts to work with HDF5 MAGIC3D output.

**Fortran subroutines**:<br>

1. out3 (out3_mpi_h5.f90)<br>
Output full 3D domain of all q to hdf5 format file<br>
Output filename format: fort.qXXXX.h5<br>
<p>
  
2. outslicehor3 (out3_mpi_h5slicehor.f90)<br>
Output data only at altitudes of interest.<br>
Set next parameters in out3_mpi_h5slicehor.f90:<br>
numberofslices - set number of slices to output<br>
numslice - set cell number to output<br>
Output filename format: fort.qhXXXX.h5<br>
<p>
  
3. outslicever3 (out3_mpi_h5slicever.f90)<br>
Output data only at altitudes of interest.<br>
Set next parameters in out3_mpi_h5slicever.f90:<br>
arraysize<br>
idarray<br>
Both parameters can be found using "calcsliceh5.m" script<br>
Output filename format: fort.qvXXXX.h5<br>
<p>
  
4. restart (restart3_mpi_hdf.f90)<br>
Subroutine to restart from needed frame. Note that output used to restart should be output from out3 (out3_mpi_h5.f90) routine
<p>

In order to use these functions please add them in clawez_mpi.f
<p>
  
Compilation on Vega assumes changes in Makefile and modules used:<br>
1. Executable MAGIC should be compiled with next modules which should be loaded on Vega:
<p>
module load szip/intel-compiler/2.1.1<br>
module load hdf5/intel-mpi/intel-compiler/1.8.19<br>
module load intel/mpi/64/2018/1.163<br>
module load intel/compiler/64/2018/18.0.1<br>
Be sure to remove other compilers before compiling (module rm ......).
<p>
  
2. Makefile change:
Updated Makefile can be found in /MAGIC3D folder
