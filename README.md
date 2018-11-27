# MGCHDF5

Fortran subroutines and MATLAB scripts to work with HDF5 MAGIC3D output. This README file is for parallel compressed version from MAGIC3Db1.10.4 version.

<b>Current versions</b>:
1. MAGIC3Dc1.10.4 - parallel output/restart routines with compression working with HDF5 development branch for MAGIC3D2018b version which is compatible with CLAWPACK 5.5.0<br>
2. MAGIC3Db1.10.4 - parallel output/restart routines with compression working with HDF5 development branch for MAGIC3D2018b version which is compatible with CLAWPACK 5.5.0 (depreciated)<br>
3. MAGIC3D1.10.4 - parallel output/restart routines with compression working with HDF5 development branch which is compatible with CLAWPACK < 5.5.0 (depreciated)<br>
3. MAGIC2D - parallel output/restart routines for 2D version of MAGIC with compression working with HDF5 development branch which is compatible with CLAWPACK < 5.5.0 <br>

<b>Important note</b>:<br>
Even though MAGIC3D2018b version uses new indexing (m,i,j,k), the output routines still saves data as (i,j,k,m) in order to be compatible with Matlab routines.

## Fortran subroutines<br>

1. out3 (out3_mpi_h5.f90)<br>
Output full 3D domain of all q to hdf5 format file<br>
Output filename format: fort.qXXXX.h5<br>
<p>
  
2. out3h (out3_mpi_h5slicehor.f90)<br>
Output data only at altitudes of interest.<br>
Set next parameters in out3_mpi_h5slicehor.f90:<br>
numberofslices - set number of slices to output<br>
numslice - set cell number to output<br>
Output filename format: fort.qhXXXX.h5<br>
<p>
  
3. out3a (out3_mpi_h5sliceair.f90)<br>
Output data for the continuous range of altitudes (e.g. airglow output).<br>
Set next parameters in out3_mpi_h5slicehorr.f90:<br>
dimsf(3) - set number of slices to output<br>
k - set range of grids in "z" direction to output<br>
Output filename format: fort.qaXXXX.h5<br>
<p>
  
4. out3v (out3_mpi_h5slicever.f90)<br>
Output data only at altitudes of interest.<br>
Set next parameters in out3_mpi_h5slicever.f90:<br>
arraysize<br>
idarray<br>
Both parameters can be found using "calcsliceh5.m" script<br>
Output filename format: fort.qvXXXX.h5<br>
<p>
  
5. restart (restart3_mpi_hdf.f90)<br>
Subroutine to restart from needed frame. Note that output used to restart should be output from out3 (out3_mpi_h5.f90) routine
<p>

<b>Notes</b>:<br>
1. In order to use these functions please add them in clawez_mpi.f<br>
2. If several different outputs are needed (for example full output for every 1 min and slice for every 1 sec) use different conditions in claw3ez_mpi.f for *if (iframe*nstepout .eq. n) then* (updated code with setting flags will be available soon).

<p>
  
## Compilation
Compilation on Vega assumes changes in Makefile and used modules:<br>
1. Executable MAGIC should be compiled with next modules:
<p>
module load hdf5/intel-mpi/intel-compiler/1.10.4<br>
module load intel/mpi/64/2018/1.163<br>
module load intel/compiler/64/2018/18.0.1<br>
Be sure to remove other compilers before compiling (module rm ......).
<p>
  
2. Makefile change:
Updated Makefile can be found in /MAGIC3D1.10.4 folder

<b>NOTE</b><br>
Be sure to compile new routines before running Makefile using:<br>
h5pfc -O3 -xhost -ipo -pthread -c out3_mpi_h5.f90<br>
h5pfc -O3 -xhost -ipo -pthread -c out3_mpi_h5slicehor.f90<br>
h5pfc -O3 -xhost -ipo -pthread -c out3_mpi_h5slicever.f90<br>
h5pfc -O3 -xhost -ipo -pthread -c out3_mpi_h5sliceair.f90<br>
h5pfc -O3 -xhost -ipo -pthread -c restart3_mpi_hdf.f90<br>
Because they are written for f90 and are compiled without all other .f routines

## Matlab script<br>
<b>NOTE</b><br>
There are 2 sets of Matlab scripts: MAGIC3Dc1.10.4 works with MATLABc1.10.4. MAGIC3Db1.10.4 and MAGIC3D1.10.4 work with MATLAB.

Matlab scripts can import the whole 3D domain (mx,my,mz,q) as it was done for earlier versions, however it may take a lot of time and resources, thus this function is depreciated. New functionality allows importing and working with only particular horizontal (in x or y) and horizontal slices.
<p>
  <b>Current scripts </b><br>
  1. calcsliceh5.m - calculates which threads should be output when out3ver output subroutine is used<br>
  2. readmagich5.m - main script to set parameters and run output routines<br>
  3. initialization.m - load MSIS profile and reads simulation parameters<br>
  4. slicing.m - contains function to retrieving data from h5 binary output files. Current version allows:<br>
    a. Retrieve vertical slice in x or y direction setting slice position in meters<br>
    b. Retrieve horizontal slice in z direction setting slice position in meters<br>
    c. Load full 3D domain of all variables<br>
  The script allows working with all 3 types of outputs (out3, out3ver, out3hor)<br>
  The script allows retrieving only particular variables (in current version: 'u','v','w','rhop','rhorp','doxp','dnit2p','dox2p')<br>
  The script allows make plots and video output<br>
  5. plotting.m - this routine calculates needed variables, does or does not make plots<br>
  6. afterframeh5.m - to calculate and plot airglow (based on afterframe.m)
