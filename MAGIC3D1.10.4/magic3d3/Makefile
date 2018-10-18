#
#               Makefile for the mpiclaw code:
#
#       To make an MPI executable, type:          make xclawmpi
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xclawmpihdf
#       (from the application sub-directory)
#
#       To combine ASCII output files (one from
#       each processor at each time level) into
#       fort.qXXXX files for use with MATLAB, type:  make catfiles
#       (in directory with fort.qXXXX.YY files)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
FFLAGS = -O3 -xhost -ipo -pthread
LFLAGS = $(FFLAGS)
F77    = h5pfc -c
LINK   = h5pfc

# HDFLIBS = -L/usr/lib64 -lz -ljpeg \
#          -L/scratch/inchinp/mylib/lib \
#          -lsz -lmfhdf -ldf

# HDFLIBS = h5pfc

HDF5 = /cm/shared/apps/hdf5/intel-mpi/intel-compiler/1.8.19
FORTRANLIB=-I$(HDF5)/include $(HDF5)/lib/libhdf5_fortran.a

HDFLIBS = $(FORTRANLIB) $(HDF5)/lib/libhdf5.a -lz -lm -static-intel
.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

OBJECTS = \
  b4step3.o \
  bc3_mpi.o \
  qinit.o \
  rpn3euX.o \
  rpt3euX.o \
  rptt3euX.o \
  fcns.o \
  setaux.o \
  setprob.o \
  viscous.o \
  conduct.o \
  ohsteady.o \
  src3.o

LIBOBJECTS = \
  mpilib/claw3ez_mpi_driver.o \
  mpilib/claw3ez_mpi.o \
  mpilib/claw3_mpi.o \
  mpilib/bc3_aux_mpi.o \
  lib/chkmth.o \
  lib/step3.o \
  lib/step3ds.o \
  lib/dimsp3.o \
  lib/flux3fw.o \
  lib/copyq3.o \
  lib/inlinelimiter.o

OUTOBJECTS = mpilib/out3_mpi.o
# OUTOBJECTS = out3_mpi_h5.o

HDFOBJECTS = out3_mpi_h5.o \
  restart3_mpi_hdf.o \
  out3_mpi_h5slicever.o \
  out3_mpi_h5slicehor.o

# mpilib/out3_mpi_hdf.o

xclawmpi: $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclawmpi 

# xclawmpihdf: $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS)
#         $(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) -o xclawmpi
# 	$(HDFLIBS)

xclawmpihdf: $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) -o xclawmpi


# Type "make catfiles" to combine fort.qXXXX.YY output files from out3_mpi
# into fort.qXXXX files (for use with CLAWPACK's MATLAB graphics routines).
# This step is not necessary with the HDF output routines.
catfiles : fort.q0000

fort.q0000 : fort.q0000.00
	mpilib/catfiles

fort.q0000.00 :
	mpilib/catfiles

### DO NOT remove this line - make depends on it ###
