#
#               Makefile for the clawpack code:
#
#       To make an executable, type:            make xclawmpi
#       (from the application sub-directory)
#
#       To make an executable that generates
#       output in HDF (version 4) format, type:   make xclawmpihdf
#       (from the application sub-directory)
#
#       To compile a single file.f type:          make file.o
#       (from the application sub-directory)
#
#       To compile the lib library files type:    make lib
#       (from the claw/clawpack/2d/lib sub-directory)
#
#       To combine all source files type:         make program
#       (from the application sub-directory)
#       (sometimes easier to debug all in one)
#
FFLAGS = -O3 -ipo -xhost
LFLAGS = $(FFLAGS)
F77    = h5pfc -c
# mpif77 -fc=ifort -c -mcmodel medium -shared-intel
LINK   = h5pfc
# mpif77 -fc=ifort
HDFLIBS = -L/usr/lib64 -lz \
          -L/scratch/pineyrob/mylib/lib -ljpeg\
          -lsz -lmfhdf -ldf

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

OBJECTS = \
  qinit.o \
  rpn2euX-ThirdOrder.o \
  rpt2euX.o \
  setprob.o \
  b4step2.o \
  bc2_mpi_tracer.o \
  bc2_diffuse.o \
  src2.o \
  setaux.o \
  viscous1.o \
  viscous2.o \
  conduct.o \
  ohdata.o \
  ohsteady.o \
  fcns.o \
  out2_mpi_h5.o
  
LIBOBJECTS = \
  mpilib/claw2ez_mpi_driver.o \
  mpilib/claw2ez_mpi.o \
  mpilib/claw2_mpi.o \
  mpilib/bc2_aux_mpi.o \
  lib/step2.o \
  lib/step2ds.o \
  lib/dimsp2.o \
  lib/flux2fw.o \
  lib/copyq2.o \
  lib/inlinelimiter.o

SOURCES = \
  qinit.f \
  rpn2euX-ThirdOrder.f \
  rpt2euX.f \
  setprob.f \
  bc2_mpi_tracer.f \
  bc2_diffuse.f \
  b4step2.f \
  src2.f \
  setaux.f \
  viscous1.f \
  viscous2.f \
  conduct.f \
  ohdata.f \
  ohsteady.f \
  fcns.f \
  out2_mpi_h5.f90 \
  mpilib/claw2ez_mpi_driver.f \
  mpilib/claw2ez_mpi.f \
  mpilib/claw2_mpi.f \
  mpilib/bc2_aux_mpi.f \
  lib/b4step2.f \
  lib/step2.f \
  lib/step2ds.f \
  lib/dimsp2.f \
  lib/flux2fw.f \
  lib/copyq2.f \
  lib/inlinelimiter.f

# OUTOBJECTS = mpilib/out2_mpi.o \
#         mpilib/restart2_mpi.o
# OUTSOURCES = mpilib/out2_mpi.f \
#        mpilib/restart2_mpi.f
  
HDFOBJECTS = out2_mpi_h5.o
# out2_mpi_hdf.o
#        mpilib/restart2_mpi_hdf.o
HDFSOURCES = out2_mpi_h5.f90
# mpilib/out2_mpi_hdf.f \
#        mpilib/restart2_mpi_hdf.f
  
xclawmpi: $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(OUTOBJECTS) -o xclawmpi
  
#xclawmpihdf: $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS)
#	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS) -o xclawmpi \
#	$(HDFLIBS)

xclawmpihdf: $(OBJECTS) $(LIBOBJECTS) $(HDFOBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) -o xclawmpi
  
program : $(SOURCES) $(OUTSOURCES)   
	cat $(SOURCES) $(OUTSOURCES) > claw2program.f
  
programhdf : $(SOURCES) $(HDFSOURCES)
	cat $(SOURCES) $(HDFSOURCES) > claw2program.f
  
# Type "make catfiles" to combine fort.qXXXX.YY output files from out2_mpi
# into fort.qXXXX files (for use with CLAWPACK's MATLAB graphics routines).
# This step is not necessary with the HDF output routines.
catfiles : fort.q0000

fort.q0000 : fort.q0000.00
	mpilib/catfiles

fort.q0000.00 :
	mpilib/catfiles

### DO NOT remove this line - make depends on it ###
