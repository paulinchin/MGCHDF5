#
#		Makefile for the clawpack code:
#
#	To compile the claw3 library files type:       make lib
#       (from the lib sub-directory)
#
#	To make an executable, type:                   make xclaw
#       (from the application sub-directory)
#
#	To compile a single file.f type:               make file.o
#       (from the application sub-directory)
#
#	To combine source all source:                  make backup
#       (date is appended)
#
#
#
FFLAGS = -O
LFLAGS =
F77 = f77 -c
F90 = f95

.f.o: ; $(F77) $(FFLAGS) $*.f -o $*.o

LIBOBJECTS = \
  $(CLAW)/clawpack/3d/lib/claw3ez.o \
  $(CLAW)/clawpack/3d/lib/bc3.o \
  $(CLAW)/clawpack/3d/lib/setprob.o \
  $(CLAW)/clawpack/3d/lib/setaux.o \
  $(CLAW)/clawpack/3d/lib/b4step3.o \
  $(CLAW)/clawpack/3d/lib/claw3.o \
  $(CLAW)/clawpack/3d/lib/step3.o \
  $(CLAW)/clawpack/3d/lib/step3ds.o \
  $(CLAW)/clawpack/3d/lib/dimsp3.o \
  $(CLAW)/clawpack/3d/lib/flux3.o \
  $(CLAW)/clawpack/3d/lib/copyq3.o \
  $(CLAW)/clawpack/3d/lib/chkmth.o \
  $(CLAW)/clawpack/3d/lib/inlinelimiter.o \
  $(CLAW)/clawpack/3d/lib/limiter.o \
  $(CLAW)/clawpack/3d/lib/philim.o \
  $(CLAW)/clawpack/3d/lib/src3.o \
  $(CLAW)/clawpack/3d/lib/rpn3.o \
  $(CLAW)/clawpack/3d/lib/rpt3.o \
  $(CLAW)/clawpack/3d/lib/rptt3.o \
  $(CLAW)/clawpack/3d/lib/cellave.o \
  $(CLAW)/clawpack/3d/lib/fdisc.o \
  $(CLAW)/clawpack/3d/lib/out3.o \
  $(CLAW)/clawpack/3d/lib/out3_hdf.o \
  $(CLAW)/clawpack/3d/lib/restart3.o \
  $(CLAW)/clawpack/3d/lib/restart3_hdf.o

SOURCES = \
  $(CLAW)/clawpack/3d/lib/claw3ez.f \
  $(CLAW)/clawpack/3d/lib/bc3.f \
  $(CLAW)/clawpack/3d/lib/setprob.f \
  $(CLAW)/clawpack/3d/lib/setaux.f \
  $(CLAW)/clawpack/3d/lib/b4step3.f \
  $(CLAW)/clawpack/3d/lib/claw3.f \
  $(CLAW)/clawpack/3d/lib/step3.f \
  $(CLAW)/clawpack/3d/lib/step3ds.f \
  $(CLAW)/clawpack/3d/lib/dimsp3.f \
  $(CLAW)/clawpack/3d/lib/flux3.f \
  $(CLAW)/clawpack/3d/lib/copyq3.f \
  $(CLAW)/clawpack/3d/lib/chkmth.f \
  $(CLAW)/clawpack/3d/lib/inlinelimiter.f \
  $(CLAW)/clawpack/3d/lib/limiter.f \
  $(CLAW)/clawpack/3d/lib/philim.f \
  $(CLAW)/clawpack/3d/lib/src3.f \
  $(CLAW)/clawpack/3d/lib/rpn3.f \
  $(CLAW)/clawpack/3d/lib/rpt3.f \
  $(CLAW)/clawpack/3d/lib/rptt3.f \
  $(CLAW)/clawpack/3d/lib/cellave.f \
  $(CLAW)/clawpack/3d/lib/fdisc.f \
  $(CLAW)/clawpack/3d/lib/out3.f \
  $(CLAW)/clawpack/3d/lib/out3_hdf.f \
  $(CLAW)/clawpack/3d/lib/restart3.f \
  $(CLAW)/clawpack/3d/lib/restart3_hdf.f

lib: $(LIBOBJECTS)

hdf: out3_hdf.o outhdf3.o

out3_hdf.o: out3_hdf.f
	$(F77) -c -I/usr/local/hdf/include out3_hdf.f

outhdf3.o: outhdf3.f90
	$(F90) -c -I/usr/local/hdf/include outhdf3.f90

xclaw: $(OBJECTS)
	$(LINK) $(LFLAGS) $(OBJECTS) $(LIBOBJECTS) -o xclaw

backup : $(SOURCES)
	cat $(SOURCES) $(IFS) > clawpack3lib.`date +%y%m%d`

clean:
	-rm -f $(LIBOBJECTS)

### DO NOT remove this line - make depends on it ###
