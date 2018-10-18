c
c     =====================================================
      subroutine restart(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,q)
c     =====================================================
c
c
c     # Set initial conditions for q.
c
      implicit double precision (a-h,o-z)
      include 'mpif.h'
c
c
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &      1-mbc:maxmz+mbc, meqn)
      character*10 fname1, fname2
      logical outt0
      common /restrt_block/ tinitial, iframe, outt0

c
c     # MPI: get number of processors and id of this process.
c
      call mpi_comm_rank(mpi_comm_world,id,ierr)

      if (id.eq.0) then
         write(*,*) '*** Restart is not supported from fort.qXXXX ',
     &        'output files from MPICLAW'
         write(*,*) '*** The HDF output version of MPICLAW supports ',
     &        'restarting'
      end if
      call mpi_finalize(ierr)
      STOP 

      return
      end
