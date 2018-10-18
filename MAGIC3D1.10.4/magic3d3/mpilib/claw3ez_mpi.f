c
c
c
c     =================================================================
      subroutine claw3ez_mpi(mx,my,mz,meqn,mbc,mwaves,maux,method,
     &      mthlim,dtv,cflv,nv,nout,outstyle,tfinal,tout,
     &      nstepout,t0,xlower,xupper,ylower,yupper,zlower,
     &      zupper,mthbc,mwork)

c     =================================================================
c
c     An easy-to-use clawpack driver routine for simple applications
c     Documentation is available at
c                 http://www.amath.washington.edu/~claw/doc.html
c
c     Author: Randall J. LeVeque
c     Version of August, 2002 --  CLAWPACK Version 4.1
c
c
      implicit double precision (a-h,o-z)

      include 'mpif.h'

      external bc3,rpn3,rpt3,rptt3,src3,b4step3

c     # MPI : These arrays are automatically allocated here.  They will be
C     # automatically deallocated when we exit this routine.
      dimension    q(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc, meqn)
      dimension  aux(1-mbc:mx+mbc, 1-mbc:my+mbc,
     &               1-mbc:mz+mbc, maux)
      dimension work(mwork)

	integer popo, timesss
	character(4) myidd
	character(15) xtring
	real*8 inpp(1:1680000)

c     # MPI : These arrays allocated in claw3ez_mpi_driver, so they are only
C     # dimensioned here
      dimension mthlim(mwaves)
      dimension method(7),dtv(5),cflv(4),nv(2),mthbc(6)
      dimension tout(100)  !! dimensioned as tout(100) in claw3ez_mpi_driver.
      logical outt0
      common /restrt_block/ tinitial, iframe, outt0

c     # MPI : Added so that we know which processor we are running on.
      common /mpi_proc_info/ np, id
c
c
c     # -----------------------------------------------------------------
c     # MPI : We have removed all code that reads in data, and does error
c     # checking on input data.  We can start right in with other stuff..
c     # -----------------------------------------------------------------

c     # MPI : Set these here so that we are sure that they are defined
      maxmx = mx
      maxmy = my
      maxmz = mz

c
c     # grid spacing
      dx = (xupper - xlower) / float(mx)
      dy = (yupper - ylower) / float(my)
      dz = (zupper - zlower) / float(mz)
c
c	print*,dx
c	print*,dy
c	print*,dz

c     # time increments between outputting solution:
      if (outstyle .eq. 1) then
         dtout = (tfinal - t0)/float(nout)
      endif

c
c     # MPI : The file fort.nplot (unit=11) was created in claw3ez_mpi_driver,
C     # so we comment out the code here.

c       write(11,1101) nout
c       write(11,1101) 1
c c
c  1101 format(i5)

c
      if (id.eq.0) then
         open(10,file='fort.info',form='formatted')
      end if

c     # -------------------------------------------------------------------
c     # MPI : All code from this point is taken directly from claw3ez, and
C     # has not been changed.

c
c     # call user's routine setprob to set any specific parameters
c     # or other initialization required.
c
      do i = 0,np-1
         if (id.eq.i) then
            call setprob
c            write(*,*) 'Processor ', i, ' done with setprob'
         end if
         call mpi_barrier(MPI_COMM_WORLD,ierr)
      end do
c
c     # set aux array:
c
      if (maux .gt. 0)  then
         call setaux(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower,ylower,
     &               zlower,dx,dy,dz,maux,aux)
      endif
c
c     # set initial conditions:
      
      iframe = 0
      outt0  = .true.
      tinitial = t0
c
      call qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,xlower,ylower,
     &           zlower,dx,dy,dz,q,maux,aux)
c
c     # Reset initial time if changed within qinit (by restart3_mpi_hdf.f)
c
      t0 = tinitial
c      outt0 = .false.
      
      if (outt0) then
c        # output initial data

c	Do full output:
         call out3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,xlower,ylower,
     &          zlower,dx,dy,dz,q,t0,iframe)

c	Do vertical slice:
c	 call outslicever3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,xlower,
c     &          ylower,zlower,dx,dy,dz,q,t0,iframe)

c	Do horizontal slice:
c         call outslicehor3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,xlower,
c     &          ylower,zlower,dx,dy,dz,q,t0,iframe)

         if (id.eq.0) write(6,601) iframe, t0
      endif   
     
	popo=1000+id

	write(myidd,'(i4)') popo
	xtring = 'bcSP'//myidd
	xtring = xtring//'.dat'

        open(unit=id,file=xtring)
        do 490 j=1,1680000
        read(id,*) inpp(j)
490     continue

c       	if(id.eq.1) then
c       	print*,inpp(1)
c		endif
c        print*,'Boundary data file IS READ'
       close(id)

c      if (id==100) then
c      print*,id
c      print*,tend
c      print*,n0
c      print*,q
c      endif
       
c
c     ----------
c     Main loop:
c     ----------
c
      tend = t0
      n0   = iframe*nstepout + 1
  
      do 100 n=n0,nout
         tstart = tend
         if (outstyle .eq. 1)  tend = tstart + dtout
         if (outstyle .eq. 2)  tend = tout(n)
         if (outstyle .eq. 3)  tend = tstart - 1.d0  !# single-step mode
c
         call claw3(maxmx,maxmy,maxmz,meqn,mwaves,mbc,mx,my,mz,
     &           q,aux,xlower,ylower,zlower,dx,dy,dz,tstart,tend,dtv,
     &           cflv,nv,method,mthlim,mthbc,
     &           work,mwork,info,bc3,rpn3,rpt3,rptt3,src3,b4step3,inpp)
c
c        # check to see if an error occured:
c        if (info .ne. 0) then
c           write(6,*) '*** ERROR in claw3 ***  info =',info
c           go to 999
c           endif
c
         dtv(1) = dtv(5)  !# use final dt as starting value on next call
c
c        # output solution at this time
c        ------------------------------
c
c        # if outstyle=1 or 2, then nstepout=1 and we output every time
c        # we reach this point, since claw1 was called for the entire time
c        # increment between outputs.
c
c        # if outstyle=3 then we only output if we have taken nstepout
c        # time steps since the last output.

c        # iframe is the frame number used to form file names in out1
         iframe = n/nstepout
         if (iframe*nstepout .eq. n) then
            call out3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,xlower,ylower,
     &            zlower,dx,dy,dz,q,tend,iframe)

c         call outslice3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,xlower,
c     &          ylower,zlower,dx,dy,dz,q,t0,iframe)

c         call outslicehor3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,xlower,
c     &          ylower,zlower,dx,dy,dz,q,t0,iframe)

            if (id.eq.0) then
               write(6,601) iframe,tend
               write(10,1010) tend,info,dtv(3),dtv(4),dtv(5),
     &              cflv(3),cflv(4),nv(2)
            end if
         end if

c
c        # formats for writing out information about this call to claw:

  601    format('CLAW3EZ: Frame ',i4,
     &         ' matlab plot files done at time t = ', d12.4)
c
 1010    format('tend =',d15.4,/,
     &       'info =',i5,/,'smallest dt =',d15.4,/,'largest dt =',
     &       d15.4,/,'last dt =',d15.4,/,'largest cfl =',
     &         d15.4,/,'last cfl =',d15.4,/,'steps taken =',i4,/)
c
  100    continue
c
  999 continue
c
      if (id.eq.0) then
         close(10)
      end if

      return
      end
