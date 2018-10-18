c
c
c =========================================================
      subroutine out3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,
     &     mz,xlower,ylower,zlower,dx,dy,dz,q,t,iframe)
c =========================================================
c
c     # Output the results for a general system of conservation laws
c     # in 3 dimensions
c
      implicit double precision (a-h,o-z)
      include 'mpif.h'
	real*8 pmod
	integer altits
	integer ff
	dimension altits(5)
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &     1-mbc:maxmz+mbc,meqn)
      dimension mtotal(3), mstart(3)
      character*100 fname1, fname2

c
      common /mpicomm/ mpi_comm_3d, lx, ly, lz, mtotal, mstart 
      common /mpi_proc_info/ np, id
c
c     # Write the results to the file fort.q<iframe>.<id>
c     # The files fort.q<iframe>.xxx can be concatenated together
c     # into a single file fort.q<iframe> using the script
c     # 'catfiles'.  The single fort.q<iframe> will be in a format
c     # suitable for the matlab script plotclaw2.m

c     # The same format is used by the amrclaw package.
c     # Here it's adapted to output just the single grid.
c
c     # first create the file name and open file.
c
	pmod = MOD(id-5,10)

	altits(1)=1
	altits(2)=51
	altits(3)=101
	altits(4)=151
	altits(5)=199

      print*,dx
	print*,dy
	print*,dz

c	if((pmod.eq.0).or.(id.eq.0)) then

        if((id.eq.0).or.(id.eq.40).or.(id.eq.120).or.(id.eq.200)) then

c	if((t.eq.0).or.(t.eq.520).or.(t.eq.640).or.(t.eq.760)
c     +  .or.(t.eq.880).or.(t.eq.1000).or.
c     +  (t.eq.1200).or.(t.eq.1500)) then



c        if((t.eq.400).or.(t.eq.520).or.(t.eq.640).or.(t.eq.760)) then

	fname1 = 'fort.qxxxx.xxx'
      fname2 = 'fort.txxxx'
      nstp = iframe
      do 55 ipos = 10, 7, -1
         idigit = mod(nstp,10)
         fname1(ipos:ipos) = char(ichar('0') + idigit)
         fname2(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
 55   continue

c     # Now add in processor id to file name.  We will cat the files later
c     # so we can read them into matlab.

      if (id == 0 .and. np > 999) then
               write(6,*) 'out2mpi : np > 999'
         stop
      endif

      nstp = id
      do ipos = 14,12,-1
         idigit = mod(nstp,10)
         fname1(ipos:ipos) = char(ichar('0') + idigit)
         nstp = nstp / 10
      enddo

      iunitq = 60 + id
      open(unit=iunitq,file=fname1,status='unknown',form='formatted')

c     # Only node 0 needs to write out the time file.
      if (id .eq. 0) then
         iunitt = 50
         open(unit=iunitt,file=fname2,status='unknown',form='formatted')
      endif

c
c     # the following parameters are used in amrclaw where there are
c     # multiple grids.   Here, all grids are assumed to be at level 1.
c     # but we may have several grids.
      ngrids = np
      mptr = id + 1
      level = 1

      write(iunitq,1001) mptr,level,mx,my,mz
 1001 format(i5,'                 grid_number',/,
     &       i5,'                 AMR_level',/,
     &       i5,'                 mx',/,
     &       i5,'                 my',/,
     &       i5,'                 mz')

      write(iunitq,1002) xlower,ylower,zlower,dx,dy,dz
 1002 format(e18.8,'    xlow', /,
     &       e18.8,'    ylow', /,
     &       e18.8,'    zlow', /,
     &       e18.8,'    dx', /,
     &       e18.8,'    dy', /,
     &       e18.8,'    dz',/)
c


       do 30 k = 1,mz
c	k = 175
c       do 30 ff =1,5
c	k=altits(ff)
         do 20 j=1,my
            do 10 i=1,mx
               do m=1,meqn
c     # exponents with more than 2 digits cause problems reading
c     # into matlab... reset tiny values to zero:
       if (dabs(q(i,j,k,m)) .lt. 1d-99) q(i,j,k,m) = 0.d0
               enddo
c
               write(iunitq,1005) (q(i,j,k,m), m=1,meqn)
 1005          format(4e24.16)
c
 10         continue
            write(iunitq,*) ' '
 20      continue
         write(iunitq,*) ' '
 30   continue
      write(iunitq,*) ' '

c     # Write out file fort.t<iframe>  We only need a single file - let node 0
c     # do the job.


      if (id == 0) then
c         print*,maxmx
c	 print*,maxmy
c	 print*,maxmz
c	 print*,meqn
c	 print*,mbc
c	 print*,mx
c	 print*,my
c	 print*,mz
c	 print*,xlower
c	 print*,ylower
c	 print*,zlower
c	 print*,dx
c	 print*,dy
c	 print*,dz
c	 print*,t
c	 print*,iframe
	
         write(iunitt,1000) t,meqn,ngrids
 1000    format(e18.8,'    time', /,
     &        i5,'                 meqn'/,
     &        i5,'                 ngrids'/,/)
      endif
c
      close(unit=iunitq)


      if (id == 0) then
         close(unit=iunitt)
      endif

	endif

      return
      end
