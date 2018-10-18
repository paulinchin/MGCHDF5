c
c     =====================================================
      subroutine qinit(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,q,maux,aux)
c     =====================================================
c
c     # Set initial conditions for q.
c
      use ifport
      use mpi
      
	  implicit double precision (a-h,o-z)
	  

c
      external ohsteady, restart
c
      dimension q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &       1-mbc:maxmz+mbc,meqn)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &       1-mbc:maxmz+mbc,maux)
c
c     # MPI : Added so that we know which processor we are running on.
      common /mpi_proc_info/ np, id
      logical outt0
      common /restrt_block/ tinitial, iframe, outt0

c
      logical loadrestart
      logical usenoise
      data loadrestart /.false./
      data usenoise /.false./
      tweakfactor=1.d-5
c
      if (loadrestart) then
      tinitial = 200.d0
      iframe = 50
         call restart(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &                xlower,ylower,zlower,dx,dy,dz,q)     
         ohstep=1.d0
         goto 10
      else
         ohstep=0.d0
      endif
c
      if (usenoise) then    
         do i = 1-mbc,mx+mbc
            call srand(i*(id+1))
            do j = 1-mbc,my+mbc
               do k = 1-mbc,mz+mbc
                  tweak = tweakfactor*2*(rand()-0.5d0)
                  q(i,j,k,1) = aux(i,j,k,1)*(1+tweak) 
                  tweakmomentum = q(i,j,k,1)/aux(i,j,k,1)
                  q(i,j,k,2) = tweakmomentum*aux(i,j,k,2)
                  q(i,j,k,3) = tweakmomentum*aux(i,j,k,3)
                  q(i,j,k,4) = tweakmomentum*aux(i,j,k,4)
                  q(i,j,k,5) = aux(i,j,k,5)-
     &               0.5d0*(aux(i,j,k,2)**2+aux(i,j,k,3)**2+
     &               aux(i,j,k,4)**2)/aux(i,j,k,1)+
     &               0.5d0*(q(i,j,k,2)**2+q(i,j,k,3)**2+
     &               q(i,j,k,4)**2)/q(i,j,k,1)
                  do m = 6, meqn
                     q(i,j,k,m) = aux(i,j,k,m)
                  enddo 
               enddo
            enddo
         enddo
      else
         do i = 1-mbc,mx+mbc
            do j = 1-mbc,my+mbc
               do k = 1-mbc,mz+mbc
                  do m = 1,meqn
                     q(i,j,k,m) = aux(i,j,k,m)
                  enddo
               enddo
            enddo
         enddo
      endif
c
c     # Initialize Photochemistry
c
c 10   continue
 10   call ohsteady(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &     xlower,ylower,zlower,dx,dy,dz,q,maux,aux,ohstep)
c
c
c
      return
      end
