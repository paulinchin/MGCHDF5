c     =======================================================
      subroutine conduct(maxmx,maxmy,maxmz,mbc,mx,my,mz,
     & xlower,ylower,zlower,dx,dy,dz,
     & meqn,q,mdiffq,maux,aux,mcoefaux,dtr)
c     =======================================================
c
      implicit double precision (a-h,o-z)
c
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)
c
c  --------------------------------------------------------
c
c     Conduction Solution for 3D Atmospheric Model Extensions
c        by J. B. Snively, 2003-2015
c
c  --------------------------------------------------------
c
c     Explicit First Order Conduction Solver
c
      dtdx2 = (dtr) / (dx*dx)
      dtdy2 = (dtr) / (dy*dy)
      dtdz2 = (dtr) / (dz*dz) 
c 
c     # Apply conduction
c
      do k=1,mz
        do j=1,my
          do i=1,mx
c         # Quantity (minus, center, plus, up, down)
            am=q(i-1,j,k,mdiffq)
            ac=q(i,j,k,mdiffq)
            ap=q(i+1,j,k,mdiffq)
            ab=q(i,j-1,k,mdiffq)
            af=q(i,j+1,k,mdiffq)
            au=q(i,j,k+1,mdiffq)
            ad=q(i,j,k-1,mdiffq)
c
c         # Streamwise direction (x)
c
            aux(i,j,k,maux) = aux(i,j,k,mcoefaux) *
     &         dtdx2*(ap - 2.0d0*ac + am)     
c
c         # Transverse direction (y)
c
            aux(i,j,k,maux) = aux(i,j,k,maux) +
     &         aux(i,j,k,mcoefaux) *
     &         dtdy2*(af - 2.0d0*ac + ab)
c
c         # Vertical direction (z)
c
            aux(i,j,k,maux) = aux(i,j,k,maux) +
     &         aux(i,j,k,mcoefaux) *
     &         dtdz2*(au - 2.0d0*ac + ad)
          enddo
        enddo
      enddo
c
c     # Reassign q+aux(maux) to q
c
      do k=1,mz
        do j=1,my
          do i=1,mx
            q(i,j,k,mdiffq)=q(i,j,k,mdiffq)+aux(i,j,k,maux)
          enddo
        enddo
      enddo
c
c
c
      return
      end
