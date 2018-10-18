c     =======================================================
      subroutine viscous(maxmx,maxmy,maxmz,mbc,mx,my,mz,
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
c     Diffusion Solution for 3D Atmospheric Model Extensions
c        by J. B. Snively, 2003-2015
c
c  --------------------------------------------------------
c
c     Explicit First Order Diffusion Solver
c
c     # Three dimensions, one algorithm...
c
      dtdx2 = (dtr) / (dx*dx)
      dtdy2 = (dtr) / (dy*dy)
      dtdz2 = (dtr) / (dz*dz) 
c
      if (mdiffq .eq. 2) then
         dtdx2 = dtdx2 * (4/3)
         dt2dxdy = (dtr/3) / (2*dx*dy)
         dt2dxdz = (dtr/3) / (2*dx*dz)
      else if (mdiffq .eq. 3) then
         dtdy2 = dtdy2 * (4/3)
         dt2dxdy = (dtr/3) / (2*dx*dy)
         dt2dxdz = (dtr/3) / (2*dy*dz)
      else
         dtdz2 = dtdz2 * (4/3)
         dt2dxdy = (dtr/3) / (2*dx*dz)
         dt2dxdz = (dtr/3) / (2*dy*dz)
      endif 
c
c     # Apply diffusion
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
            if (mdiffq .eq. 2) then
              ad1a=q(i+1,j+1,k,3)
              ad2a=q(i-1,j+1,k,3)
              ad3a=q(i-1,j-1,k,3)
              ad4a=q(i+1,j-1,k,3)
              ad1b=q(i+1,j,k+1,4)
              ad2b=q(i-1,j,k+1,4)
              ad3b=q(i-1,j,k-1,4)
              ad4b=q(i+1,j,k-1,4)
            else if (mdiffq .eq. 3) then
              ad1a=q(i+1,j+1,k,2)
              ad2a=q(i-1,j+1,k,2)
              ad3a=q(i-1,j-1,k,2)
              ad4a=q(i+1,j-1,k,2)
              ad1b=q(i,j+1,k+1,4)
              ad2b=q(i,j-1,k+1,4)
              ad3b=q(i,j-1,k-1,4)
              ad4b=q(i,j+1,k-1,4)
            else
              ad1a=q(i+1,j,k+1,2)
              ad2a=q(i-1,j,k+1,2)
              ad3a=q(i-1,j,k-1,2)
              ad4a=q(i+1,j,k-1,2)
              ad1b=q(i,j+1,k+1,3)
              ad2b=q(i,j-1,k+1,3)
              ad3b=q(i,j-1,k-1,3)
              ad4b=q(i,j+1,k-1,3)
            endif
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
c
c         # First cross terms ('xy')
c
            aux(i,j,k,maux) = aux(i,j,k,maux) + 
     &         aux(i,j,k,mcoefaux) * dt2dxdy * 
     &         (ad1a-ad2a+ad3a-ad4a)
c
c         # Second cross terms ('xz')
c
            aux(i,j,k,maux) = aux(i,j,k,maux) +
     &         aux(i,j,k,mcoefaux) * dt2dxdz * 
     &         (ad1b-ad2b+ad3b-ad4b)
c
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
