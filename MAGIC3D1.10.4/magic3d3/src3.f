c
c     =======================================================
      subroutine src3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &                xlower,ylower,zlower,dx,dy,dz,q,
     &                maux,aux,t,dt,mthbc)
c     =======================================================
c
c
      implicit double precision (a-h,o-z)
      external gaussian, forcefnz, viscous, 
     &         conduct, bc3, ohsteady
c
      common /param/ gamma, gamma1, grav, dens, pres, 
     &                  speed, dtherm0, dmolec0
      common /forcing/ omega, amplitude, propx, propz, 
     &                  w_n, w_o, vsrcx, forcemth
      common /wind/ omega1, amp0, amp1, prop1z, z1width, z1pos,
     &                  z1Gpos, phi, windscale1
      common /gauss/ xpos, ypos, zpos, xwidth, ywidth,
     &               zwidth, tcenter, twidth
      common /options/ loadprof, nspecadv
c
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, meqn) 
      dimension    aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux) 
c
      dimension mthbc(6)
c
c  --------------------------------------------------------
c
c     Source Terms for Atmospheric Model Extensions
c        by J. B. Snively, 2003-2018
c
c  --------------------------------------------------------
c
c
c     # Apply viscous diffusion and thermal conduction
c 
  100 do i=1-mbc,mx+mbc
         do j=1-mbc,my+mbc
            do k=1-mbc,mz+mbc
            gamma=aux(i,j,k,19)
            gamma1=gamma-1
c        # Convert energy to temperature
            q(i,j,k,5) = gamma1* 
     &         (q(i,j,k,5)-0.5d0*(1/q(i,j,k,1))*(q(i,j,k,2)**2
     &        + q(i,j,k,3)**2 + q(i,j,k,4)**2))
     &        /(q(i,j,k,1)*aux(i,j,k,18))
            auxtemp = gamma1*
     &         (aux(i,j,k,5)-0.5d0*(1/aux(i,j,k,1))
     &        *(aux(i,j,k,2)**2 + aux(i,j,k,3)**2 
     &        + aux(i,j,k,4)**2))
     &        /(aux(i,j,k,1)*aux(i,j,k,18))
            q(i,j,k,5) = q(i,j,k,5)-auxtemp
c        # Subtract Wind from u, convert momentum to velocity
            q(i,j,k,2) = q(i,j,k,2)/q(i,j,k,1)
     &        -aux(i,j,k,2)/aux(i,j,k,1)
c        # Convert momentum to velocity
            q(i,j,k,3) = q(i,j,k,3)/q(i,j,k,1)
     &        -aux(i,j,k,3)/aux(i,j,k,1)
c        # Convert momentum to velocity
            q(i,j,k,4) = q(i,j,k,4)/q(i,j,k,1)
     &        -aux(i,j,k,4)/aux(i,j,k,1)
            enddo
         enddo
      enddo
c
c     # Run diffusion and conduction algorithms
c
c     # Calculate fractional time step
      dxyzmin=min(dx,dy,dz)
      difmax=max(aux(1,1,mz,11),aux(1,1,mz,12))
      difcfl=0.1d0
      ndtr=max(int((dt*difmax/(dxyzmin*dxyzmin))/difcfl),1)
      dtr=dt/ndtr
c
c     # Substep...
c
      do n=1,ndtr
c
c     # Calculate intermediate time
         tr=t+n*dtr      
c
c     # X-Velocity Viscosity
         call viscous(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower,
     &      ylower,zlower,dx,dy,dz,meqn,q,2,maux,aux,11,dtr)
c     # Y-velocity Viscosity
         call viscous(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower, 
     &      ylower,zlower,dx,dy,dz,meqn,q,3,maux,aux,11,dtr)
c     # Z-velocity Viscosity
         call viscous(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower, 
     &      ylower,zlower,dx,dy,dz,meqn,q,4,maux,aux,11,dtr)
c     # Thermal Conduction
         call conduct(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower,
     &      ylower,zlower,dx,dy,dz,meqn,q,5,maux,aux,12,dtr)
c
c    # Apply boundary conditions
c
c     Disable lower boundary when calling viscous substeps
c
      mthbc_save_5=mthbc(5)
      mthbc(5)=-1
c
      call bc3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &         xlower,ylower,zlower,dx,dy,dz,
     &         q,maux,aux,tr,dtr,mthbc)
c
      mthbc(5)=mthbc_save_5
c
      enddo
c
c    # Reassemble q variables
c
      do i=1-mbc,mx+mbc
         xcell = xlower + (i-0.5d0)*dx
         do j=1-mbc,my+mbc
            ycell = ylower + (j-0.5d0)*dy            
            do k=1-mbc,mz+mbc
            zcell = zlower + (k-0.5d0)*dz            
            gamma = aux(i,j,k,19)
            gamma1= gamma-1
c        # Update time-dependent wind field
            if (omega1 .NE. 0) then
               aux(i,j,k,5) = aux(i,j,k,5)-
     &            0.5d0*((aux(i,j,k,2))**2)/aux(i,j,k,1)
               scale = (aux(i,j,1,1)/aux(i,j,k,1))**0.5d0
               aux(i,j,k,2) = aux(i,j,k,2)-aux(i,j,k,1)
     &            *windfnx(xcell,ycell,zcell,t-dt,scale)
               aux(i,j,k,2) = aux(i,j,k,2)+aux(i,j,k,1)
     &            *windfnx(xcell,ycell,zcell,t,scale)
               aux(i,j,k,5) = aux(i,j,k,5)+     
     &            0.5d0*((aux(i,j,k,2))**2)/aux(i,j,k,1)
            endif
c        # Convert velocity to momentum
            q(i,j,k,2) = q(i,j,k,2)*q(i,j,k,1) +
     &              (aux(i,j,k,2)/aux(i,j,k,1))*q(i,j,k,1)
            q(i,j,k,3) = q(i,j,k,3)*q(i,j,k,1) +
     &              (aux(i,j,k,3)/aux(i,j,k,1))*q(i,j,k,1)
            q(i,j,k,4) = q(i,j,k,4)*q(i,j,k,1) +
     &              (aux(i,j,k,4)/aux(i,j,k,1))*q(i,j,k,1)
c        # Convert temperature perturbation to energy density
            auxtemp = gamma1*
     &         (aux(i,j,k,5)-0.5d0*(1/aux(i,j,k,1))*
     &         (aux(i,j,k,2)**2 + aux(i,j,k,3)**2 + 
     &          aux(i,j,k,4)**2))/(aux(i,j,k,1)*aux(i,j,k,18))
            q(i,j,k,5) = (auxtemp + q(i,j,k,5))* 
     &          q(i,j,k,1)*aux(i,j,k,18)/gamma1 +
     &          0.5d0*(1/q(i,j,k,1))*(q(i,j,k,2)**2 + 
     &          q(i,j,k,3)**2 + q(i,j,k,4)**2)
            enddo
         enddo
      enddo
c
c
c     # OH Airglow Chemistry
c
      call ohsteady(maxmx,maxmy,maxmz,meqn,mbc,mx,my,my,
     &     xlower,ylower,zlower,dx,dy,dz,q,maux,aux,dt)
c
c
c     # Oscillatory Forcing Function
c 
  200 if (forcemth.EQ.0) goto 300
  250 do i=1,mx
        xcell = xlower + (i-0.5d0)*dx - t*vsrcx
        do j=1,my
          ycell = ylower + (j-0.5d0)*dy
          do k=1,mz
            zcell = zlower + (k-0.5d0)*dz  
c       # Vertical Forcing
            bf=dt*amplitude*
     &            gaussian(xcell,ycell,zcell,t)*
     &            forcefnz(xcell,ycell,zcell,t)
            q(i,j,k,4) = q(i,j,k,4)+bf*q(i,j,k,1)
          enddo
        enddo
      enddo  
c
c
c     # Correct Negative Densities
c     # (This shouldn't happen, but with layers it might.)
c
      do m = 6, meqn
        do k = 1-mbc, mz+mbc
          do j = 1-mbc, my+mbc
            do i = 1-mbc, mx+mbc
              if (q(i,j,k,m).LT.1.d-14) q(i,j,k,m)=1.d-14
            enddo
          enddo
        enddo
      enddo
c
c      
  300 return
      end
