c
c
c     ==================================================================
      subroutine setaux(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower,ylower,
     &                  zlower,dx,dy,dz,maux,aux)
c     ==================================================================
c
c     # set auxiliary arrays 
c
c     
      implicit double precision (a-h,o-z)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 
     &               1-mbc:maxmz+mbc, maux)

      common /param/ gamma, gamma1, grav, dens, pres,
     &                  speed, dtherm0, dmolec0
      common /wind/ omega1, amp0, amp1, prop1z, z1width, z1pos,
     &                  z1Gpos, phi, windscale1
      common /gauss/ xpos, ypos, zpos, xwidth, ywidth,
     &               zwidth, tcenter, twidth
      common /options/ loadprof, nspecadv
c
c     # MPI : Added so that we know which processor we are running on.
      common /mpi_proc_info/ np, id
c
c
c  --------------------------------------------------------
c
c     Aux. Initialization for Atmospheric Model Extensions
c        by J. B. Snively, 2003-2013
c
c  --------------------------------------------------------
c
c     # Open profile if desired
c
      if (loadprof.EQ.1) then
         open(unit=7,file='profile.data',
     &            status='old',form='formatted')
      endif
c
c     # Initialize auxiliary variables
      do 15 i=1-mbc,mx+mbc
         xcell = xlower + (i-0.5d0)*dx
         do 15 j=1-mbc,my+mbc
         ycell = ylower + (j-0.5d0)*dy
         rewind(7)
            do 15 k=1-mbc,mz+mbc
            zcell = zlower + (k-0.5d0)*dz
c        # Isothermal Atmosphere          
c         
            if (loadprof.EQ.0) then
               aux(i,j,k,1) = density(zcell)
c           # Analytical wind function (if desired)
               if ((amp0 .EQ. 0) .AND. (amp1 .EQ. 0)) then
                  aux(i,j,k,2) = 0
               else
                  scale=(aux(i,j,1,1)/aux(i,j,k,1))**0.5d0
                  aux(i,j,k,2) = aux(i,j,k,1)*
     &               windfnx(xcell,ycell,zcell,0,scale)
               endif
               aux(i,j,k,3) = 0
               aux(i,j,k,4) = 0
               pres1= pressure(zcell)
               velx = aux(i,j,k,2)/aux(i,j,k,1)
               vely = 0
               velz = 0
               aux(i,j,k,5) = pres1/gamma1+
     &            .5*(velx**2+vely**2+velz**2)*aux(i,j,k,1)
c           # Reserve the following variables
               aux(i,j,k,6) = 0
               aux(i,j,k,7) = 0
               aux(i,j,k,8) = 0
               aux(i,j,k,9) = 0
               aux(i,j,k,10)= 0
c           # Default thermodynamic parameters
               aux(i,j,k,19) = gamma
               aux(i,j,k,18) = 287
c           # Set viscosity / conductivity
               dmolecular=dmolec0*(aux(i,j,1,1)/aux(i,j,k,1))
               dthermal=dtherm0*(aux(i,j,1,1)/aux(i,j,k,1))
               aux(i,j,k,11)=dmolecular
               aux(i,j,k,12)=dthermal*aux(i,j,k,19)
c
            else
c
c        # Otherwise, load vertical profile from profile.data
c
c           # Require exact height alignment
 10            read (7,*) height,dox,dnit2,dox2,aux(i,j,k,1),
     &             aux(i,j,k,5),dhyd,dne,windx,windy
               if (height.lt.(zcell/1000)) then
                  print*,height,'.LT.',(zcell/1000)
                  goto 10
               elseif (height.gt.(zcell/1000)) then
                  print*,height,'.GT.',(zcell/1000)
                  stop
               endif
c           # R / M    (approximately diatomic, gamma=constant)
               aux(i,j,k,18) = 8.31d3/
     &              ((dox*16+dnit2*28+dox2*32)/
     &               (dox+dnit2+dox2))
c           # Gamma
               aux(i,j,k,19) = (7*(dox2+dnit2)+5*dox)/
     &                         (5*(dox2+dnit2)+3*dox)
               temp  =aux(i,j,k,5)
               RoverM=aux(i,j,k,18)
               gamma =aux(i,j,k,19)
               gamma1=gamma-1

c           # Convert to MKS
               aux(i,j,k,1) = (dox*16.d0+dox2*32.d0+dnit2*28.d0)*1.d3
               aux(i,j,k,1) = aux(i,j,k,1)/(6.022d23)
c           # Analytical wind function (if desired)
               if ((amp0 .EQ. 0) .AND. (amp1 .EQ. 0)) then
                   aux(i,j,k,2) = aux(i,j,k,1)*windx
		   aux(i,j,k,3) = aux(i,j,k,1)*windy
               else
                  windx=amp0
                  scale=(aux(i,j,1,1)/aux(i,j,k,1))**0.5d0
                  aux(i,j,k,2) = 
     &            aux(i,j,k,1)*windfnx(xcell,ycell,zcell,0,scale)
                  aux(i,j,k,2) = aux(i,j,k,2)+windx*aux(i,j,k,1)
               endif
c              aux(i,j,k,3) = 0
               aux(i,j,k,4) = 0
               aux(i,j,k,5) = aux(i,j,k,1)*aux(i,j,k,5)*
     &                aux(i,j,k,18)/gamma1+
     &                .5*(aux(i,j,k,2)**2)/aux(i,j,k,1)+
     &                .5*(aux(i,j,k,3)**2)/aux(i,j,k,1)
c           # Reserve the following variables
               aux(i,j,k,6) = dox*16.d0*1.d3/(6.022d23*aux(i,j,k,1))
               aux(i,j,k,7) = dox2*32.d0*1.d3/(6.022d23*aux(i,j,k,1))
               aux(i,j,k,8) = dne/(dox+dox2+dnit2)
               aux(i,j,k,9) = dhyd/(dox+dox2+dnit2)
c           # Reserve for Ozone
               aux(i,j,k,10)= 0.d0
c           # Calculate diffusivities from MSIS Profile
               Cp=1/(1-1/gamma)*RoverM
               dmolec=(3.43*dnit2+4.03*dox2+3.90*dox)*1d-7
               dmolec=dmolec*(temp**(0.69))/(dox+dnit2+dox2)
               dtherm=(56.0*(dnit2+dox2)+75.9*dox)*1d-5
               dtherm=dtherm*(temp**(0.69))/(dox+dnit2+dox2)
c           # Assign values of diffusivity
               dmolecular=dmolec*(1/aux(i,j,k,1))
               dthermal=dtherm*(1/(aux(i,j,k,1)*Cp))
               aux(i,j,k,11)=dmolecular
               aux(i,j,k,12)=dthermal*aux(i,j,k,19)
            endif
c
c
c           # Set "virtual" gravity force 
c
            if (k.GT.(1-mbc)) then
               auxpl = (aux(i,j,k-1,19)-1.d0)
     &          * (aux(i,j,k-1,5)-0.5d0*(aux(i,j,k-1,2)**2
     &          + aux(i,j,k-1,3)**2 + aux(i,j,k-1,4)**2)/aux(i,j,k-1,1))
               auxpr = (aux(i,j,k,19)-1.d0)
     &          * (aux(i,j,k,5)-0.5d0*(aux(i,j,k,2)**2
     &          + aux(i,j,k,3)**2 + aux(i,j,k,4)**2)/aux(i,j,k,1))
               aux(i,j,k,17) = (auxpr-auxpl)/
     &           (.5d0*(aux(i,j,k,1)+aux(i,j,k-1,1)))
            endif
c
c
  15       continue
c
      if (id.EQ.0) then
        print*,'Max Diffusivities: ', aux(1,1,mz,11),aux(1,1,mz,12)
        dxyzmin=min(dx,dy,dz)
        difmax=max(aux(1,1,mz,11),aux(1,1,mz,12))
        difcfl=0.10d0
        gamma=aux(1,1,mz,19)
        gamma1=gamma-1
        pres = gamma1*(aux(1,1,mz,5)-0.5d0*(aux(1,1,mz,2)**2
     &       + aux(1,1,mz,3)**2 + aux(1,1,mz,4)**2)/aux(1,1,mz,1))
        dt=min(dx,dy,dz)/sqrt(gamma*pres/aux(1,1,mz,1))
        ndtr=max(int((dt*difmax/(dxyzmin*dxyzmin))/difcfl),1)
        print*,'Diffusive Burden: ', ndtr
      endif
c
c     # Close profile.data
c
      if (loadprof.EQ.1) then
         close(7)
      endif
c
      return
      end
