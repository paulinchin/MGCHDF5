c
c
c     =======================================================
      subroutine ohsteady(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,
     &     xlower,ylower,zlower,dx,dy,dz,q,maux,aux,ohstep)
c     =======================================================
c
c
c
c  --------------------------------------------------------
c
c     Hydroxyl Solution for Atmospheric Model Extensions
c        by J. B. Snively, 2003-2017
c
c  --------------------------------------------------------

c
      implicit double precision (a-h,o-z)
c
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, meqn)
      dimension    aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)
      dimension b(0:9),a8(0:9),a9(1:9,0:8),suma9(1:9),
     &          a10(0:9),aa(1:9,1:6)
c
      data b/0,0,0,0,0,0,0.08,0.17,0.27,0.48/
      data a8/3.9,10.5,25,25,25,25,25,25,25,25/
      data ((a9(i,j),j=0,8),i=1,9)/
     *2,0,0,0,0,0,0,0,0,
     *0,4,0,0,0,0,0,0,0,
     *0,1,7,0,0,0,0,0,0,
     *0,1,2,10,0,0,0,0,0,
     *0,1,2,6,16,0,0,0,0,
     *1,1,3,6,11,22,0,0,0,
     *4,6,9,12,16,23,32,0,0,
     *4,6,8,10,14,19,25,33,0,
     *28,29,31,32,34,36,38,40,42/
      data suma9/2,4,8,13,25,44,102,119,310/
      data a10/0,0.58,1.0,1.7,3.0,5.2,9.1,16,70,48/
      data ((aa(i,j),j=1,6),i=1,9)/
     *22.74,0,0,0,0,0,
     *30.43,15.42,0,0,0,0,
     *28.12,40.33,2.032,0,0,0,
     *20.30,69.77,7.191,0.299,0,0,
     *11.05,99.42,15.88,1.315,0.051,0,
     *4,125.6,27.94,3.479,0.274,0.010,
     *2.34,145.1,42.91,7.165,0.847,0.063,
     *8.6,154.3,59.98,12.68,2.007,0.230,
     *23.72,148.9,78.64,19.94,4.053,0.620/
c
c     Find partial steady-state solution for
c     OH Airglow vibrational states v=0-9
c     [e.g., Adler-Golden, 1997]
c
      if (ohstep.EQ.0) then 
        hydmax = 0.d0
        kmin = int(80*(1000/dz))
        kmax = int(90*(1000/dz))
       do k = kmin, kmax
          dox   = 6.022d23*(q(1,1,k,6)/16.d0)*q(1,1,k,1)*1.d-3
          dnit2 = 6.022d23*((1.d0-q(1,1,k,6)-q(1,1,k,7))
     &            /28.d0)*q(1,1,k,1)*1.d-3
          dox2  = 6.022d23*(q(1,1,k,7)/32.d0)*q(1,1,k,1)*1.d-3
          dhyd  = q(1,1,k,9)*(dox+dox2+dnit2) 
          if (dhyd.GT.hydmax) then 
            hydmax  = dhyd
            hydzmax = zlower + (k-0.5d0)*dz
          endif
        enddo
      endif
      kmin = int(75*(1000/dz)) 
      kmax = int(125*(1000/dz))
      do i = 1, mx
        do j = 1, my
          do k = kmin, kmax
            zcell = zlower + (k-0.5d0)*dz
            if (q(i,j,k,6).LT.1.d-14) q(i,j,k,6)=1.d-14
            if (q(i,j,k,7).LT.1.d-14) q(i,j,k,7)=1.d-14
            if (q(i,j,k,9).LT.1.d-14) q(i,j,k,9)=1.d-14
            if (q(i,j,k,10).LT.1.d-14) q(i,j,k,10)=1.d-14
            dox   = 6.022d23*(q(i,j,k,6)/16.d0)*q(i,j,k,1)*1.d-3
            dnit2 = 6.022d23*((1.d0-q(i,j,k,6)-q(i,j,k,7))
     &              /28.d0)*q(i,j,k,1)*1.d-3
            dox2  = 6.022d23*(q(i,j,k,7)/32.d0)*q(i,j,k,1)*1.d-3   
            dhyd  = q(i,j,k,9)*(dox+dox2+dnit2)
            dox3  = q(i,j,k,10)*(dox+dox2+dnit2)
            gamma = aux(i,j,k,19)
            gamma1= gamma-1.d0
c        # Convert energy to temperature
            tempk = gamma1*(q(i,j,k,5)-0.5d0*(q(i,j,k,2)**2+
     &      q(i,j,k,3)**2+q(i,j,k,4)**2)/q(i,j,k,1))/
     &      (q(i,j,k,1)*aux(i,j,k,18))
            rk1   = (1.4d-10)*dexp(-470/tempk)
            rk3N2 = (5.7d-34)*((300/tempk)**(2.62))
            rk3O2 = (5.96d-34)*((300/tempk)**(2.37))
            rk8 = 1d-11
            rk9 = 1d-13
            rk10= 1d-14
c        # First Step only...
          if (ohstep.EQ.0) then
c        # Approx. Ozone Steady-State Profile
            dhyd1 = hydmax*dexp(-(zcell-hydzmax)/8620)
            rL = dhyd1*rk1
            rP = dox*dox2*(dox2*rk3O2+dnit2*rk3N2)
            dox3 = rP/rL
            q(i,j,k,10) = dox3/(dox+dox2+dnit2)
          endif
c        # Calculate [OH(v=9)]
            rP=rk1*b(9)*dhyd*dox3
            rL=rk8*a8(9)*dox+rk9*suma9(9)*dox2+rk10*a10(9)*dnit2
     &         +aa(9,1)+aa(9,2)+aa(9,3)+aa(9,4)+aa(9,5)+aa(9,6)
            aux(i,j,k,29)=rP/rL
c        # Calculate [OH(v=8)]
            rP=rk1*b(8)*dhyd*dox3+aux(i,j,k,29)*aa(9,1)
     &         +rk9*(a9(9,8)*aux(i,j,k,29))*dox2
     &         +rk10*a10(9)*dnit2*aux(i,j,k,29)
            rL=rk8*a8(8)*dox+rk9*suma9(8)*dox2+rk10*a10(8)*dnit2
     &         +aa(8,1)+aa(8,2)+aa(8,3)+aa(8,4)+aa(8,5)+aa(8,6)
            aux(i,j,k,28)=rP/rL
c        # Calculate [OH(v=7)]
            rP=rk1*b(7)*dhyd*dox3+aux(i,j,k,29)*aa(9,2)
     &         +aux(i,j,k,28)*aa(8,1)+rk9*(a9(9,7)*aux(i,j,k,29)
     &         +a9(8,7)*aux(i,j,k,28))*dox2
     &         +rk10*a10(8)*dnit2*aux(i,j,k,28)
            rL=rk8*a8(7)*dox+rk9*suma9(7)*dox2+rk10*a10(7)*dnit2
     &         +aa(7,1)+aa(7,2)+aa(7,3)+aa(7,4)+aa(7,5)+aa(7,6)
            aux(i,j,k,27)=rP/rL
c        # Calculate [OH(v=6)]
            rP=rk1*b(6)*dhyd*dox3+aux(i,j,k,29)*aa(9,3)
     &         +aux(i,j,k,28)*aa(8,2)+aux(i,j,k,27)*aa(7,1)
     &         +rk9*(a9(9,6)*aux(i,j,k,29)+a9(8,6)*aux(i,j,k,28)
     &         +a9(7,6)*aux(i,j,k,27))*dox2
     &         +rk10*a10(7)*dnit2*aux(i,j,k,27)
            rL=rk8*a8(6)*dox+rk9*suma9(6)*dox2+rk10*a10(6)*dnit2
     &         +aa(6,1)+aa(6,2)+aa(6,3)+aa(6,4)+aa(6,5)+aa(6,6)
            aux(i,j,k,26)=rP/rL
c        # Calculate [OH(v=5)]
            rP=aux(i,j,k,29)*aa(9,4)+aux(i,j,k,28)*aa(8,3)
     &         +aux(i,j,k,27)*aa(7,2)+aux(i,j,k,26)*aa(6,1)
     &         +rk9*(a9(9,5)*aux(i,j,k,29)+a9(8,5)*aux(i,j,k,28)
     &         +a9(7,5)*aux(i,j,k,27)+a9(6,5)*aux(i,j,k,26))*dox2
     &         +rk10*a10(6)*dnit2*aux(i,j,k,26)  
            rL=rk8*a8(5)*dox+rk9*suma9(5)*dox2+rk10*a10(5)*dnit2
     &         +aa(5,1)+aa(5,2)+aa(5,3)+aa(5,4)+aa(5,5)
            aux(i,j,k,25)=rP/rL
c        # Calculate [OH(v=4)]
            rP=aux(i,j,k,29)*aa(9,5)+aux(i,j,k,28)*aa(8,4)
     &         +aux(i,j,k,27)*aa(7,3)+aux(i,j,k,26)*aa(6,2)
     &         +aux(i,j,k,25)*aa(5,1)+rk9*(a9(9,4)*aux(i,j,k,29)
     &         +a9(8,4)*aux(i,j,k,28)+a9(7,4)*aux(i,j,k,27)
     &         +a9(6,4)*aux(i,j,k,26)+a9(5,4+1)*aux(i,j,k,25))*dox2
     &         +rk10*a10(5)*dnit2*aux(i,j,k,25)
            rL=rk8*a8(4)*dox+rk9*suma9(4)*dox2+rk10*a10(4)*dnit2
     &         +aa(4,1)+aa(4,2)+aa(4,3)+aa(4,4)
            aux(i,j,k,24)=rP/rL
c        # Calculate [OH(v=3)]
            rP=aux(i,j,k,29)*aa(9,6)+aux(i,j,k,28)*aa(8,5)
     &         +aux(i,j,k,27)*aa(7,4)+aux(i,j,k,26)*aa(6,3)
     &         +aux(i,j,k,25)*aa(5,2)+aux(i,j,k,24)*aa(4,1)
     &         +rk9*(a9(9,3)*aux(i,j,k,29)+a9(8,3)*aux(i,j,k,28)
     &         +a9(7,3)*aux(i,j,k,27)+a9(6,3)*aux(i,j,k,26)
     &         +a9(5,3)*aux(i,j,k,25)+a9(4,3)*aux(i,j,k,24))*dox2
     &         +rk10*a10(4)*dnit2*aux(i,j,k,24)
            rL=rk8*a8(3)*dox+rk9*suma9(3)*dox2+rk10*a10(3)*dnit2
     &         +aa(3,1)+aa(3,2)+aa(3,3)
            aux(i,j,k,23)=rP/rL
c        # Calculate [OH(v=2)]
            rP=aux(i,j,k,28)*aa(8,6)+aux(i,j,k,27)*aa(7,5)
     &         +aux(i,j,k,26)*aa(6,4)+aux(i,j,k,25)*aa(5,3)
     &         +aux(i,j,k,24)*aa(4,2)+aux(i,j,k,23)*aa(3,1)
     &         +rk9*(a9(9,2)*aux(i,j,k,29)+a9(8,2)*aux(i,j,k,28)
     &         +a9(7,2)*aux(i,j,k,27)+a9(6,2)*aux(i,j,k,26)
     &         +a9(5,2)*aux(i,j,k,25)+a9(4,2)*aux(i,j,k,24)
     &         +a9(3,2)*aux(i,j,k,23))*dox2
     &         +rk10*a10(3)*dnit2*aux(i,j,k,23)
            rL=rk8*a8(2)*dox+rk9*suma9(2)*dox2
     &         +rk10*a10(2)*dnit2+(aa(2,1)+aa(2,2))
            aux(i,j,k,22)=rP/rL
c        # Calculate [OH(v=1)]
            rP=aux(i,j,k,27)*aa(7,6)+aux(i,j,k,26)*aa(6,5)
     &         +aux(i,j,k,25)*aa(5,4)+aux(i,j,k,24)*aa(4,3)
     &         +aux(i,j,k,23)*aa(3,2)+aux(i,j,k,22)*aa(2,1)
     &         +rk9*(a9(9,1)*aux(i,j,k,29)+a9(8,1)*aux(i,j,k,28)
     &         +a9(7,1)*aux(i,j,k,27)+a9(6,1)*aux(i,j,k,26)
     &         +a9(5,1)*aux(i,j,k,25)+a9(4,1)*aux(i,j,k,24)
     &         +a9(3,1)*aux(i,j,k,23)+a9(2,1)*aux(i,j,k,22))*dox2
     &         +rk10*a10(2)*dnit2*aux(i,j,k,22)
            rL=rk8*a8(1)*dox+rk9*suma9(1)*dox2
     &         +rk10*a10(1)*dnit2+aa(1,1)
            aux(i,j,k,21)=rP/rL
c        # Calculate [OH(v=0)]
            rP=aux(i,j,k,26)*aa(6,6)+aux(i,j,k,25)*aa(5,5)
     &         +aux(i,j,k,24)*aa(4,4)+aux(i,j,k,23)*aa(3,3)
     &         +aux(i,j,k,22)*aa(2,2)+aux(i,j,k,21)*aa(1,1)
     &         +rk9*(a9(9,0)*aux(i,j,k,29)+a9(8,0)*aux(i,j,k,28)
     &         +a9(7,0)*aux(i,j,k,27)+a9(6,0)*aux(i,j,k,26)
     &         +a9(5,0)*aux(i,j,k,25)+a9(4,0)*aux(i,j,k,24)
     &         +a9(3,0)*aux(i,j,k,23)+a9(2,0)*aux(i,j,k,22)
     &         +a9(1,0)*aux(i,j,k,21))*dox2
     &         +rk10*a10(1)*dnit2*aux(i,j,k,21)
            rL=rk8*a8(0)*dox
            aux(i,j,k,20)=rP/rL
c        # Set Hydrogen Steady-State
            rP = 0.d0
c        # Production and loss rates
            do mvib = 0, 9
               rP = rP+aux(i,j,k,20+mvib)*dox*rk8*a8(mvib)
            enddo
            rL = dhyd*dox3*rk1
            if (ohstep.EQ.0) then
c        #  ... define initial tendency ...
c        # Hydrogen
              aux(i,j,k,30) = rP-rL
c        # Ozone
              rP = dox*dox2*(dox2*rk3O2+dnit2*rk3N2)
c        #  ... define initial tendency ...
              aux(i,j,k,31) = rP-rL
            else
c        #  ... time step P and L for H and O3 ...
c        # Hydrogen
              dhyd = dhyd+ohstep*((rP-rL)-aux(i,j,k,30))
              q(i,j,k,9) = dhyd/(dox+dox2+dnit2)
c        # Ozone
              rP   = dox*dox2*(dox2*rk3O2+dnit2*rk3N2)
              dox3 = dox3+ohstep*((rP-rL)-aux(i,j,k,31))
              q(i,j,k,10) = dox3/(dox+dox2+dnit2)
            endif 
          enddo
        enddo
      enddo
c
c
c
      return
      end

