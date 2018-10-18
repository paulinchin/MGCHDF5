c     ============================================
      subroutine b4step3(maxmx,maxmy,maxmz,mbc,mx,my,mz,meqn,q,
     &            xlower,ylower,zlower,dx,dy,dz,t,dt,maux,aux)
c     ============================================
c
c     # called from claw3 before each call to step3.
c     # use to set time-dependent aux arrays or perform other tasks
c     # which must be done every time step.
c
c
      implicit double precision (a-h,o-z)
      dimension    q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, meqn)
      dimension  aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 
     &               1-mbc:maxmz+mbc, maux)
c
      logical varyR
      logical varyGamma
      logical varyMu
      data varyR     /.false./
      data varyGamma /.false./
      data varyMu    /.false./
c
c     
      do 15 k = 1-mbc, mz+mbc
        do 15 j = 1-mbc, my+mbc
          do 15 i = 1-mbc, mx+mbc
c
c         # Calculate Number Densities
            if (varyR.OR.varyGamma.OR.varyMu) then
              dox   = 6.022d23*(q(i,j,k,6)/16.d0)*q(i,j,k,1)*1.d-3
              dox2  = 6.022d23*(q(i,j,k,7)/32.d0)*q(i,j,k,1)*1.d-3
              dnit2 = 6.022d23*((1.d0-q(i,j,k,6)-q(i,j,k,7))/28.d0)
     &              * q(i,j,k,1)*1.d-3
            end if
c
c         # Specific Gas Contant R
            if (varyR) then
              aux(i,j,k,18) = 8.31d3 /
     &          ((dox*16.d0+dnit2*28.d0+dox2*32.d0) /
     &           (dox+dnit2+dox2))
            end if
c
c         # Ratio of Specific Heats Gamma
c         # Calculate similar to Hickey and Walterscheid ...
            if (varyGamma) then
              aux(i,j,k,19) = (7*(dox2+dnit2)+5*dox)/
     &                        (5*(dox2+dnit2)+3*dox)
            end if
c
c         # Ratio of Specific Heats Gamma
c         # Calculate similar to Hickey and Walterscheid ...
            if (varyGamma) then
              aux(i,j,k,19) = (7*(dox2+dnit2)+5*dox)/
     &                        (5*(dox2+dnit2)+3*dox)
            end if
c     
c         # Molecular Diffusivities / Conductivities
c         # Calculate from MSIS Profile
            if (varyMu) then
              temp = (aux(i,j,k,19)-1.d0)*(q(i,j,k,4)-
     &            0.5d0*(q(i,j,k,2)**2 + q(i,j,k,3)**2)/q(i,j,k,1))
     &            / (q(i,j,k,1)*aux(i,j,k,18))
              Cp = 1.d0/(1.d0-1.d0/aux(i,j,k,19))*aux(i,j,k,18)
              dmolec = (3.43*dnit2+4.03*dox2+3.90*dox)*1.d-7
              dmolec = dmolec*(temp**(0.69))/(dox+dnit2+dox2)
              dtherm = (56.0*(dnit2+dox2)+75.9*dox)*1.d-5
              dtherm = dtherm*(temp**(0.69))/(dox+dnit2+dox2)
              aux(i,j,k,11) = dmolec*(1.d0/q(i,j,k,1))
              aux(i,j,k,12) = (dtherm*(1.d0/(q(i,j,k,1)*Cp)))
     &                 * aux(i,j,k,19)
            end if
c
  15  continue
c
      return
      end
