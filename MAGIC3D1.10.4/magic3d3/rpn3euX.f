c
c     ==================================================================
      subroutine rpn3(ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &		      auxl,auxr,maux,dtdx,fwave,s,amdq,apdq)
c     ==================================================================
c
c  --------------------------------------------------------
c
c     Riemann Solution for Atmospheric Model Extensions
c        by J. B. Snively, 2003-2017,
c        based on work by LeVeque, Langseth, et al.
c
c  --------------------------------------------------------
c
c     # Roe-solver for the Euler equations, f-wave, with conservative tracers.
c     # Flux limiting is applied with this routine, so that only 3 total waves 
c     # are needed, including for tracers.
c
c     # Tracers are solved conservatively based on edge (Roe) velocities.
c     # The number of tracers = meqn-5 (Each lumped into wave 2)
c     # The first tracers nspecadv are solved via simple advection,
c     # the next tracers meqn-5-nspecadv are solved via continuity.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixyz=1
c     #                               the y-direction if ixyz=2.
c     #                               the z-direction if ixyz=3.
c     # On output, fwave contains the f-waves, s the speeds,
c     # and amdq, apdq the left-going and right-going flux differences,
c     # respectively.  
c
c     # A limiter is applied for each tracer species individually;
c     # this is because limiting them together does not guarantee stability.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision (a-h,o-z)
c
      dimension fwave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)
      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension auxl(1-mbc:maxm+mbc, maux)
      dimension auxr(1-mbc:maxm+mbc, maux)
c
c     # Now, to enable CFL-dependent limiting... (modified Clawpack)
c
      dimension  dtdx(1-mbc:maxm+mbc)
c
c     local arrays and common block (comroe is passed to rpt3eu)
c     ------------
      parameter (maxmrp = 2300)
      dimension delta(1-mbc:maxm+mbc,5)
c
      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
      common /comroe/ u2v2w2(-1:maxmrp),
     &       u(-1:maxmrp),v(-1:maxmrp),w(-1:maxmrp),enth(-1:maxmrp),
     &       a(-1:maxmrp),g1a2(-1:maxmrp),euv(-1:maxmrp)
c
      common /options/ loadprof, nspecadv
      logical efix, laxliu, thirdorder
c
      data efix /.true./       !# use entropy fix for transonic rarefactions
      data laxliu /.true./     !# use more-expensive Lax-Liu limiting
      data thirdorder /.true./ !# use third order CFL-dependent method
c
c     # CFL Superbee beta-theta limiter parameters
c     # (Theta = 1, Beta = 2/3 Recommended)
c     
      theta = 1.0d0
      beta  = 0.666666666666666d0
c
c
c     # check hard-coded limits and settings
c
      if (-1.GT.1-mbc .OR. maxmrp .LT. maxm+mbc) then
	 write(6,*) 'need to increase maxmrp in rpn3eu'
         write(6,*) 'maxm+mbc=',maxm+mbc
	 stop
	 endif
c
c     # check that mwaves=3, because this solver combines waves at u(i)
c
      if (mwaves .NE. 3) then
         write(6,*) 'need mwaves=3 for this Riemann solver'
         stop
         endif
c
c     # set mu to point to  the component of the system that corresponds
c     # to momentum in the direction of this slice, mv and mw to the
c     # orthogonal momentums:
c
      if (ixyz .EQ. 1) then
	  mu = 2
	  mv = 3
          mw = 4
      else if (ixyz .EQ. 2) then
	  mu = 3
	  mv = 4
          mw = 2
      else
          mu = 4
          mv = 2
          mw = 3
      endif
c
c     # note that notation for u,v, and w reflects assumption that the
c     # Riemann problems are in the x-direction with u in the normal
c     # direction and v and w in the orthogonal directions, but with the
c     # above definitions of mu, mv, and mw the routine also works with
c     # ixyz=2 and ixyz = 3
c     # and returns, for example, f0 as the Godunov flux g0 for the
c     # Riemann problems u_t + g(u)_y = 0 in the y-direction.
c
c
c     # Compute the Roe-averaged variables needed in the Roe solver.
c     # These are stored in the common block comroe since they are
c     # later used in routine rpt3eu to do the transverse wave
c     # splitting.
c
c
c     # Loop over rest ...
c
      do 20 i = 2-mbc, mx+mbc
c        # Check for negative densities, however unlikely...
         if (qr(i-1,1) .le. 0.d0 .or. ql(i,1) .le. 0.d0) then
            write(*,*) i, mu, mv, mw
            write(*,990) (qr(i-1,j),j=1,5)
            write(*,990) (ql(i,j),j=1,5)
 990        format(5e12.4)
            if (ixyz .eq. 1) 
     &         write(6,*) '*** rho .le. 0 in x-sweep at ',i,jcom,kcom
            if (ixyz .eq. 2) 
     &         write(6,*) '*** rho .le. 0 in y-sweep at ',icom,i,kcom
            if (ixyz .eq. 3) 
     &         write(6,*) '*** rho .le. 0 in z-sweep at ',icom,jcom,i
            write(6,*) 'stopped with rho < 0...'
           stop
         endif
	 rhsqrtl = dsqrt(qr(i-1,1))
	 rhsqrtr = dsqrt(ql(i,1))
	 rhsq2   = rhsqrtl + rhsqrtr
         gamma   = (rhsqrtr*auxl(i,19)+rhsqrtl*auxr(i-1,19))/rhsq2
         gamma1  = gamma-1.d0
	 pl = (auxr(i-1,19)-1.d0)*(qr(i-1,5) - 0.5d0*(qr(i-1,mu)**2 +
     &	       qr(i-1,mv)**2 + qr(i-1,mw)**2)/qr(i-1,1))
	 pr = (auxl(i,19)-1.d0)*(ql(i,5) - 0.5d0*(ql(i,mu)**2 +
     &	       ql(i,mv)**2 + ql(i,mw)**2)/ql(i,1))
	 u(i) = (qr(i-1,mu)/rhsqrtl + ql(i,mu)/rhsqrtr) / rhsq2 
	 v(i) = (qr(i-1,mv)/rhsqrtl + ql(i,mv)/rhsqrtr) / rhsq2
	 w(i) = (qr(i-1,mw)/rhsqrtl + ql(i,mw)/rhsqrtr) / rhsq2
	 enth(i) = (((qr(i-1,5)+pl)/rhsqrtl
     &	         + (ql(i,5)+pr)/rhsqrtr)) / rhsq2
	 u2v2w2(i) = u(i)**2 + v(i)**2 + w(i)**2
         a2 = gamma1*(enth(i) - .5d0*u2v2w2(i))
c        # Further error checking, also unlikely...

c	print*,'ARE WE JERE?'
c	if(id==0) then
c	print*,'WHAT IS ABOUT A2?'
c	print*,a2
c	end if

         if (a2 .le. 0.d0) then
            if (ixyz .eq. 1) 
     &         write(6,*) '*** a2 .le. 0 in x-sweep at ',i,jcom,kcom
            if (ixyz .eq. 2) 
     &         write(6,*) '*** a2 .le. 0 in y-sweep at ',icom,i,kcom
            if (ixyz .eq. 3) 
     &         write(6,*) '*** a2 .le. 0 in z-sweep at ',icom,jcom,i
            write(6,*) 'stopped with a2 < 0...'
c           stop
         endif
         a(i) = dsqrt(a2)
	 g1a2(i) = gamma1 / a2
	 euv(i) = enth(i) - u2v2w2(i)
c
c
c     # Now split the jump in q1d at each interface into waves
c
c     # Find b1 thru b5, the coefficients of the 5 eigenvectors:
         ur = ql(i,mu)/ql(i,1)
         ul = qr(i-1,mu)/qr(i-1,1)
         delta(i,1) = ql(i,mu) - qr(i-1,mu)
         delta(i,2) = (ql(i,mu)*ur + pr) - (qr(i-1,mu)*ul + pl)
         delta(i,3) = ur*ql(i,mv) - ul*qr(i-1,mv)
         delta(i,4) = ur*ql(i,mw) - ul*qr(i-1,mw)
         delta(i,5) = (ql(i,5)+pr)*ur - (qr(i-1,5)+pl)*ul
c     # Apply gravity for vertical direction ...
         if (ixyz.eq.3) then
c        # Based on cell centers
            delta(i,2) = delta(i,2) - 
     &                0.5d0*(ql(i,1)+qr(i-1,1))*auxl(i,17)
            delta(i,5) = delta(i,5) -
     &                0.5d0*(ql(i,mu)+qr(i-1,mu))*auxl(i,17)
         endif
         b4 = g1a2(i) * (euv(i)*delta(i,1)
     &      + u(i)*delta(i,2) + v(i)*delta(i,3) + w(i)*delta(i,4)
     &      - delta(i,5))
         b2 = delta(i,3) - v(i)*delta(i,1)
         b3 = delta(i,4) - w(i)*delta(i,1)
         b5 = (delta(i,2) + (a(i)-u(i))*delta(i,1)  
     &      - a(i)*b4) / (2.d0*a(i))
         b1 = delta(i,1) - b4 - b5
c
c     # Compute the waves.
c     # Note that the 2-, 3- and 4-wave travel at the same speed,
c     # and are combined. They are limited separately, however -- 
c     # This is a critical correction for stability.
c     # The 5-wave is stored in fwave(.,.,3). 
c     # Tracers are in fwave(.,.,2).
c
c     # Acoustic:
         fwave(i,1,1)  = b1
         fwave(i,mu,1) = b1*(u(i)-a(i))
         fwave(i,mv,1) = b1*v(i)
         fwave(i,mw,1) = b1*w(i)
         fwave(i,5,1)  = b1*(enth(i) - u(i)*a(i))
         s(i,1) = u(i)-a(i)
c
c     # Shear + Entropy (+Tracers):
         fwave(i,1,2)  = b4
         fwave(i,mu,2) = b4*u(i)
         fwave(i,mv,2) = b2+b4*v(i)
         fwave(i,mw,2) = b3+b4*w(i)
         fwave(i,5,2)  = b2*v(i)+b3*w(i)+b4*0.5d0*u2v2w2(i)
         s(i,2) = u(i)
c
c     # Acoustic:
         fwave(i,1,3)  = b5
         fwave(i,mu,3) = b5*(u(i)+a(i))
         fwave(i,mv,3) = b5*v(i)
         fwave(i,mw,3) = b5*w(i)
         fwave(i,5,3)  = b5*(enth(i)+u(i)*a(i))
         s(i,3) = u(i)+a(i)
c
c     # Additional nspec waves added for tracers or mixing ratios;
c     # assume conserved/advected species coupled at Roe velocity.
c
c     # Specify fluxes/fwaves for meqn-5 species:
         do 15 m = 6, meqn
c        # Total tracers = meqn-5
c        # Total advected tracers = nspecadv
c        # Total conserved tracers = meqn-5-nspecadv
c        # Set neighboring fwaves to zero
            fwave(i,m,1) = 0.d0
            fwave(i,m,3) = 0.d0
c        # Select conservative vs. pure advection
            if ((m-5) .LE. nspecadv) then
c           # Set fwave without velocity differencing
               fwave(i,m,2) = s(i,2)*(ql(i,m) - qr(i-1,m))
            else
c           # Set flux difference at cell-center velocity:
               fwave(i,m,2) = ur*ql(i,m)-ul*qr(i-1,m)
            endif
   15    continue
c
   20 continue
c
c
c     # compute flux differences amdq and apdq.
c     ---------------------------------------
c
      if (efix) goto 110
c
c     # no entropy fix
c     ----------------
c
c     # amdq = SUM fwave   over left-going waves
c     # apdq = SUM fwave   over right-going waves
c
      do 100 m = 1, meqn
         do 100 i = 2-mbc, mx+mbc
            amdq(i,m) = 0.d0
            apdq(i,m) = 0.d0
            do 90 mws = 1, 3
               if (s(i,mws) .LT. 0.d0) then
                  amdq(i,m) = amdq(i,m) + fwave(i,m,mws)
               else
                  apdq(i,m) = apdq(i,m) + fwave(i,m,mws)
               endif
   90       continue
  100 continue
      go to 300
c
c-----------------------------------------------------
c
  110 continue
c
c     # With entropy fix
c     ------------------
c
c    # compute flux differences amdq and apdq.
c    # First compute amdq as sum of fwave for left going waves.
c    # Incorporate entropy fix by adding a modified fraction of wave
c    # if s should change sign.
c
      do 200 i = 2-mbc, mx+mbc
c
c        # check 1-wave:
c        ---------------
c
         gamma  = auxr(i-1,19)
         gamma1 = gamma-1.d0
	 rhoim1 = qr(i-1,1)
	 pim1   = gamma1*(qr(i-1,5) - 0.5d0*(qr(i-1,mu)**2
     &          + qr(i-1,mv)**2 + qr(i-1,mw)**2) / rhoim1)
	 cim1   = dsqrt(gamma*pim1/rhoim1)
	 s0     = qr(i-1,mu)/rhoim1 - cim1     !# u-c in left state (cell i-1)
c
c
c     # check for fully supersonic case:
         if (s0.GE.0.d0 .AND. s(i,1).GT.0.d0)then
c        # everything is right-going
            do 115 m = 1, meqn
               amdq(i,m) = 0.d0
  115       continue
            go to 200
         endif
c
         rho1   = qr(i-1,1) + fwave(i,1,1)/s(i,1)
         rhou1  = qr(i-1,mu) + fwave(i,mu,1)/s(i,1)
         rhov1  = qr(i-1,mv) + fwave(i,mv,1)/s(i,1)
         rhow1  = qr(i-1,mw) + fwave(i,mw,1)/s(i,1)
         en1    = qr(i-1,5) + fwave(i,5,1)/s(i,1)
         p1     = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2 +
     &                rhow1**2)/rho1)
         c1     = dsqrt(gamma*p1/rho1)
         s1     = rhou1/rho1 - c1  !# u-c to right of 1-wave
         if (s0.LT.0.d0 .AND. s1.GT.0.d0) then
c        # transonic rarefaction in the 1-wave
	    sfract = s0 * (s1-s(i,1)) / (s1-s0)
	 else if (s(i,1) .LT. 0.d0) then
c	 # 1-wave is leftgoing
	    sfract = s(i,1)
         else
c	 # 1-wave is rightgoing
            sfract = 0.d0   !# this shouldn't happen since s0 < 0
         endif
         do 120 m = 1, meqn
            amdq(i,m) = sfract*fwave(i,m,1)/s(i,1)
  120    continue
c
c        # check 2-wave:
c        ---------------
c
         if (s(i,2) .GE. 0.d0) go to 200  !# 2-,3- and 4- waves are rightgoing
         do 140 m = 1, meqn
            amdq(i,m) = amdq(i,m) + fwave(i,m,2)
  140    continue
c
c        # check 5-wave:
c        ---------------
c
         gamma  = auxl(i,19)
         gamma1 = gamma-1.d0
         rhoi   = ql(i,1)
         pi     = gamma1*(ql(i,5) - 0.5d0*(ql(i,mu)**2
     &          + ql(i,mv)**2 + ql(i,mw)**2) / rhoi)
         ci     = dsqrt(gamma*pi/rhoi)
         s3     = ql(i,mu)/rhoi + ci     !# u+c in right state  (cell i)
c
         rho2   = ql(i,1) - fwave(i,1,3)/s(i,3)
         rhou2  = ql(i,mu) - fwave(i,mu,3)/s(i,3)
         rhov2  = ql(i,mv) - fwave(i,mv,3)/s(i,3)
         rhow2  = ql(i,mw) - fwave(i,mw,3)/s(i,3)
         en2    = ql(i,5) - fwave(i,5,3)/s(i,3)
         p2     = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2 +
     &                    rhow2**2)/rho2)
         c2     = dsqrt(gamma*p2/rho2)
         s2     = rhou2/rho2 + c2   !# u+c to left of 5-wave
         if (s2 .LT. 0.d0 .AND. s3.GT.0.d0 ) then
c        # transonic rarefaction in the 5-wave
            sfract = s2 * (s3-s(i,3)) / (s3-s2) 
         else if (s(i,3) .LT. 0.d0) then
c        # 5-wave is leftgoing
            sfract = s(i,3)
         else
c        # 5-wave is rightgoing
            go to 200
         endif
c
         do 160 m = 1, meqn
            amdq(i,m) = amdq(i,m) + sfract*fwave(i,m,3)/s(i,3)
  160    continue
  200 continue
c
c     # compute the rightgoing flux differences:
c     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
c
      do 220 m = 1, meqn
         do 220 i = 2-mbc, mx+mbc
            df = 0.d0
            do 210 mws = 1, 3
               df = df + fwave(i,m,mws)
  210       continue
         apdq(i,m) = df - amdq(i,m)
  220 continue
c
c
  300 continue
c
c
c     # Apply flux limiters to species individually:
c     # Reduce compute time by sharing wave 3; ensure
c     # stability by separating limiter -- must also set
c     # mthlim = 0 or limiters will be applied twice.
c     ---------------------------------------------------
c
      do 320 m = 6, meqn
         fwavec = 0.d0
         do 310 i = 0, mx+1
            fwavel = fwavec
            fwavec = fwave(i,m,2)
            fwaver = fwave(i+1,m,2)
            if (i .EQ. 0) goto 310
            if (fwavec .EQ. 0.d0) goto 310
            if (s(i,2) .GT. 0.d0) then
               r = fwavel/fwavec
               CFL2 = s(i,2)*dtdx(i)
            else
               r = fwaver/fwavec
               CFL2 = -s(i,2)*dtdx(i+1)
            endif
c           # Apply limiter to fwave
            if ((m-5).LE.nspecadv) then
c           # Advection Equation gets same limiter as Entropy Wave
               if (thirdorder) then
c              # CFL-Superbee Beta = 2/3
                  wlimitr = cflsupbee(theta,beta,r,CFL2)
               else
c              # MC Limiter (CFL independent, Good Choice)
                  c = (1.d0 + 1.d0*r)/2.d0
                  wlimitr = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r))
               endif
            else
c           # Continuity Equation gets minmod, as a safer choice
               wlimitr = dmax1(0.d0,dmin1(1.d0, r))
            endif
c           # Apply limiter to fwave
            fwave(i,m,2) = wlimitr * fwave(i,m,2)
  310    continue
  320 continue
c
c
c     # Apply Lax-Liu style limiters to other fwaves:
c     ---------------------------------------------------
c
      b1 = 0.d0
      b2 = 0.d0
      b3 = 0.d0
      b4 = 0.d0
      b5 = 0.d0
      do 330 i = 0, mx+1
         if (.NOT.(laxliu)) then
c        # Store left coefficients from previous
            b1l = b1
            b2l = b2
            b3l = b3
            b4l = b4
            b5l = b5
c        # Store right coefficients from fwaves
            b1r = fwave(i+1,1,1)
            b4r = fwave(i+1,1,2)
            b2r = fwave(i+1,mv,2)-b4r*v(i+1)
            b3r = fwave(i+1,mw,2)-b4r*w(i+1)
            b5r = fwave(i+1,1,3)
         endif
c     # Obtain coefficients from fwaves
         b1 = fwave(i,1,1)
         b4 = fwave(i,1,2)
         b2 = fwave(i,mv,2)-b4*v(i)
         b3 = fwave(i,mw,2)-b4*w(i)
         b5 = fwave(i,1,3)
         if (i .EQ. 0) goto 330
c
c     # Acoustic Wave 1
         if (b1 .NE. 0.d0) then
            Id1 = int(i-dsign(1.d0,s(i,1)))
c        # Solve coefficient, find ratio r1
            if (laxliu) then
               b4Id = g1a2(i) * (euv(i)*delta(Id1,1)
     &              + u(i)*delta(Id1,2) + v(i)*delta(Id1,3) 
     &              + w(i)*delta(Id1,4) - delta(Id1,5))
               b5Id = (delta(Id1,2) + (a(i)-u(i))*delta(Id1,1) 
     &              - a(i)*b4Id) / (2.d0*a(i))
               b1Id = delta(Id1,1) - b4Id - b5Id
            else  !# Classic Clawpack Limiter
               if (s(i,1) .GT. 0.d0) then
                  b1Id = b1l
               else
                  b1Id = b1r
               endif
            endif
            r1 = b1Id/b1
            if (thirdorder) then
               CFL1 = dabs(s(i,1)*dtdx(max(i,Id1)))
c           # SuperPower Limiter for Wave 1:
               wlimiter1 = suppow(r1,CFL1)
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r1)/2.d0
               wlimiter1 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r1))
            endif
         endif
c
c     # Shear Wave 2
         if (b2 .NE. 0.d0) then
            Id2 = int(i-dsign(1.d0,s(i,2)))
c        # Solve coefficient, find ratio r2
            if (laxliu) then
               b2Id = delta(Id2,3) - v(i)*delta(Id2,1)
            else  !# Classic Clawpack Limiter
               if (s(i,2) .GT. 0.d0) then
                  b2Id = b2l
               else
                  b2Id = b2r
               endif
            endif
            r2 = b2Id/b2
            if (thirdorder) then
            CFL2 = dabs(s(i,2)*dtdx(max(i,Id2)))
c           # CFL-Superbee beta=2/3 Limiter for Wave 2:
               wlimiter2 = cflsupbee(theta,beta,r2,CFL2)
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r2)/2.d0
               wlimiter2 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r2))
            endif
         endif
c
c     # Shear Wave 3
         if (b3 .NE. 0.d0) then
            Id3 = int(i-dsign(1.d0,s(i,2)))
c        # Solve coefficient, find ratio r3
            if (laxliu) then
               b3Id = delta(Id3,4) - w(i)*delta(Id3,1)
            else  !# Classic Clawpack Limiter
               if (s(i,2) .GT. 0.d0) then
                  b3Id = b3l
               else
                  b3Id = b3r
               endif
            endif
            r3 = b3Id/b3
            if (thirdorder) then
               CFL3 = dabs(s(i,2)*dtdx(max(i,Id3)))
c           # CFL-Superbee beta=2/3 Limiter for Wave 3:
               wlimiter3 = cflsupbee(theta,beta,r3,CFL3)
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r3)/2.d0
               wlimiter3 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r3))
            endif
         endif
c
c     # Entropy Wave 4
         if (b4 .NE. 0.d0) then
            Id4 = int(i-dsign(1.d0,s(i,2)))
c        # Solve coefficient, find ratio r4
            if (laxliu) then
               b4Id = g1a2(i) * (euv(i)*delta(Id4,1)
     &              + u(i)*delta(Id4,2) + v(i)*delta(Id4,3)
     &              + w(i)*delta(Id4,4) - delta(Id4,5))
            else  !# Classic Clawpack Limiter
               if (s(i,2) .GT. 0.d0) then
                  b4Id = b4l
               else
                  b4Id = b4r
               endif
            endif
            r4 = b4Id/b4
            if (thirdorder) then
               CFL4 = dabs(s(i,2)*dtdx(max(i,Id4)))
c           # CFL-Superbee beta=2/3 Limiter for Wave 4:
               wlimiter4 = cflsupbee(theta,beta,r4,CFL4)
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r4)/2.d0
               wlimiter4 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r4))
            endif
         endif
c
c     # Acoustic Wave 5
         if (b5 .NE. 0.d0) then
            Id5 = int(i-dsign(1.d0,s(i,3)))
c        # Solve coefficient, find ratio r5
            if (laxliu) then
               b4Id = g1a2(i) * (euv(i)*delta(Id5,1)
     &              + u(i)*delta(Id5,2) + v(i)*delta(Id5,3)
     &              + w(i)*delta(Id5,4) - delta(Id5,5))
               b5Id = (delta(Id5,2) + (a(i)-u(i))*delta(Id5,1) 
     &              - a(i)*b4Id) / (2.d0*a(i))
            else  !# Classic Clawpack Limiter
               if (s(i,3) .GT. 0.d0) then
                  b5Id = b5l
               else
                  b5Id = b5r
               endif
            endif      
            r5 = b5Id/b5
            if (thirdorder) then
               CFL5 = dabs(s(i,3)*dtdx(max(i,Id5)))
c           # CFL-Superbee beta=2/3 Limiter for Wave 5:
               wlimiter5 = suppow(r5,CFL5)
            else
c           # MC Limiter (CFL independent, Good Choice)
               c = (1.d0 + 1.d0*r5)/2.d0
               wlimiter5 = dmax1(0.d0,dmin1(c, 2.d0, 2.d0*r5))
            endif
         endif
c        
c     # Recompute the limited fwaves, collapsing into three ...
c
c     # Acoustic Waves 1 and 5:
         do m = 1, 5
            fwave(i,m,1) = fwave(i,m,1)*wlimiter1
            fwave(i,m,3) = fwave(i,m,3)*wlimiter5
         enddo
c     # Shear + Entropy Waves 2, 3, and 4:
         b2ltd = wlimiter2*b2
         b3ltd = wlimiter3*b3
         b4ltd = wlimiter4*b4
         fwave(i,1,2)  = b4ltd
         fwave(i,mu,2) = b4ltd*u(i)
         fwave(i,mv,2) = b2ltd+b4ltd*v(i)
         fwave(i,mw,2) = b3ltd+b4ltd*w(i)
         fwave(i,5,2)  = b2ltd*v(i)+b3ltd*w(i) + b4ltd*0.5d0*u2v2w2(i)
c
  330 continue

  500 continue

      return
      end
c
c     =====================================================
      double precision function cflsupbee(theta,beta,r,CFL)
c     =====================================================

      implicit double precision (a-h,o-z)

c     # CFL-Superbee Limiter:
      s1 = theta*2.d0/dmax1(0.001d0,CFL)
      s2 = (1.d0 + CFL) / 3.d0
      phimax = theta * 2.d0 / (1.d0 - dmin1(0.999d0,CFL))
      ultra = dmax1(0.d0,dmin1(s1*r,phimax))
      c1 = 1.d0 + (s2 - beta/2.d0) * (r-1.d0)
      c2 = 1.d0 + (s2 + beta/2.d0) * (r-1.d0)
      cflsupbee = dmax1(0.d0, dmin1(ultra,dmax1(c1,c2)))

      return
      end
c
c     =====================================================
      double precision function suppow(r,CFL)
c     =====================================================

      implicit double precision (a-h,o-z)

c     # Superpower Limiter:
      s2 = (1.d0 + CFL) / 3.d0
      s3 = 1.d0 - s2
      if (r.LE.1.d0) then
         pp = dabs(2.d0/CFL)*2.d0*s3
      else
         pp = dabs(2.d0/(1.d0-CFL))*2.d0*s2
      endif
      rabs = dabs(r)
      rfrac = dabs((1.d0-rabs)/(1.d0+rabs))
      signum = dmax1(0.d0,dsign(1.d0,r))
      suppow = signum * (s3+s2*r) * (1.d0-rfrac**pp)

      return
      end

