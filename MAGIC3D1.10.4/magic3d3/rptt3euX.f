c
c     ==================================================================
      subroutine rptt3  (ixyz,icoor,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,maux,ilr,impt,
     &                  bsasdq,cmbsasdq,cpbsasdq)
c     ==================================================================
c
c     # Riemann solver in the transverse direction for the
c     # Euler equations.
c     #
c     # Uses Roe averages and other quantities which were
c     # computed in rpn3eu and stored in the common block comroe.
c     #
c     #
c     # On input,
c
c     #    ql,qr is the data along some one-dimensional slice, as in rpn3
c     #         This slice is
c     #             in the x-direction if ixyz=1,
c     #             in the y-direction if ixyz=2, or
c     #             in the z-direction if ixyz=3.
c
c     #    bsasdq is an array of flux differences that result from a
c     #         transverse splitting (a previous call to rpt3).
c     #         This stands for B^* A^* \Dq but could represent any of
c     #         6 possibilities, e.g.  C^* B^* \Dq, as specified by ixyz
c     #         and icoor (see below).
c     #         Moreover, each * represents either + or -, as specified by
c     #         imp and impt.
c
c     #    ixyz indicates the direction of the original Riemann solve,
c     #         called the x-like direction in the table below:
c
c     #               x-like direction   y-like direction   z-like direction
c     #      ixyz=1:        x                  y                  z
c     #      ixyz=2:        y                  z                  x
c     #      ixyz=3:        z                  x                  y
c
c     #    icoor indicates direction in which the transverse solve should
c     #         be performed.
c     #      icoor=2: split in the y-like direction.
c     #      icoor=3: split in the z-like direction.
c     #      icoor=3: split in the z-like direction.
c
c     #    For example,
c     #        ixyz=1, icoor=3 means bsasdq=B^*A^*\Dq, and should be
c     #                        split in z into
c     #                           cmbsasdq = C^-B^*A^*\Dq,
c     #                           cpbsasdq = C^+B^*A^*\Dq.
c     #
c     #        ixyz=2, icoor=3 means bsasdq=C^*B^*\Dq, and should be
c     #                        split in x into
c     #                           cmbsasdq = A^-C^*B^*\Dq,
c     #                           cpbsasdq = A^+C^*B^*\Dq.
c
c
      implicit real*8(a-h,o-z)
      dimension       ql(1-mbc:maxm+mbc, meqn)
      dimension       qr(1-mbc:maxm+mbc, meqn)
      dimension   bsasdq(1-mbc:maxm+mbc, meqn)
      dimension cmbsasdq(1-mbc:maxm+mbc, meqn)
      dimension cpbsasdq(1-mbc:maxm+mbc, meqn)
      dimension     aux1(1-mbc:maxm+mbc, maux, 3)
      dimension     aux2(1-mbc:maxm+mbc, maux, 3)
      dimension     aux3(1-mbc:maxm+mbc, maux, 3)
c
      dimension waveb(meqn,3),sb(3)
      parameter (maxmrp = 2300)
      integer tflag

      common /comroe/ u2v2w2(-1:maxmrp),
     &       u(-1:maxmrp),v(-1:maxmrp),w(-1:maxmrp),enth(-1:maxmrp),
     &       a(-1:maxmrp),g1a2(-1:maxmrp),euv(-1:maxmrp)
c
      if (-3.gt.1-mbc .or. maxmrp .lt. maxm+mbc) then
	 write(6,*) 'need to increase maxmrp in rp3t'
	 stop
	 endif
c
      if(ixyz .eq. 1)then
	  mu = 2
	  mv = 3
          mw = 4
      else if(ixyz .eq. 2)then
	  mu = 3
	  mv = 4
          mw = 2
      else
          mu = 4
          mv = 2
          mw = 3
      endif
c
c     # Solve Riemann problem in the second coordinate direction
c
      if( icoor .eq. 2 )then

	 do 20 i = 2-mbc, mx+mbc
            a4 = g1a2(i) * (euv(i)*bsasdq(i,1)
     &             + u(i)*bsasdq(i,mu) + v(i)*bsasdq(i,mv)
     &             + w(i)*bsasdq(i,mw) - bsasdq(i,5))
	    a2 = bsasdq(i,mu) - u(i)*bsasdq(i,1)
            a3 = bsasdq(i,mw) - w(i)*bsasdq(i,1)
	    a5 = (bsasdq(i,mv) + (a(i)-v(i))*bsasdq(i,1) - a(i)*a4)
     &              / (2.d0*a(i))
	    a1 = bsasdq(i,1) - a4 - a5
c
            waveb(1,1)  = a1
            waveb(mu,1) = a1*u(i)
            waveb(mv,1) = a1*(v(i)-a(i))
            waveb(mw,1) = a1*w(i)
            waveb(5,1)  = a1*(enth(i) - v(i)*a(i))
	    sb(1) = v(i) - a(i)
c
            waveb(1,2)  = a4
            waveb(mu,2) = a2 + u(i)*a4
            waveb(mv,2) = v(i)*a4
            waveb(mw,2) = a3 + w(i)*a4
            waveb(5,2)  = a4*0.5d0*u2v2w2(i) + a2*u(i) + a3*w(i)
	    sb(2) = v(i)
c
            waveb(1,3)  = a5
            waveb(mu,3) = a5*u(i)
            waveb(mv,3) = a5*(v(i)+a(i))
            waveb(mw,3) = a5*w(i)
            waveb(5,3)  = a5*(enth(i)+v(i)*a(i))
	    sb(3) = v(i) + a(i)
c
            do m = 6, meqn
               waveb(m,1) = 0.d0
               waveb(m,2) = bsasdq(i,m)
               waveb(m,3) = 0.d0
            enddo 
c
	 do 25 m = 1, meqn
	    cmbsasdq(i,m) = 0.d0
	    cpbsasdq(i,m) = 0.d0
	    do 25 mws = 1, 3
	       cmbsasdq(i,m) = cmbsasdq(i,m)
     &			   + dmin1(sb(mws), 0.d0) * waveb(m,mws)
	       cpbsasdq(i,m) = cpbsasdq(i,m)
     &			   + dmax1(sb(mws), 0.d0) * waveb(m,mws)
 25         continue
c
 20      continue
c
      else
c
c        # Solve Riemann problem in the third coordinate direction
c
	 do 30 i = 2-mbc, mx+mbc
            a4 = g1a2(i) * (euv(i)*bsasdq(i,1)
     &             + u(i)*bsasdq(i,mu) + v(i)*bsasdq(i,mv)
     &             + w(i)*bsasdq(i,mw) - bsasdq(i,5))
	    a2 = bsasdq(i,mu) - u(i)*bsasdq(i,1)
            a3 = bsasdq(i,mv) - v(i)*bsasdq(i,1)
	    a5 = (bsasdq(i,mw) + (a(i)-w(i))*bsasdq(i,1) - a(i)*a4)
     &              / (2.d0*a(i))
	    a1 = bsasdq(i,1) - a4 - a5
c
            waveb(1,1)  = a1
            waveb(mu,1) = a1*u(i)
            waveb(mv,1) = a1*v(i)
            waveb(mw,1) = a1*(w(i) - a(i))
            waveb(5,1)  = a1*(enth(i) - w(i)*a(i))
	    sb(1) = w(i) - a(i)
c
            waveb(1,2)  = a4
            waveb(mu,2) = a2 + u(i)*a4
            waveb(mv,2) = a3 + v(i)*a4
            waveb(mw,2) = w(i)*a4
            waveb(5,2)  = a4*0.5d0*u2v2w2(i) + a2*u(i) + a3*v(i)
	    sb(2) = w(i)
c
            waveb(1,3)  = a5
            waveb(mu,3) = a5*u(i)
            waveb(mv,3) = a5*v(i)
            waveb(mw,3) = a5*(w(i)+a(i))
            waveb(5,3)  = a5*(enth(i)+w(i)*a(i))
	    sb(3) = w(i) + a(i)
c
            do m = 6, meqn
               waveb(m,1) = 0.d0
               waveb(m,2) = bsasdq(i,m)
               waveb(m,3) = 0.d0
            enddo
c
	 do 35 m = 1, meqn
	    cmbsasdq(i,m) = 0.d0
	    cpbsasdq(i,m) = 0.d0
	    do 35 mws = 1, 3
	       cmbsasdq(i,m) = cmbsasdq(i,m)
     &			   + dmin1(sb(mws), 0.d0) * waveb(m,mws)
	       cpbsasdq(i,m) = cpbsasdq(i,m)
     &			   + dmax1(sb(mws), 0.d0) * waveb(m,mws)
 35         continue
c
 30      continue
c
      endif
c
      return
      end
