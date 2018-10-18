c
c
c
c     ==================================================================
      subroutine dimsp3(maxm,maxmx,maxmy,maxmz,meqn,mwaves,mbc,mx,my,
     &                  mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cfl,
     &                  cflv,qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,
     &                  dtdz1d,aux1,aux2,aux3,maux,work,mwork,
     &                  rpn3,rpt3,rptt3)
c     ==================================================================
c
c     # Take one time step, updating q, using dimensional
c     # splitting. Only first order Godunov splitting is
c     # included. The main reason why not also Strang splitting 
c     # is included, is that it is not clear how to include
c     # boundary conditions in this case.
c
      implicit real*8(a-h,o-z)
      external rpn3,rpt3,rptt3
      dimension qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 
     &          1-mbc:maxmz+mbc, meqn)
      dimension qnew(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 
     &          1-mbc:maxmz+mbc, meqn)
      dimension  q1d(1-mbc:maxm+mbc, meqn)
      dimension cflv(4)
      dimension qadd(1-mbc:maxm+mbc, meqn)
      dimension fadd(1-mbc:maxm+mbc, meqn)
      dimension gadd(1-mbc:maxm+mbc, meqn, 2, -1:1)
      dimension hadd(1-mbc:maxm+mbc, meqn, 2, -1:1)
      dimension aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 
     &              1-mbc:maxmz+mbc, *)
      dimension aux1(1-mbc:maxm+mbc, maux, 3)
      dimension aux2(1-mbc:maxm+mbc, maux, 3)
      dimension aux3(1-mbc:maxm+mbc, maux, 3)
      dimension dtdx1d(1-mbc:maxm+mbc)
      dimension dtdy1d(1-mbc:maxm+mbc)
      dimension dtdz1d(1-mbc:maxm+mbc)
      dimension method(7),mthlim(mwaves)
      dimension work(mwork)
c
c     # Take one full time step in the x-direction
c
      idir = 1
      call step3ds(maxm,maxmx,maxmy,maxmz,meqn,mwaves,mbc,mx,my,
     &             mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cflx,
     &             qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d,
     &             aux1,aux2,aux3,maux,work,mwork,
     &             rpn3,rpt3,rptt3,idir)
c
      if (cflx .gt. cflv(1)) then
c        # Abort if the Courant number was too large in x-sweep
         cfl = cflx
         return
         endif
c
c     # Take one full time step in the y-direction
c
      idir = 2
c
      call step3ds(maxm,maxmx,maxmy,maxmz,meqn,mwaves,mbc,mx,my,
     &             mz,qnew,qnew,aux,dx,dy,dz,dt,method,mthlim,cfly,
     &             qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d,
     &             aux1,aux2,aux3,maux,work,mwork,
     &             rpn3,rpt3,rptt3,idir)
c
      if (cfly .gt. cflv(1)) then
c        # Abort if the Courant number was too large in y-sweep
         cfl = cfly
         return
         endif
c
c     # Take one full time step in the z-direction
c
      idir = 3
      call step3ds(maxm,maxmx,maxmy,maxmz,meqn,mwaves,mbc,mx,my,
     &             mz,qnew,qnew,aux,dx,dy,dz,dt,method,mthlim,cflz,
     &             qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d,
     &             aux1,aux2,aux3,maux,work,mwork,
     &             rpn3,rpt3,rptt3,idir)
c
      cfl = dmax1(cflx,cfly,cflz)
c
      return
      end
