c
c     ==================================================================
      subroutine step3(maxm,maxmx,maxmy,maxmz,meqn,mwaves,mbc,mx,my,
     &                 mz,qold,qnew,aux,dx,dy,dz,dt,method,mthlim,cfl,
     &                 qadd,fadd,gadd,hadd,q1d,dtdx1d,dtdy1d,dtdz1d,
     &                 aux1,aux2,aux3,maux,work,mwork,rpn3,rpt3,rptt3)
c     ==================================================================
c
c     # Take one time step, updating q.
c     # On entry, qold and qnew should be identical and give the
c     #    initial data for this step
c     # On exit, qnew returns values at the end of the time step.
c     #    qold is unchanged.
c
c     # qadd is used to return increments to q from flux3.
c     # fadd, gadd and hadd are used to return flux increments from flux3.
c     # See the flux3 documentation for more information.
c
c
      implicit real*8(a-h,o-z)
      external rpn3,rpt3,rptt3
      dimension qold(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &          1-mbc:maxmz+mbc, meqn)
      dimension qnew(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &          1-mbc:maxmz+mbc, meqn)
      dimension  q1d(1-mbc:maxm+mbc, meqn)
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


      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom

c
c
c     # partition work array into pieces needed for local storage in
c     # flux3 routine.  Find starting index of each piece:
c
      i0wave     = 1
      i0s        = i0wave     + (maxm+2*mbc)*meqn*mwaves
      i0amdq     = i0s        + (maxm+2*mbc)*mwaves
      i0apdq     = i0amdq     + (maxm+2*mbc)*meqn
      i0cqxx     = i0apdq     + (maxm+2*mbc)*meqn
      i0bmamdq   = i0cqxx     + (maxm+2*mbc)*meqn
      i0bmapdq   = i0bmamdq   + (maxm+2*mbc)*meqn
      i0bpamdq   = i0bmapdq   + (maxm+2*mbc)*meqn
      i0bpapdq   = i0bpamdq   + (maxm+2*mbc)*meqn
      i0cmamdq   = i0bpapdq   + (maxm+2*mbc)*meqn
      i0cmapdq   = i0cmamdq   + (maxm+2*mbc)*meqn
      i0cpamdq   = i0cmapdq   + (maxm+2*mbc)*meqn
      i0cpapdq   = i0cpamdq   + (maxm+2*mbc)*meqn
      i0cmamdq2  = i0cpapdq   + (maxm+2*mbc)*meqn
      i0cmapdq2  = i0cmamdq2  + (maxm+2*mbc)*meqn
      i0cpamdq2  = i0cmapdq2  + (maxm+2*mbc)*meqn
      i0cpapdq2  = i0cpamdq2  + (maxm+2*mbc)*meqn
      i0bmcqxxp  = i0cpapdq2  + (maxm+2*mbc)*meqn
      i0bmcqxxm  = i0bmcqxxp  + (maxm+2*mbc)*meqn
      i0bpcqxxp  = i0bmcqxxm   + (maxm+2*mbc)*meqn
      i0bpcqxxm  = i0bpcqxxp   + (maxm+2*mbc)*meqn
      i0cmcqxxp  = i0bpcqxxm   + (maxm+2*mbc)*meqn
      i0cmcqxxm  = i0cmcqxxp   + (maxm+2*mbc)*meqn
      i0cpcqxxp  = i0cmcqxxm   + (maxm+2*mbc)*meqn
      i0cpcqxxm  = i0cpcqxxp   + (maxm+2*mbc)*meqn
      i0bmcmamdq = i0cpcqxxm   + (maxm+2*mbc)*meqn
      i0bmcmapdq = i0bmcmamdq + (maxm+2*mbc)*meqn
      i0bpcmamdq = i0bmcmapdq + (maxm+2*mbc)*meqn
      i0bpcmapdq = i0bpcmamdq + (maxm+2*mbc)*meqn
      i0bmcpamdq = i0bpcmapdq + (maxm+2*mbc)*meqn
      i0bmcpapdq = i0bmcpamdq + (maxm+2*mbc)*meqn
      i0bpcpamdq = i0bmcpapdq + (maxm+2*mbc)*meqn
      i0bpcpapdq = i0bpcpamdq + (maxm+2*mbc)*meqn
      iused      = i0bpcpapdq + (maxm+2*mbc)*meqn - 1
c
      if (iused.gt.mwork) then
c        # This shouldn't happen due to checks in claw3
         write(6,*) '*** not enough work space in step3'
         write(6,*) '*** iused = ', iused, '   mwork =',mwork
         stop
      endif
c


      mcapa = method(6)
      maux = method(7)
      cfl = 0.d0
      dtdx = dt/dx
      dtdy = dt/dy
      dtdz = dt/dz
c
      if (mcapa.eq.0) then
c        # no capa array:
         do 5 i=1-mbc,maxm+mbc
            dtdx1d(i) = dtdx
            dtdy1d(i) = dtdy
            dtdz1d(i) = dtdz
    5       continue
         endif
c
c
c     # perform x-sweeps
c     ==================
c

      do 50 k = 0,mz+1
         do 50 j = 0,my+1

            do 20 m=1,meqn
               do 20 i = 1-mbc, mx+mbc
c                 # copy data along a slice into 1d array:
                  q1d(i,m) = qold(i,j,k,m)
   20          continue
c
         if (mcapa.gt.0)  then
           do 23 i = 1-mbc, mx+mbc
               dtdx1d(i) = dtdx / aux(i,j,k,mcapa)
   23      continue
         endif
c
         if (maux .gt. 0)  then
             do 22 ma=1,maux
               do 22 i = 1-mbc, mx+mbc
                 aux1(i,ma,1) = aux(i,j-1,k-1,ma)
                 aux1(i,ma,2) = aux(i,j-1,k,ma)
                 aux1(i,ma,3) = aux(i,j-1,k+1,ma)
                 aux2(i,ma,1) = aux(i,j,k-1,ma)
                 aux2(i,ma,2) = aux(i,j,k,ma)
                 aux2(i,ma,3) = aux(i,j,k+1,ma)
                 aux3(i,ma,1) = aux(i,j+1,k-1,ma)
                 aux3(i,ma,2) = aux(i,j+1,k,ma)
                 aux3(i,ma,3) = aux(i,j+1,k+1,ma)
   22          continue
           endif
c
c           # Store the value of j and k along this slice in the common block
c           # comxyt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
            jcom = j
            kcom = k
c
c           # compute modifications fadd, gadd and hadd to fluxes along
c           # this slice:
c
            call flux3(1,maxm,meqn,mwaves,mbc,mx,
     &                 q1d,dtdx1d,dtdy,dtdz,aux1,aux2,aux3,maux,
     &                 method,mthlim,qadd,fadd,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),
     &                 work(i0bmcqxxp),work(i0bpcqxxp),
     &                 work(i0bmcqxxm),work(i0bpcqxxm),
     &                 work(i0cmcqxxp),work(i0cpcqxxp),
     &                 work(i0cmcqxxm),work(i0cpcqxxm),
     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cfl = dmax1(cfl,cfl1d)
c
c           # update qnew by flux differencing.
c           # (rather than maintaining arrays f, g and h for the total fluxes,
c           # the modifications are used immediately to update qnew
c           # in order to save storage.)
c
            if(mcapa. eq. 0)then
c
            do 30 m=1,meqn
               do 30 i=1,mx
                  qnew(i,j,k,m) = qnew(i,j,k,m) + qadd(i,m)
     &                        - dtdx * (fadd(i+1,m) - fadd(i,m))
     &                        - dtdy * (gadd(i,m,2,0) - gadd(i,m,1,0))
     &                        - dtdz * (hadd(i,m,2,0) - hadd(i,m,1,0))
                  qnew(i,j-1,k,m)   = qnew(i,j-1,k,m)
     &                              - dtdy * gadd(i,m,1,0)
     &                              - dtdz * ( hadd(i,m,2,-1)
     &                                      -   hadd(i,m,1,-1) )
                  qnew(i,j-1,k-1,m) = qnew(i,j-1,k-1,m)
     &                              - dtdy * gadd(i,m,1,-1)
     &                              - dtdz * hadd(i,m,1,-1)
                  qnew(i,j,k-1,m)   = qnew(i,j,k-1,m)
     &                              - dtdy * ( gadd(i,m,2,-1)
     &                                     -   gadd(i,m,1,-1) )
     &                              - dtdz * hadd(i,m,1,0)
                  qnew(i,j+1,k-1,m) = qnew(i,j+1,k-1,m)
     &                              + dtdy * gadd(i,m,2,-1)
     &                              - dtdz * hadd(i,m,1,1)
                  qnew(i,j+1,k,m)   = qnew(i,j+1,k,m)
     &                              + dtdy * gadd(i,m,2,0)
     &                              - dtdz * ( hadd(i,m,2,1)
     &                                     -   hadd(i,m,1,1) )
                  qnew(i,j+1,k+1,m) = qnew(i,j+1,k+1,m)
     &                              + dtdy * gadd(i,m,2,1)
     &                              + dtdz * hadd(i,m,2,1)
                  qnew(i,j,k+1,m)   = qnew(i,j,k+1,m)
     &                              - dtdy * ( gadd(i,m,2,1)
     &                                     -   gadd(i,m,1,1) )
     &                              + dtdz * hadd(i,m,2,0)
                  qnew(i,j-1,k+1,m) = qnew(i,j-1,k+1,m)
     &                              - dtdy * gadd(i,m,1,1)
     &                              + dtdz * hadd(i,m,2,-1)
c
   30          continue
c
            else
c              # with capa array
               do 40 m=1,meqn
                  do 40 i=1,mx
                     qnew(i,j,k,m) = qnew(i,j,k,m) + qadd(i,m)
     &                        - (dtdx * (fadd(i+1,m) - fadd(i,m))
     &                        +  dtdy * (gadd(i,m,2,0) - gadd(i,m,1,0))
     &                        +  dtdz * (hadd(i,m,2,0) - hadd(i,m,1,0)))
     &                        / aux(i,j,k,mcapa)

                     qnew(i,j-1,k,m)   = qnew(i,j-1,k,m)
     &                                 - (dtdy * gadd(i,m,1,0)
     &                                 +  dtdz * ( hadd(i,m,2,-1)
     &                                        -   hadd(i,m,1,-1) ))
     &                                 / aux(i,j-1,k,mcapa)
                     qnew(i,j-1,k-1,m) = qnew(i,j-1,k-1,m)
     &                                 - (dtdy * gadd(i,m,1,-1)
     &                                 +  dtdz * hadd(i,m,1,-1))
     &                                 / aux(i,j-1,k-1,mcapa)
                     qnew(i,j,k-1,m)   = qnew(i,j,k-1,m)
     &                                 - (dtdy * ( gadd(i,m,2,-1)
     &                                        -   gadd(i,m,1,-1) )
     &                                 +  dtdz * hadd(i,m,1,0))
     &                                 / aux(i,j,k-1,mcapa)
                     qnew(i,j+1,k-1,m) = qnew(i,j+1,k-1,m)
     &                                 + (dtdy * gadd(i,m,2,-1)
     &                                 -  dtdz * hadd(i,m,1,1))
     &                                 / aux(i,j+1,k-1,mcapa)
                     qnew(i,j+1,k,m)   = qnew(i,j+1,k,m)
     &                                 + (dtdy * gadd(i,m,2,0)
     &                                 - dtdz * ( hadd(i,m,2,1)
     &                                        -   hadd(i,m,1,1) ))
     &                                 / aux(i,j+1,k,mcapa)
                     qnew(i,j+1,k+1,m) = qnew(i,j+1,k+1,m)
     &                                 + (dtdy * gadd(i,m,2,1)
     &                                 +  dtdz * hadd(i,m,2,1))
     &                                 / aux(i,j+1,k+1,mcapa)
                     qnew(i,j,k+1,m)   = qnew(i,j,k+1,m)
     &                                 - (dtdy * ( gadd(i,m,2,1)
     &                                        -   gadd(i,m,1,1) )
     &                                 -  dtdz * hadd(i,m,2,0))
     &                                 / aux(i,j,k+1,mcapa)
                     qnew(i,j-1,k+1,m) = qnew(i,j-1,k+1,m)
     &                                 - (dtdy * gadd(i,m,1,1)
     &                                 -  dtdz * hadd(i,m,2,-1))
     &                                 / aux(i,j-1,k+1,mcapa)
c
c
   40          continue
c
            endif
c
   50    continue
   51    continue

c
c
c     # perform y sweeps
c     ==================
c
c
      do 100 k = 0, mz+1
         do 100 i = 0, mx+1
c
            do 70 m=1,meqn
               do 70 j = 1-mbc, my+mbc
c                 # copy data along a slice into 1d array:
                  q1d(j,m) = qold(i,j,k,m)
   70          continue
c
         if (mcapa.gt.0)  then
           do 71 j = 1-mbc, my+mbc
               dtdy1d(j) = dtdy / aux(i,j,k,mcapa)
   71      continue
         endif
c
         if (maux .gt. 0)  then
             do 72 ma=1,maux
               do 72 j = 1-mbc, my+mbc
                 aux1(j,ma,1) = aux(i-1,j,k-1,ma)
                 aux1(j,ma,2) = aux(i,j,k-1,ma)
                 aux1(j,ma,3) = aux(i+1,j,k-1,ma)
                 aux2(j,ma,1) = aux(i-1,j,k,ma)
                 aux2(j,ma,2) = aux(i,j,k,ma)
                 aux2(j,ma,3) = aux(i+1,j,k,ma)
                 aux3(j,ma,1) = aux(i-1,j,k+1,ma)
                 aux3(j,ma,2) = aux(i,j,k+1,ma)
                 aux3(j,ma,3) = aux(i+1,j,k+1,ma)
   72          continue
         endif
c
c           # Store the value of i and k along this slice in the common block
c           # comxyzt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
            icom = i
            kcom = k
c
c           # compute modifications fadd, gadd and hadd to fluxes along this
c           # slice:
c
            call flux3(2,maxm,meqn,mwaves,mbc,my,
     &                 q1d,dtdy1d,dtdz,dtdx,aux1,aux2,aux3,maux,
     &                 method,mthlim,qadd,fadd,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),
     &                 work(i0bmcqxxp),work(i0bpcqxxp),
     &                 work(i0bmcqxxm),work(i0bpcqxxm),
     &                 work(i0cmcqxxp),work(i0cpcqxxp),
     &                 work(i0cmcqxxm),work(i0cpcqxxm),
     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cfl = dmax1(cfl,cfl1d)
c
c           # update qnew by flux differencing.
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the g-fluxes
c           # gadd - modifies the h-fluxes
c           # hadd - modifies the f-fluxes
c
            if( mcapa.eq. 0)then
c               # no capa array.  Standard flux differencing:
            do 80 m=1,meqn
               do 80 j=1,my
                  qnew(i,j,k,m) = qnew(i,j,k,m) + qadd(j,m)
     &                        - dtdy * (fadd(j+1,m) - fadd(j,m))
     &                        - dtdz * (gadd(j,m,2,0) - gadd(j,m,1,0))
     &                        - dtdx * (hadd(j,m,2,0) - hadd(j,m,1,0))
                  qnew(i,j,k+1,m)   = qnew(i,j,k+1,m)
     &                              + dtdz * gadd(j,m,2,0)
     &                              - dtdx * ( hadd(j,m,2,1)
     &                                     -   hadd(j,m,1,1) )
                  qnew(i+1,j,k+1,m) = qnew(i+1,j,k+1,m)
     &                              + dtdz * gadd(j,m,2,1)
     &                              +  dtdx * hadd(j,m,2,1)
                  qnew(i+1,j,k,m)   = qnew(i+1,j,k,m)
     &                              - dtdz * ( gadd(j,m,2,1)
     &                                     -   gadd(j,m,1,1) )
     &                              + dtdx * hadd(j,m,2,0)
                  qnew(i+1,j,k-1,m) = qnew(i+1,j,k-1,m)
     &                              - dtdz * gadd(j,m,1,1)
     &                              + dtdx * hadd(j,m,2,-1)
                  qnew(i,j,k-1,m)   = qnew(i,j,k-1,m)
     &                              - dtdz * gadd(j,m,1,0)
     &                              - dtdx * ( hadd(j,m,2,-1)
     &                                     -   hadd(j,m,1,-1) )
                  qnew(i-1,j,k-1,m) = qnew(i-1,j,k-1,m)
     &                              - dtdz * gadd(j,m,1,-1)
     &                              - dtdx * hadd(j,m,1,-1)
                  qnew(i-1,j,k,m)   = qnew(i-1,j,k,m)
     &                              - dtdz * ( gadd(j,m,2,-1)
     &                                     -   gadd(j,m,1,-1) )
     &                              - dtdx * hadd(j,m,1,0)
                  qnew(i-1,j,k+1,m) = qnew(i-1,j,k+1,m)
     &                              + dtdz * gadd(j,m,2,-1)
     &                              -  dtdx*hadd(j,m,1,1)
c
   80          continue
c
             else
c
c              #with capa array.
c
                do 85 m=1,meqn
                   do 85 j=1,my
                      qnew(i,j,k,m) = qnew(i,j,k,m) + qadd(j,m)
     &                        - (dtdy * (fadd(j+1,m) - fadd(j,m))
     &                        +  dtdz * (gadd(j,m,2,0) - gadd(j,m,1,0))
     &                        +  dtdx * (hadd(j,m,2,0) - hadd(j,m,1,0)))
     &                        / aux(i,j,k,mcapa)
                      qnew(i,j,k+1,m)   = qnew(i,j,k+1,m)
     &                                 + (dtdz * gadd(j,m,2,0)
     &                                 -  dtdx * ( hadd(j,m,2,1)
     &                                        -   hadd(j,m,1,1) ))
     &                                 / aux(i,j,k+1,mcapa)
                      qnew(i+1,j,k+1,m) = qnew(i+1,j,k+1,m)
     &                                  + (dtdz * gadd(j,m,2,1)
     &                                  +  dtdx * hadd(j,m,2,1))
     &                                  / aux(i+1,j,k+1,mcapa)
                      qnew(i+1,j,k,m)   = qnew(i+1,j,k,m)
     &                                  - (dtdz * ( gadd(j,m,2,1)
     &                                          -   gadd(j,m,1,1) )
     &                                  -  dtdx * hadd(j,m,2,0) )
     &                                  / aux(i+1,j,k,mcapa)
                      qnew(i+1,j,k-1,m) = qnew(i+1,j,k-1,m)
     &                                  - (dtdz * gadd(j,m,1,1)
     &                                  -  dtdx * hadd(j,m,2,-1))
     &                                  / aux(i+1,j,k-1,mcapa)
                      qnew(i,j,k-1,m)   = qnew(i,j,k-1,m)
     &                                  - (dtdz * gadd(j,m,1,0)
     &                                  +  dtdx * ( hadd(j,m,2,-1)
     &                                         -   hadd(j,m,1,-1) ))
     &                                  / aux(i,j,k-1,mcapa)
                      qnew(i-1,j,k-1,m) = qnew(i-1,j,k-1,m)
     &                                  - (dtdz * gadd(j,m,1,-1)
     &                                  +  dtdx * hadd(j,m,1,-1))
     &                                  / aux(i-1,j,k-1,mcapa)
                      qnew(i-1,j,k,m)   = qnew(i-1,j,k,m)
     &                                  - (dtdz * ( gadd(j,m,2,-1)
     &                                         -   gadd(j,m,1,-1) )
     &                                  +  dtdx * hadd(j,m,1,0))
     &                                  / aux(i-1,j,k,mcapa)
                      qnew(i-1,j,k+1,m) = qnew(i-1,j,k+1,m)
     &                                  + (dtdz * gadd(j,m,2,-1)
     &                                  -  dtdx*hadd(j,m,1,1))
     &                                  / aux(i-1,j,k+1,mcapa)
c

   85              continue
c
            endif
c
  100    continue
  101    continue
c
c
c
c     # perform z sweeps
c     ==================
c
c
      do 150 j = 0, my+1
         do 150 i = 0, mx+1
c
            do 110 m=1,meqn
               do 110 k = 1-mbc, mz+mbc
c                 # copy data along a slice into 1d array:
                  q1d(k,m) = qold(i,j,k,m)
 110           continue
c
         if (mcapa.gt.0)  then
           do 130 k = 1-mbc, mz+mbc
               dtdz1d(k) = dtdz / aux(i,j,k,mcapa)
 130       continue
         endif
c
         if (maux .gt. 0)  then
             do 131 ma=1,maux
               do 131 k = 1-mbc, mz+mbc
                 aux1(k,ma,1) = aux(i-1,j-1,k,ma)
                 aux1(k,ma,2) = aux(i-1,j,k,ma)
                 aux1(k,ma,3) = aux(i-1,j+1,k,ma)
                 aux2(k,ma,1) = aux(i,j-1,k,ma)
                 aux2(k,ma,2) = aux(i,j,k,ma)
                 aux2(k,ma,3) = aux(i,j+1,k,ma)
                 aux3(k,ma,1) = aux(i+1,j-1,k,ma)
                 aux3(k,ma,2) = aux(i+1,j,k,ma)
                 aux3(k,ma,3) = aux(i+1,j+1,k,ma)
  131          continue
           endif
c
c           # Store the value of i and j along this slice in the common block
c           # comxyzt in case it is needed in the Riemann solver (for
c           # variable coefficient problems)
c
            icom = i
            jcom = j
c
c           # compute modifications fadd, gadd and hadd to fluxes along this
c           # slice:
c
            call flux3(3,maxm,meqn,mwaves,mbc,mz,
     &                 q1d,dtdz1d,dtdx,dtdy,aux1,aux2,aux3,maux,
     &                 method,mthlim,qadd,fadd,gadd,hadd,cfl1d,
     &                 work(i0wave),work(i0s),work(i0amdq),
     &                 work(i0apdq),work(i0cqxx),
     &                 work(i0bmamdq),work(i0bmapdq),
     &                 work(i0bpamdq),work(i0bpapdq),
     &                 work(i0cmamdq),work(i0cmapdq),
     &                 work(i0cpamdq),work(i0cpapdq),
     &                 work(i0cmamdq2),work(i0cmapdq2),
     &                 work(i0cpamdq2),work(i0cpapdq2),
     &                 work(i0bmcqxxp),work(i0bpcqxxp),
     &                 work(i0bmcqxxm),work(i0bpcqxxm),
     &                 work(i0cmcqxxp),work(i0cpcqxxp),
     &                 work(i0cmcqxxm),work(i0cpcqxxm),
     &                 work(i0bmcmamdq),work(i0bmcmapdq),
     &                 work(i0bpcmamdq),work(i0bpcmapdq),
     &                 work(i0bmcpamdq),work(i0bmcpapdq),
     &                 work(i0bpcpamdq),work(i0bpcpapdq),
     &                 rpn3,rpt3,rptt3)
c
            cfl = dmax1(cfl,cfl1d)
c
c           # update qnew by flux differencing.
c           # Note that the roles of the flux updates are changed.
c           # fadd - modifies the h-fluxes
c           # gadd - modifies the f-fluxes
c           # hadd - modifies the g-fluxes
c
            if(mcapa .eq. 0)then
c
c              #no capa array. Standard flux differencing:
c
            do 120 m=1,meqn
               do 120 k=1,mz
                  qnew(i,j,k,m) = qnew(i,j,k,m) + qadd(k,m)
     &                        - dtdz * (fadd(k+1,m) - fadd(k,m))
     &                        - dtdx * (gadd(k,m,2,0) - gadd(k,m,1,0))
     &                        - dtdy * (hadd(k,m,2,0) - hadd(k,m,1,0))
                  qnew(i,j+1,k,m)   = qnew(i,j+1,k,m)
     &                              - dtdx * ( gadd(k,m,2,1)
     &                                     -   gadd(k,m,1,1) )
     &                              + dtdy * hadd(k,m,2,0)
                  qnew(i+1,j+1,k,m) = qnew(i+1,j+1,k,m)
     &                              + dtdx * gadd(k,m,2,1)
     &                              + dtdy * hadd(k,m,2,1)
                  qnew(i+1,j,k,m)   = qnew(i+1,j,k,m)
     &                              + dtdx * gadd(k,m,2,0)
     &                              - dtdy * ( hadd(k,m,2,1)
     &                                     -   hadd(k,m,1,1) )
                  qnew(i+1,j-1,k,m) = qnew(i+1,j-1,k,m)
     &                              + dtdx * gadd(k,m,2,-1)
     &                              - dtdy * hadd(k,m,1,1)
                  qnew(i,j-1,k,m)   = qnew(i,j-1,k,m)
     &                              - dtdx * ( gadd(k,m,2,-1)
     &                                     -   gadd(k,m,1,-1) )
     &                              - dtdy * hadd(k,m,1,0)
                  qnew(i-1,j-1,k,m) = qnew(i-1,j-1,k,m)
     &                              - dtdx * gadd(k,m,1,-1)
     &                              - dtdy * hadd(k,m,1,-1)
                  qnew(i-1,j,k,m)   = qnew(i-1,j,k,m)
     &                              - dtdx * gadd(k,m,1,0)
     &                              - dtdy * ( hadd(k,m,2,-1)
     &                                     -   hadd(k,m,1,-1) )
                  qnew(i-1,j+1,k,m) = qnew(i-1,j+1,k,m)
     &                              - dtdx * gadd(k,m,1,1)
     &                              + dtdy * hadd(k,m,2,-1)
c
 120           continue
c
            else
c
c              # with capa array
c
               do 145 m=1,meqn
                  do 145 k=1,mz
                     qnew(i,j,k,m) = qnew(i,j,k,m) + qadd(k,m)
     &                        - (dtdz * (fadd(k+1,m) - fadd(k,m))
     &                        +  dtdx * (gadd(k,m,2,0) - gadd(k,m,1,0))
     &                        +  dtdy * (hadd(k,m,2,0) - hadd(k,m,1,0)))
     &                        / aux(i,j,k,mcapa)
                     qnew(i,j+1,k,m) = qnew(i,j+1,k,m)
     &                               - (dtdx * ( gadd(k,m,2,1)
     &                                      -   gadd(k,m,1,1) )
     &                               -  dtdy * hadd(k,m,2,0))
     &                               / aux(i,j+1,k,mcapa)
                     qnew(i+1,j+1,k,m) = qnew(i+1,j+1,k,m)
     &                                 + (dtdx * gadd(k,m,2,1)
     &                                 +  dtdy * hadd(k,m,2,1))
     &                                 / aux(i+1,j+1,k,mcapa)
                     qnew(i+1,j,k,m)   = qnew(i+1,j,k,m)
     &                                 + (dtdx * gadd(k,m,2,0)
     &                                 - dtdy * ( hadd(k,m,2,1)
     &                                        -   hadd(k,m,1,1) ))
     &                                 / aux(i+1,j,k,mcapa)
                     qnew(i+1,j-1,k,m) = qnew(i+1,j-1,k,m)
     &                                 + (dtdx * gadd(k,m,2,-1)
     &                                 -  dtdy * hadd(k,m,1,1))
     &                                 / aux(i+1,j-1,k,mcapa)
                     qnew(i,j-1,k,m)   = qnew(i,j-1,k,m)
     &                                 - (dtdx * ( gadd(k,m,2,-1)
     &                                        -   gadd(k,m,1,-1) )
     &                                 +  dtdy * hadd(k,m,1,0))
     &                                 / aux(i,j-1,k,mcapa)
                     qnew(i-1,j-1,k,m) = qnew(i-1,j-1,k,m)
     &                                 - (dtdx * gadd(k,m,1,-1)
     &                                 +  dtdy * hadd(k,m,1,-1))
     &                                 / aux(i-1,j-1,k,mcapa)
                     qnew(i-1,j,k,m)   = qnew(i-1,j,k,m)
     &                                 - (dtdx * gadd(k,m,1,0)
     &                                 + dtdy * ( hadd(k,m,2,-1)
     &                                       -   hadd(k,m,1,-1) ))
     &                                 / aux(i-1,j,k,mcapa)
                     qnew(i-1,j+1,k,m) = qnew(i-1,j+1,k,m)
     &                                 - (dtdx * gadd(k,m,1,1)
     &                                 -  dtdy * hadd(k,m,2,-1))
     &                                 / aux(i-1,j+1,k,mcapa)
c
 145              continue
c
         endif
c
 150  continue
  151 continue
c
c
      return
      end
