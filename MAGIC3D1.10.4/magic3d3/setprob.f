c
c
c     ============================================
      subroutine setprob
c     ============================================
c
      implicit double precision (a-h,o-z)
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
c
      open(unit=7,file='setprob.data',status='old',form='formatted')
      rewind(7)
c
      read(7,*) loadprof      
      read(7,*) nspecadv
      read(7,*) gamma     
      read(7,*) grav
      read(7,*) dens
      read(7,*) pres
      read(7,*) dmolec0
      read(7,*) dtherm0
      read(7,*) forcemth
      read(7,*) omega
      read(7,*) amplitude
      read(7,*) propx
      read(7,*) propz
      read(7,*) xpos
      read(7,*) ypos
      read(7,*) zpos
      read(7,*) xwidth
      read(7,*) ywidth
      read(7,*) zwidth
      read(7,*) tcenter
      read(7,*) twidth
      read(7,*) vsrcx
      read(7,*) amp0
      read(7,*) amp1
      read(7,*) omega1
      read(7,*) prop1z
      read(7,*) phi
      read(7,*) z1width
      read(7,*) z1pos
      read(7,*) z1Gpos
      read(7,*) windscale1
c           
      if (loadprof .EQ. 0) then 
         gamma1 = gamma - 1.d0
         speed = dsqrt(pres*gamma/dens)
         w_n = dsqrt(((grav**2)/(speed**2))*(gamma1))
         w_o = dsqrt((gamma**2)*(grav**2)/(4*(speed**2)))
      endif
c
      close(7)
c
      return
      end
