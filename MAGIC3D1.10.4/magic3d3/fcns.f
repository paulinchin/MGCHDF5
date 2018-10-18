c   
c     
c     =====================================================
      double precision function gravity(z)
c     =====================================================
c
      implicit double precision (a-h,o-z)
      
      gravity=-3.989e14/((6.38d6+z)**2)
            
      return
      end     
c 
c     =====================================================
      double precision function pressure(z)
c     =====================================================
c      
      implicit double precision (a-h,o-z)
      common /param/ gamma, gamma1, grav, dens, pres, 
     &                  speed, dtherm0, dmolec0
          
      pressure=pres*dexp(-z*gamma*grav/(speed**2))

      return
      end
c
c     =====================================================
      double precision function density(z)
c     =====================================================
c    
      implicit double precision (a-h,o-z)  
      common /param/ gamma, gamma1, grav, dens, pres, 
     &                  speed, dtherm0, dmolec0
               
      density=dens*dexp(-z*gamma*grav/(speed**2))
      
      return
      end      
c      
c     =====================================================
      double precision function gaussian(x,y,z,t)
c     =====================================================
c     
      implicit double precision (a-h,o-z)
      common /gauss/ xpos, ypos, zpos, xwidth, ywidth, 
     &               zwidth, tcenter, twidth
      gaussian=dexp(-.5*((x-xpos)**2/xwidth**2))*
     &         dexp(-.5*((y-ypos)**2/ywidth**2))*
     &         dexp(-.5*((z-zpos)**2/zwidth**2))*
     &         dexp(-.5*((t-tcenter)**2/twidth**2))
           
      return
      end  
c      
c     =====================================================
      double precision function forcefnx(x,y,z,t)
c     =====================================================
c     
      implicit double precision (a-h,o-z)
      common /param/ gamma, gamma1, grav, dens, pres, 
     &                  speed, dtherm0, dmolec0
      common /forcing/ omega, amplitude, propx, propz, 
     &                  w_n, w_o, forcemth
      common /gauss/ xpos, ypos, zpos, xwidth, ywidth,
     &               zwidth, tcenter, twidth
      
      forcefnx=cos(omega*(t-tcenter))*
     &         cos(propx*(x-xpos)-propz*(z-zpos))
             
      return
      end        
      
c
c     =====================================================
      double precision function forcefnz(x,y,z,t)
c     =====================================================
c      
      implicit double precision (a-h,o-z)
      common /param/ gamma, gamma1, grav, dens, pres, 
     &                  speed, dtherm0, dmolec0
      common /forcing/ omega, amplitude, propx, propz,
     &                  w_n, w_o, forcemth      
      common /gauss/ xpos, ypos, zpos, xwidth, ywidth,
     &               zwidth, tcenter, twidth

      forcefnz=cos(omega*(t-tcenter)-
     &             propx*(x-xpos)-propz*(y-ypos))

      return
      end  
c
c     =====================================================
      double precision function windfnx(x,y,z,t,scale)
c     =====================================================
c      
      implicit double precision (a-h,o-z)
      common /wind/ omega1, amp0, amp1, prop1z, z1width, z1pos,
     &                  z1Gpos, phi, windscale1
     
      windfnx = amp1
      windfnx = windfnx*cos(prop1z*(z-z1pos)-omega1*t-phi)
      windfnx = windfnx*dexp(-.5*(((z-z1Gpos)**2)/z1width**2))
      if (windscale1 .EQ. 1) then
         windfnx = windfnx*scale
      endif
      windfnx = windfnx + amp0
      
      return
      end
