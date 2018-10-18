c
c
c     ==================================================================
      subroutine copyq3(maxmx,maxmy,maxmz,meqn,mbc,mx,my,mz,q1,q2)
c     ==================================================================
c
c     # copy the contents of q1 into q2
c
      implicit real*8(a-h,o-z)
      dimension q1(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 
     &             1-mbc:maxmz+mbc,meqn)
      dimension q2(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, 
     &             1-mbc:maxmz+mbc,meqn)
c
c
      do 10 m=1,meqn
         do 10 k = 1-mbc, mz+mbc
            do 10 j = 1-mbc, my+mbc
               do 10 i = 1-mbc, mx+mbc
                  q2(i,j,k,m) = q1(i,j,k,m)
   10       continue
c
      return
      end
