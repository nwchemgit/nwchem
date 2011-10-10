      subroutine argos_cafe_binit(r,m1,m2)
c
      implicit none
c
      integer m1,m2
      real*8 r(m1,m2)
c
      integer i,j
c
      do 1 j=1,m2
      do 2 i=1,m1
      r(i,j)=0.0d0
    2 continue
    1 continue
c
      return
      end
c $Id$
