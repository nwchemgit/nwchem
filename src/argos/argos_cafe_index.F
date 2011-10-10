      subroutine argos_cafe_index(ip,np)
c
      implicit none
c
      integer np
      integer ip(0:np,2)
c
      integer i,j
c
      ip(0,1)=0
      do 1 i=1,np
      j=iabs(ip(i,1))
      if(j.gt.0) ip(j,1)=i
    1 continue
      do 2 i=1,np
      if(ip(i,1).lt.0) ip(i,1)=ip(i-1,1)
    2 continue
c
      return
      end
c $Id$
