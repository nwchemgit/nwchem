      integer function ishft(i,j)
*
* $Id: linux_shift.F,v 1.3 1997-10-31 20:45:33 d3e129 Exp $
*
      implicit none
      integer i,j
c
c     routine to return i shifted by j bits
c
      if(j.eq.0) ishft=i
      if(j.gt.0) ishft=lshift(i,j)
      if(j.lt.0) ishft=rshift(i,-j)
      return
      end
