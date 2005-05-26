#ifdef GCC4
#else
      integer function ishft(i,j)
*
* $Id: linux_shift.F,v 1.4 2005-05-26 00:23:57 edo Exp $
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
#endif
