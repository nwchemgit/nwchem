#if ( __GNUC__ >= 4) || defined(PGLINUX)
#else
      integer function ishft(i,j)
*
* $Id$
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
