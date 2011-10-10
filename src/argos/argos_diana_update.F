      subroutine argos_diana_update(x,m,ifr,ito)
c
      implicit none
c
      integer m,ifr,ito
      real*8 x(m,3)
c
#include "global.fh"
#include "msgids.fh"
c
      integer i,j
c
      if(ifr.gt.1) then
      do 1 j=1,3
      do 2 i=1,ifr-1
      x(i,j)=0.0d0
    2 continue
    1 continue
      endif
      if(ito.lt.m) then
      do 3 j=1,3
      do 4 i=ito+1,m
      x(i,j)=0.0d0
    4 continue
    3 continue
      endif
c
      call ga_dgop(mag_d12,x,3*m,'+')
c
      return
      end
c $Id$
