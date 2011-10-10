      subroutine argos_diana_evalangl(xtrj,tval,itval)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "argos_diana_common.fh"
c
      real*8 argos_diana_angle
      external argos_diana_angle
c
      integer itval
      real*8 xtrj(ntrj,3,nfrdim),tval(ntval,nfrdim)
c
      real*8 wtrj(3)
c
      integer i,k
c
      do 1 i=1,nangl
      do 2 k=1,nfrme
      tval(itval+i-1,k)=argos_diana_angle(xtrj(1,1,k),ntrj,wtrj,1,3,
     + iangl(i,1,1),iangl(i,1,2),iangl(i,2,1),iangl(i,2,2),
     + iangl(i,3,1),iangl(i,3,2))
    2 continue
    1 continue
c
      return
      end
c $Id$
