      subroutine argos_diana_evaltors(xtrj,tval,itval)
c
      implicit none
c
#include "msgids.fh"
#include "global.fh"
#include "argos_diana_common.fh"
c
      real*8 argos_diana_torsion
      external argos_diana_torsion
c
      integer itval
      real*8 xtrj(ntrj,3,nfrdim),tval(ntval,nfrdim)
c
      real*8 wtrj(3)
c
      integer i,k
c
      do 1 i=1,ntors
      do 2 k=1,nfrme
      tval(itval+i-1,k)=argos_diana_torsion(xtrj(1,1,k),ntrj,wtrj,1,3,
     + itors(i,1,1),itors(i,1,2),itors(i,2,1),itors(i,2,2),
     + itors(i,3,1),itors(i,3,2),itors(i,4,1),itors(i,4,2))
    2 continue
    1 continue
c
      return
      end
c $Id$
