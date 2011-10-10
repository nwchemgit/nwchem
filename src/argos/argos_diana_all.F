      subroutine argos_diana_all(isel,ival,wt,numsel)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "msgids.fh"
#include "rtdb.fh"
c
      integer isel(nsa)
      real*8 wt(nsa)
      integer i,ival,numsel
c
      do 1 i=1,nsa
      isel(i)=ival
    1 continue
c
      numsel=0
      if(ival.ne.0) numsel=nsa
      do 2 i=1,nsa
      wt(i)=1.0d0
      if(isel(i).eq.0) wt(i)=0.0d0
    2 continue
c
      return
      end
c $Id$
