      subroutine argos_diana_wselect(isel,xs,iwsel,xw,range,ndxw)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
#include "msgids.fh"
c
      integer isel(msa),iwsel(mwm),ndxw(mwm)
      real*8 xs(msa,3),xw(mwm,mwa,3)
c
      integer i,j,k
      real*8 range,range2
c
      range2=range*range
c
c      if(me.eq.0) then
      do 1 i=1,nwm
      if(iwsel(i).eq.0) then
      do 2 j=1,nsa
      if(isel(j).ne.0) then
      do 3 k=1,nwa
      if((xw(i,k,1)-xs(j,1))**2+(xw(i,k,2)-xs(j,2))**2+
     + (xw(i,k,3)-xs(j,3))**2.lt.range2) then
      iwsel(i)=1
      goto 1
      endif
    3 continue
      endif
    2 continue
      endif
    1 continue
      nwsel=0
      do 4 i=1,nwm
      if(iwsel(i).gt.0) then
      nwsel=nwsel+1
      ndxw(nwsel)=i
      endif
    4 continue
c
      if(np.gt.1) then
      call ga_igop(mag_d14,nwsel,1,'+')
      endif
c      endif
c
      return
      end
c $Id$
