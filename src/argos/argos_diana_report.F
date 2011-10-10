      subroutine argos_diana_report(sgmnam,xrms,isel)
c
      implicit none
c
#include "argos_diana_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      integer isel(nsa)
      real*8 xrms(nsa)
c
      integer i,n
      real*8 arms,bfact
c
c     report results
c
      bfact=128.0d0*atan(1.0d0)*atan(1.0d0)/3.0d0
c
      if(me.eq.0) then
      write(lfnrms,1000)
 1000 format('analysis')
      endif
      if(np.gt.1) then
      call ga_dgop(mag_d12,xrms,nsa,'+')
      endif
      if(me.eq.0) then
      write(lfnrms,1001) (sgmnam(i),i,isel(i),
     + sqrt(xrms(i)/dble(ndata)),bfact*xrms(i)/dble(ndata),i=1,nsa)
 1001 format(a16,i10,i5,2f12.6)
      arms=0.0d0
      n=0
      if(me.eq.0) write(lfnrms,1000)
      do 1 i=1,nsa
      if(isel(i).ne.0) then
      arms=arms+xrms(i)
      n=n+1
      endif
      if(i.eq.nsa) goto 2
      if(sgmnam(i)(11:16).ne.sgmnam(i+1)(11:16)) goto 2
      goto 1
    2 continue
      if(ndata*n.gt.0) then
      write(lfnrms,1002) sgmnam(i)(1:5),sgmnam(i)(11:16),
     + sqrt(arms/dble(ndata*n)),bfact*arms/dble(ndata*n)
 1002 format(a5,a6,2f12.6)
      endif
      arms=0.0d0
      n=0
    1 continue
      close(unit=lfnrms)
      endif
      return
      end
c $Id$
