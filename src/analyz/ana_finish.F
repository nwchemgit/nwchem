      subroutine ana_finish(sgmnam)
c
c $Id: ana_finish.F,v 1.1 1998-12-29 18:01:02 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      character*16 sgmnam(nsa)
c
      call ana_report(sgmnam,dbl_mb(i_xrms))
c
      if(.not.ma_pop_stack(l_wt))
     + call errquit('Failed to deallocate wt',0)
      if(.not.ma_pop_stack(l_xrms))
     + call errquit('Failed to deallocate xrms',0)
      if(.not.ma_pop_stack(l_xref))
     + call errquit('Failed to deallocate xref',0)
c
      if(me.eq.0) write(*,'(a)') 'Analysis arrays deallocated'
c
      return
      end
      subroutine ana_report(sgmnam,xrms)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      real*8 xrms(nsa)
c
      integer i,j
c
c     report results
c
      if(me.eq.0) write(lfnana,1000)
 1000 format('Analysis')
      if(np.gt.1) then
      call ga_dgop(mag_d03,xrms,nsa,'+')
      endif
      if(me.eq.0) then
      write(lfnana,1001) (sgmnam(i),i,
     + sqrt(xrms(i))/dble(ndata),i=1,nsa)
 1001 format(a16,i5,f12.6)
      close(unit=lfnana)
      endif
      return
      end
