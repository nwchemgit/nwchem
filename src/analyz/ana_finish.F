      subroutine ana_finish(sgmnam,ltask)
c
c $Id: ana_finish.F,v 1.14 2000-07-19 18:26:37 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "mafdecls.fh"
c
      character*16 sgmnam(nsa)
      logical ltask
c
      call ana_report(sgmnam,dbl_mb(i_xrms),int_mb(i_isel))
c
c
      if(me.eq.0) then
      if(lcop.gt.0) then
      if(fmtcop.eq.'eci') write(lfncop,'(a)') 'End'
      close(unit=lfncop)
      write(*,'(a)') ' Closing copy file '
      endif
      if(lsup.gt.0)  then
      if(fmtsup.eq.'eci') write(lfnsup,'(a)') 'End'
      close(unit=lfnsup)
      write(*,'(a)') ' Closing super file '
      endif
      endif
c
c     deallocate memory
c
      if(ltask) then
      if(.not.ma_pop_stack(l_wsel))
     + call errquit('Failed to deallocate wsel',0)
      if(.not.ma_pop_stack(l_wdat))
     + call errquit('Failed to deallocate wdat',0)
      if(.not.ma_pop_stack(l_idat))
     + call errquit('Failed to deallocate idat',0)
      if(.not.ma_pop_stack(l_xp))
     + call errquit('Failed to deallocate xp',0)
      if(.not.ma_pop_stack(l_xadat))
     + call errquit('Failed to deallocate xadat',0)
      if(.not.ma_pop_stack(l_qdat))
     + call errquit('Failed to deallocate qdat',0)
      if(.not.ma_pop_stack(l_xdat))
     + call errquit('Failed to deallocate xdat',0)
      if(.not.ma_pop_stack(l_snam))
     + call errquit('Failed to deallocate sgmnam',0)
      endif
c
      if(.not.ma_verify_allocator_stuff()) print*,'Oops'
c
      if(.not.ma_pop_stack(l_wt))
     + call errquit('Failed to deallocate wt',0)
      if(.not.ma_pop_stack(l_xrms))
     + call errquit('Failed to deallocate xrms',0)
      if(.not.ma_pop_stack(l_xref))
     + call errquit('Failed to deallocate xref',0)
c
      if(.not.ma_pop_stack(l_tag))
     + call errquit('Failed to deallocate tag',0)
      if(.not.ma_pop_stack(l_isel))
     + call errquit('Failed to deallocate isel',0)
c
      return
      end
      subroutine ana_report(sgmnam,xrms,isel)
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "msgids.fh"
c
      character*16 sgmnam(nsa)
      integer isel(nsa)
      real*8 xrms(nsa)
c
      integer i,n
      real*8 arms
c
c     report results
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
     + sqrt(xrms(i)/dble(ndata)),i=1,nsa)
 1001 format(a16,i10,i5,f12.6)
      arms=0.0d0
      n=0
      if(me.eq.0) write(lfnrms,1000)
      do 1 i=1,nsa
      arms=arms+xrms(i)
      n=n+1
      if(i.eq.nsa) goto 2
      if(sgmnam(i)(11:16).ne.sgmnam(i+1)(11:16)) goto 2
      goto 1
    2 continue
      write(lfnrms,1002) sgmnam(i)(1:5),sgmnam(i)(11:16),
     + sqrt(arms/dble(ndata*n))
 1002 format(a5,a6,f12.6)
      arms=0.0d0
      n=0
    1 continue
      close(unit=lfnrms)
      endif
      return
      end
