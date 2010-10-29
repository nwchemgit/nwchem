      subroutine dia_finish(sgmnam,ltask)
c
c $Id$
c
      implicit none
c
#include "dia_common.fh"
#include "mafdecls.fh"
c
      character*16 sgmnam(nsa)
      logical ltask
c
c      call dia_report(sgmnam,dbl_mb(i_xrms),int_mb(i_isel))
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
c      call ma_summarize_allocated_blocks()
      if(ltask) then
      if(.not.ma_pop_stack(l_atmn))
     + call md_abort('Failed to deallocate atmn',0)
      if(.not.ma_pop_stack(l_ndxw))
     + call md_abort('Failed to deallocate ndxw',0)
      if(.not.ma_pop_stack(l_wsel))
     + call md_abort('Failed to deallocate wsel',0)
      if(.not.ma_pop_stack(l_qwdat))
     + call md_abort('Failed to deallocate qwdat',0)
      if(.not.ma_pop_stack(l_wdat))
     + call md_abort('Failed to deallocate wdat',0)
      if(.not.ma_pop_stack(l_idat))
     + call md_abort('Failed to deallocate idat',0)
      if(.not.ma_pop_stack(l_xp))
     + call md_abort('Failed to deallocate xp',0)
      if(.not.ma_pop_stack(l_xadat))
     + call md_abort('Failed to deallocate xadat',0)
      if(.not.ma_pop_stack(l_qdat))
     + call md_abort('Failed to deallocate qdat',0)
      if(.not.ma_pop_stack(l_xdat))
     + call md_abort('Failed to deallocate xdat',0)
      if(.not.ma_pop_stack(l_snam))
     + call md_abort('Failed to deallocate sgmnam',0)
      endif
c
      if(.not.ma_verify_allocator_stuff()) print*,'Oops'
c
      if(.not.ma_pop_stack(l_owt))
     + call md_abort('Failed to deallocate owt',0)
      if(.not.ma_pop_stack(l_swt))
     + call md_abort('Failed to deallocate swt',0)
      if(.not.ma_pop_stack(l_wt))
     + call md_abort('Failed to deallocate wt',0)
      if(.not.ma_pop_stack(l_xrms))
     + call md_abort('Failed to deallocate xrms',0)
      if(.not.ma_pop_stack(l_xref))
     + call md_abort('Failed to deallocate xref',0)
c
      if(.not.ma_pop_stack(l_val))
     + call md_abort('Failed to deallocate val',0)
      if(.not.ma_pop_stack(l_tag))
     + call md_abort('Failed to deallocate tag',0)
      if(.not.ma_pop_stack(l_osel))
     + call md_abort('Failed to deallocate osel',0)
      if(.not.ma_pop_stack(l_ssel))
     + call md_abort('Failed to deallocate ssel',0)
      if(.not.ma_pop_stack(l_isel))
     + call md_abort('Failed to deallocate isel',0)
      if(.not.ma_pop_stack(l_hbnd))
     + call md_abort('Failed to deallocate hbnd',0)
c
      return
      end
      subroutine dia_report(sgmnam,xrms,isel)
c
      implicit none
c
#include "dia_common.fh"
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
