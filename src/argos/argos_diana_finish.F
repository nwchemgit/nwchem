      subroutine argos_diana_finish(sgmnam,ltask)
c
c $Id$
c
      implicit none
c
#include "argos_diana_common.fh"
#include "mafdecls.fh"
c
      character*16 sgmnam(nsa)
      logical ltask
c
c      call argos_diana_report(sgmnam,dbl_mb(i_xrms),int_mb(i_isel))
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
