c***********************************************************************
      subroutine property(rtdb)
c***********************************************************************
c
c>>>  Inquire rtdb for properties to be evaluated.
c
      implicit none
      integer rtdb              ! [input]
      integer nbofile
      logical status

#include "context.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
c
c     Push context down to prop.
c
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      if(ga_nodeid().eq.0)
     $     call util_print_centered(LUout,'Properties',40,.true.)

      if (.not. context_push('prop'))
     $     call errquit('prop: context_push failed', 0)

      if (rtdb_get(rtdb,'prop:nbofile',MT_INT, 1, nbofile)) then
        if(nbofile.eq.1)call wnbofile(rtdb)
      endif 

      if (.not. context_pop('prop')) 
     &     call errquit('prop: context_pop failed',0)
      end

