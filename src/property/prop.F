c***********************************************************************
      logical function property(rtdb)
c***********************************************************************
c
c>>>  Inquire rtdb for properties to be evaluated.
c
*
* $Id: prop.F,v 1.10 2003-03-27 18:55:44 bert Exp $
*
      implicit none
      integer rtdb              ! [input]
      integer nbofile
      logical status
      logical  hnd_property
      external hnd_property
#include "context.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "global.fh"
#include "stdio.fh"
c
c     Push context down to prop.
c
      call ga_sync()

      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
*      if(ga_nodeid().eq.0)
*     $     call util_print_centered(LUout,'Properties',40,.true.)
c
c     start ecce property output module
c
      call ecce_print_module_entry('task property')
c
c     hondo property calculation routines 
c
      status = hnd_property(rtdb)
c
c     nbo input file
c
      if (rtdb_get(rtdb,'prop:nbofile',MT_INT,1,nbofile)) then
        if(nbofile.eq.1) call wnbofile(rtdb)
      endif 
c
c     finish ecce property output module
c
      call ecce_print_module_exit('task property','ok')
      property = .true.
      return
      end
