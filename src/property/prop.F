c***********************************************************************
      logical function property(rtdb)
c***********************************************************************
c
c>>>  Inquire rtdb for properties to be evaluated.
c
*
* $Id: prop.F,v 1.7 1997-06-24 16:56:32 mdupuis Exp $
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
      status = rtdb_parallel(.true.) ! Broadcast reads to all processes
      if(ga_nodeid().eq.0)
     $     call util_print_centered(LUout,'Properties',40,.true.)
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
      property = .true.
      return
      end
