      subroutine ana_task
c
c $Id: ana_task.F,v 1.1 1999-01-29 18:10:43 d3j191 Exp $
c
      implicit none
c
#include "ana_common.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      logical ana_rdfram
      external ana_rdfram
c
      integer nsai,msai
      real*8 rms0,rms1
c
      integer ifr,ito,num
      integer nframe
c
      logical ignore_status
c
c
c     set sequential access to rtdb
c     and ignore current state
c
      ignore_status = rtdb_parallel(.false.)
c
      call ana_rtdbin()
c
      if(ltrj.gt.0) call ana_trajct
      if(lplt.gt.0) call ana_pltgrd
c
c     set parallel access to rtdb
c
      ignore_status = rtdb_parallel(.true.)
c
      return
      end

