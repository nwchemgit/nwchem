c
c     qmd_property: wrapper for property calculations
c
      subroutine qmd_property(rtdb)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
c
      integer rtdb               ! rtdb
c
      character*30 pname
c
      logical task_property
      external task_property
c
c     preliminaries
      pname = "qmd_property: "
c
      if(.not.task_property(rtdb))
     &   call errquit(pname//'task_property',0,CALC_ERR)
c
      return
      end
