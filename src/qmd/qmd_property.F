c
c     qmd_property: wrapper for property calculations
c
      logical function qmd_property(rtdb)
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
      logical qmd_hasproperty
      external qmd_hasproperty
c
      logical task_property
      external task_property
c
c     preliminaries
      pname = "qmd_property: "
c
      if (qmd_hasproperty(rtdb)) then
        qmd_property = task_property(rtdb)
      end if 
c
      qmd_property = .true.
c
      return
      end
c
c     qmd_hasproperty: check if theory level supports properties
c
      logical function qmd_hasproperty(rtdb)
c
      implicit none
c
#include "rtdb.fh"
#include "errquit.fh"
#include "stdio.fh"
#include "global.fh"
c
      integer rtdb               ! rtdb
c
      character*32 theory
      character*30 pname
c
c     preliminaries
      pname = "qmd_hasproperty: "
c
c     get theory
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $ call errquit(pname//'theory not specified',0, RTDB_ERR)
c
      qmd_hasproperty = .false.
      if (theory .eq. 'dft') qmd_hasproperty = .true.
c
      return
      end
