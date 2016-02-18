c
c     qmd_energy: wrapper for energy
c
      subroutine qmd_energy(rtdb,switchtheory)
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
      character(*) switchtheory
c
      character*255 origtheory
      character*30 pname
c
      logical task_energy
      external task_energy
c
c     preliminaries
      pname = "qmd_energy: "
c
      if(.not.rtdb_cget(rtdb,'task:theory',1,origtheory))
     $   call errquit(pname//'could not get theory',0,RTDB_ERR)
c     Change theory and calculate energy
      if(.not.rtdb_cput(rtdb,'task:theory',1,switchtheory))
     &   call errquit(pname//'failed to change theory',0,RTDB_ERR)
      if(.not.task_energy(rtdb))
     &   call errquit(pname//'task_energy',0,CALC_ERR)
c     Switch back to original theory
      if(.not.rtdb_cput(rtdb,'task:theory',1,origtheory))
     &   call errquit(pname//'failed to change theory',0,RTDB_ERR)
c
      return
      end
