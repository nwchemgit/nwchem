c
c     task quantum molecular dynamics
c
      logical function task_qmd(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "geom.fh"
#include "stdio.fh"
c
      integer rtdb
c
      logical status
      double precision cpu,wall
c
      character*32 theory
      character*84 tag
      character*30 pname
c
      logical qmd_driver
      external qmd_driver
c
      pname = "task_qmd: "
c
c     start clock
      cpu  = util_cpusec()
      wall = util_wallsec()
c
c     invoke the qmd driver
c     if (ga_nodeid().eq.0) write(luout,*) "calling qmd_driver"
      status = qmd_driver(rtdb)
c
c     stop clock
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (status) then
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit(pname//'failed storing cputime',0,RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit(pname//'failed storing walltime',0,RTDB_ERR)
      endif
c         
      task_qmd = status
c
      end
c $Id$
