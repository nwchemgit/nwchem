c
c     Task handler for real time TDDFT
c
      logical function task_rt_tddft(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
#include "global.fh"
c
      integer rtdb
c
      logical status
      double precision cputime, walltime
      logical dft_rt_tddft
      logical sodft_rt_tddft
      external dft_rt_tddft
      external sodft_rt_tddft
      character*32 theory
c
c     ! Time at the start
      cputime  = util_cpusec()
      walltime = util_wallsec()
c
c     ! Perform the task
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     &  call errquit('task:rt_tddft: theory not specified',0, RTDB_ERR)
      
      if (theory.eq.'dft') then
         status = dft_rt_tddft(rtdb)
      elseif (theory.eq.'sodft') then
         status = sodft_rt_tddft(rtdb)
      else
         status = .false.
      end if
      
c
c     ! Time at the end
      cputime  = util_cpusec() - cputime
      walltime = util_wallsec() - walltime
c
      task_rt_tddft = status
c
      end
c $Id$
