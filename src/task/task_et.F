C> \ingroup task
C> @{
C>
C> \brief Interface to a charge transfer probability
C>
C> \return Return .true. if successfull, and .false. otherwise.
C>
      logical function task_et(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "stdio.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c
      logical status
      double precision cpu, wall, Vba
      logical task_et_doit
      external task_et_doit
      character*32 theory
c
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     $     call errquit('task:energy: theory not specified',0, RTDB_ERR)

      if ((theory.ne.'scf').and.(theory.ne.'dft')) then
         write(LuOut,*) 'Electron Transfer has been',
     &  ' implemented only for HF and DFT. NWChem will skip this task!.'
         task_et = .false.
         return
      endif
c
      status = task_et_doit(rtdb,Vba)
c
      if (.not. rtdb_put(rtdb, 'task:et', mt_dbl, 1, Vba))
     $     call errquit('task_et: failed storing Vba',0, RTDB_ERR)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $     call errquit('task_et: failed to store status',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_et: failed storing cputime',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_et: failed storing walltime',0, RTDB_ERR)
c
#ifndef SKIP_ECCE
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call geom_print_rtdb_ecce(rtdb)
         call ecce_print1('Vba', mt_dbl, Vba, 1)
         call ecce_print_module_exit('task et', 'ok')
      else
         call ecce_print_module_exit('task et', 'failed')
      endif
#endif
c         
      task_et = status
c
      end
C> @}
c===========================================================================
      logical function task_et_doit(rtdb,Vba)
      implicit none
c
c     This function actually calculates Vba
*
* $Id$
c
*
#include "rtdb.fh"
#include "mafdecls.fh"
c#include "inp.fh"
c#include "util.fh"
c#include "global.fh"
c#include "apiP.fh"
c#include "bas.fh"
c#include "geom.fh"
c#include "cscf.fh"
c
c     Generic NWChem interface to compute Vba. SCF of DFT determinants
c     can be used.
c
c     RTDB output parameters
c     ----------------------
c     task:status (logical)- T/F for success/failure
c     if (status) then
c     .  task:et (real)   - Vba (electron transfer matrix element)
c     .  task:cputime (real)  - cpu time to execute the task
c     .  task:walltime (real) - wall time to execute the task
c
c     Also returns status through the function value
c
      integer rtdb
      logical status
      external et_calc, et_calc_fmo
      logical  et_calc, et_calc_fmo
      double precision Vba
      integer fmo
c
      task_et_doit = .false.
      if (.not.rtdb_get(rtdb,'et:fmo',MT_INT,1,fmo)) then
         status = et_calc(rtdb,Vba)
      else
         status = et_calc_fmo(rtdb, Vba)
      endif
c
      task_et_doit = status
c
      return
      end
