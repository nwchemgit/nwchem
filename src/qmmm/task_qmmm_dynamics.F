      logical function task_qmmm_dynamics(rtdb)
*
* $Id: task_qmmm_dynamics.F,v 1.4 2005-08-31 00:00:09 marat Exp $
*
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
      integer rtdb
c
      logical status
      double precision cpu,wall
c
      character*84 tag
      character*30 region(3)
      character*30 pname
      integer nregion
c
      pname = "task_qmmm_dynamics"

c
c
c      call util_print_push() 
c      call util_print_rtdb_load(rtdb, 'qmmm')
c
      call ecce_print_module_entry('task dynamics')
c
      task_qmmm_dynamics = .false.
c
      cpu  = util_cpusec()
      wall = util_wallsec()

c     Set output parameters as if failed
c
      status = rtdb_delete(rtdb, 'task:dynamics')
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_dynamics: failed to invalidate status',0,
     &       RTDB_ERR)
c
c     region definitions
c     ------------------
      tag ="qmmm:region"
      if (.not.rtdb_get(rtdb,tag(1:inp_strlen(tag))//"_n",
     >                 mt_int,1,nregion)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
      if(nregion.gt.3) 
     >      call errquit(pname//"too many regions",0,0)
      if (.not.rtdb_cget(rtdb,tag,nregion,region)) 
     >      call errquit(pname//tag,0,RTDB_ERR)
c
c     define set of active atoms
c     --------------------------
      call qmmm_cons_reset()
      call qmmm_cons_set(rtdb,"fix","solute")
      call qmmm_cons_set(rtdb,"fix","solvent")
      call qmmm_cons_set(rtdb,"free",region(1))
c
      call qmmm_bq_data_update_active()
c
      call qmmm_set_interface("mm")
c
c     Actually do the deed
c
      call md_md()
      status = .true.
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - wall
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, status))
     $   call errquit('task_dynamics: failed to set status',0, RTDB_ERR)
c
      if (status) then
c
         if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $        call errquit('task_dynamics: failed storing cputime',0,
     &       RTDB_ERR)
         if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $        call errquit('task_dynamics: failed storing walltime',0,
     &       RTDB_ERR)
c
      endif
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      if (status) then
         call ecce_print_module_exit('task dynamics', 'ok')
      else
         call ecce_print_module_exit('task dynamics', 'failed')
      endif
c         
      task_qmmm_dynamics = status
c
c      call util_print_pop()
c
      end

