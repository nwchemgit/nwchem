C> \ingroup task
C> @{
C>
C> \brief Generic interface to a frequency calculation including IR intensities
C>
C> \return Return .true. if the calculation was successfull, and .false. 
C> otherwise.
C>
      logical function task_freq(rtdb)
c     
c     Central difference calculation of the hessian using
c     the generic energy/gradient interface.  Uses a routine inside
c     stepper to do the finite difference ... this needs to be
c     cleaned up to be independent of stepper.
c
c     Also will be hooked up to analytic methods as they are available.
c
c     Since this routine is directly invoked by application modules
c     no input is processed in this routine.
c
c     RTDB input parameters
c     ---------------------
c     task:theory
c
c     RTDB output parameters
c     ----------------------
c     task:hessian file name (string) - name of file containing hessian
c     task:status (logical)      - T/F on success/failure
c     task:cputime
c     task:walltime
c
*
* $Id$
*
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "stdio.fh"
#include "global.fh"
#include "inp.fh"
#include "util.fh"
c
      logical task_hessian
      external task_hessian
c     
      integer rtdb !< [Input] The RTDB handle
c     
      integer nat, geom
      logical status
      character*(nw_max_path_len) filehess
c     
      double precision cpu, wall
c 
      logical      ignore
      logical o_reuse
c    
      call ecce_print_module_entry('task frequencies')
      cpu  = util_cpusec()
      wall = util_wallsec()
c
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .false.))
     $     call errquit('task_freq: failed to invalidate status',0,
     &       RTDB_ERR)
      if (ga_nodeid().eq.0 .and.
     $    util_print('task_freq', print_low)) then
        write(LuOut,*)
        write(LuOut,*)
        call util_print_centered(LuOut,
     $      'NWChem Nuclear Hessian and Frequency Analysis',
     $      40,.true.)
        write(LuOut,*)
      endif
*
      if (rtdb_get(rtdb,'vib:reuse',mt_log,1,ignore)) then
        o_reuse = ignore
      else
        o_reuse = .false.
      endif
*
      if (.not.(o_reuse)) then
        status = task_hessian(rtdb)
        if (.not.status) call errquit
     &      ('task_freq: task_hessian failed',911, CALC_ERR)
      else
        if (ga_nodeid().eq.0)
     &    call util_print_centered(LuOut,
     &        'reusing previously generated Hessian',
     &        40,.true.)
        status = .true.
      endif
*
      ignore = rtdb_parallel(.false.)
      if ((ga_nodeid()).eq.0) then
        if (o_reuse) then
          if (rtdb_cget(rtdb,'vib:reuse_hessian_file',1,filehess)) then
            write(LuOut,*)' re-using hessian in file ',
     &          filehess(1:inp_strlen(filehess))
          else
* in case of manual restart and renaming of hess file
            call util_file_name('hess',  .false., .false.,filehess)
          endif
        else
          if (.not. rtdb_cget(rtdb, 'task:hessian file name', 1,
     $        filehess)) call errquit
     $        ('task_freq: failed reading hessian filename from rtdb',0,
     &       RTDB_ERR)
        endif
c
c     create/load reference geometry
c
        if (.not.geom_create(geom,'geometry')) call errquit
     $      ('task_freq:geom_create failed?',1, GEOM_ERR)
        if (.not.geom_rtdb_load(rtdb,geom,'geometry'))
     $      call errquit
     $      ('task_freq:geom_rtdb_load failed?',2, GEOM_ERR)
        if (.not. geom_ncent(geom,nat)) call errquit
     $      ('task_freq:geom_ncent failed?',3, GEOM_ERR)
        if (.not. geom_destroy(geom)) call errquit
     $      ('task_freq:geom_destroy failed?',911, GEOM_ERR)
        call vib_vib(rtdb,filehess,.true.,
     $      0,.false.,0,.false.,nat)
      endif
      call ga_sync()
      ignore = rtdb_parallel(.true.)
c
      cpu  = util_cpusec() - cpu
      wall = util_wallsec() - cpu
c
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cpu))
     $     call errquit('task_freq: failed storing cputime',0, RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, wall))
     $     call errquit('task_freq: failed storing walltime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     $     call errquit('task_freq: failed to set status',0,
     &       RTDB_ERR)
c
c
      call ecce_print1('cpu time', mt_dbl, cpu, 1)
      call ecce_print1('wall time', mt_dbl, wall, 1)
      task_freq = status
      if (task_freq) then
        call ecce_print_module_exit('task frequencies', 'ok')
      else
        call ecce_print_module_exit('task frequencies', 'failed')
      endif
c
c
      end
C> @}
