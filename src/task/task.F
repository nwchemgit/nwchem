      subroutine task(rtdb)
*
* $Id: task.F,v 1.38 2000-01-06 03:11:13 d3g681 Exp $
*
      implicit none
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
      integer rtdb              ! [input] data base handle
c
c     Called by ALL processes.  task_input has read the
c     task directive and put stuff into the database.  Get the
c     data out and invoke the desired action.
c
c     If the operation is in the list of those supported by generic 
c     routines (currently energy, gradient, hessian, optimize), then 
c     the generic routine is called.  Otherwise, a match is attemped 
c     for a specialized routine.  If no operation is specified 
c     and no specialized routine located, then it is assumed that
c     a generic energy calculation is required.
c
c     This needs extending to accomodate QM/MM and other mixed methods
c     by having both MM and QM pieces specified (e.g., task md dft).
c
      character*30 operation
      character*32 theory
      logical status, ignore
      character*80 errmsg
      double precision task_cpu, task_wall
      logical oprint_task_times, odebug
c
      logical rjhtest, raktest, atwtest, task_check_sum, mdmain, cphf
      logical nwargos, property, dplot, gapss, task_shell, moints_full
      logical analyz,esp,prepar,nwmd, kgdtest, jvltest, jvltest2
      logical task_energy, task_gradient, task_optimize, task_freq
      logical task_hessian, task_saddle, task_lst,task_dynamics
      logical task_qmmm, task_python, task_nbo,task_thermo
      logical task_version
      logical task_pspw
      logical task_drdy
c
      external rjhtest, raktest, atwtest, task_check_sum, mdmain, cphf
      external nwargos, property, dplot, gapss, task_shell, moints_full
      external analyz,esp,prepar,nwmd
      external task_energy, task_gradient, task_optimize, task_freq
      external task_hessian, task_saddle, task_lst,task_dynamics
      external task_version, task_python, task_nbo,task_thermo
      external kgdtest
      external task_pspw
      external task_drdy
c
      oprint_task_times = util_print('task time',print_default)
      task_wall = -util_wallsec()
      task_cpu  = -util_cpusec()
*
      if (rtdb_get(rtdb, 'task:debug', mt_log, 1, odebug))
     $     call util_debug(rtdb)
*
      if (.not.rtdb_get(rtdb,'task:QMMM',mt_log,1,task_qmmm))
     &    task_qmmm = .false.
c      
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = ' '
c
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     + call errquit('task: no task input for theory?',0)
c
      if (.not. rtdb_get(rtdb, 'task:ignore', mt_log, 1, ignore))
     $     ignore = .false.
c
c for qmmm start in nwargos.  nwargos will call appropriate task routines
c
      if (task_qmmm) then
        if (theory.eq.'nwmd'.or.theory.eq.'md'.or.theory.eq.'nwargos')
     + call errquit('task:theory must be a QM method for QMMM',911)
*        status = nwargos(rtdb)
        status = nwmd(rtdb)
        if ((ga_nodeid().eq.0) .and. (.not. status)) then
          errmsg = ' '
          write(errmsg,'(''QMMM '',a,'' failed '')')
     &        theory(1:inp_strlen(theory))
          if (ga_nodeid().eq.0) then
            if (ignore) then
              write(6,*)' !! warning: ',errmsg(1:inp_strlen(errmsg))
              call util_flush(6)
            else
              call errquit(errmsg, 0)
            endif
          endif
        endif
        call ga_sync()
        goto 99999
      endif
c
c Look for speciallized routine tasks here 
c
      if (operation .eq. ' ') then ! Look for specialized routine
         if (inp_compare(.false.,'rjhtest', theory)) then
            status = rjhtest(rtdb)
         else if (inp_compare(.false.,'jvltest', theory)) then
            status = jvltest(rtdb)
         else if (inp_compare(.false.,'marga', theory)) then
            status = jvltest2(rtdb)
         else if (inp_compare(.false.,'python', theory)) then
            status = task_python(rtdb)
         else if (inp_compare(.false.,'raktest', theory)) then
            status = raktest(rtdb)
         else if (inp_compare(.false.,'kgdtest', theory)) then
            status = kgdtest(rtdb)
         else if (inp_compare(.false.,'version', theory)) then
            status = task_version(rtdb)
         else if (inp_compare(.false.,'atwtest', theory)) then
            status = atwtest(rtdb)
         else if (inp_compare(.false.,'rtdbprint', theory)) then
            status = rtdb_print(rtdb,.true.)
         else if (inp_compare(.false.,'checksum', theory)) then
            status = task_check_sum(rtdb)
         else if (inp_compare(.false.,'nwargos', theory)) then
            status = nwargos(rtdb)
         else if (inp_compare(.false.,'md', theory)) then
            status = nwmd(rtdb)
         else if (inp_compare(.false.,'md_ideaz', theory)) then
            status = mdmain(rtdb)
         else if (inp_compare(.false.,'cphf', theory)) then
            status = cphf(rtdb)
         else if (inp_compare(.false.,'property', theory)) then
            status = property(rtdb)
         else if (inp_compare(.false.,'drdy', theory).or.
     &            inp_compare(.false.,'drdy_nwchem', theory)) then
            status = task_drdy(rtdb)
         else if (inp_compare(.false.,'nbo', theory)) then
            status = task_nbo(rtdb)
         else if (inp_compare(.false.,'gapss', theory)) then
            status = gapss(rtdb) ! This could/should be in task_energy
         else if (inp_compare(.false.,'shell', theory)) then
            status = task_shell(rtdb)
         else if (inp_compare(.false.,'moints', theory)) then
            status = moints_full(rtdb)
         else if (inp_compare(.false.,'dplot', theory)) then
            status = dplot(rtdb)
         elseif(inp_compare(.false.,'esp',theory)) then
          status=esp(rtdb)
         elseif(inp_compare(.false.,'prepare',theory)) then
          status=prepar(rtdb)
         elseif(inp_compare(.false.,'analysis',theory)) then
          status=analyz(rtdb)
         elseif(inp_compare(.false.,'analyze',theory)) then
          status=analyz(rtdb)
         else
            operation = 'energy' ! No matching specialized operation
            goto 100            !  execute the default generic operation
         endif
         if ((ga_nodeid().eq.0) .and. (.not. status)) then
            errmsg = ' '
            write(errmsg,'(a,'' failed '')')theory(1:inp_strlen(theory))
            if (ga_nodeid().eq.0) then
               if (ignore) then
                  write(6,*)' !! warning: ',errmsg(1:inp_strlen(errmsg))
                  call util_flush(6)
               else
                  call errquit(errmsg, 0)
               endif
            endif
         endif
         call ga_sync()
         goto 99999
      endif
c
c     *** task_pspw handles pspw specific "operations" ***
 100  if (inp_compare(.false.,'pspw',theory)) then
        status = task_pspw(rtdb)
      else if (inp_compare(.false., 'energy', operation)) then
         status = task_energy(rtdb)
      else if (inp_compare(.false., 'gradient', operation)) then
         status = task_gradient(rtdb)
      else if (inp_compare(.false., 'hessian', operation)) then
         status = task_hessian(rtdb)
      else if (inp_compare(.false., 'optimize', operation)) then
         status = task_optimize(rtdb)
      else if (inp_compare(.false., 'saddle', operation)) then
         status = task_saddle(rtdb)
      else if (inp_compare(.false., 'lst', operation)) then
         status = task_lst(rtdb)
      else if (inp_compare(.false., 'freq', operation)) then
         status = task_freq(rtdb)
      else if (inp_compare(.false.,'dynamics',operation)) then
         status=task_dynamics(rtdb)
      else if(inp_compare(.false.,'thermodynamics',operation)) then
         status=task_thermo(rtdb)
      else
         call errquit('task: unknown generic operation',0)
      endif
      if ((ga_nodeid().eq.0) .and. (.not. status)) then
         errmsg = ' '
         write(errmsg,'(a,'' '',a,'' failed '')')
     $      theory(1:inp_strlen(theory)), 
     $      operation(1:inp_strlen(operation))
         if (ga_nodeid().eq.0) then
            if (ignore) then
               write(6,*)' !! warning: ',errmsg(1:inp_strlen(errmsg))
               call util_flush(6)
            else
               call errquit(errmsg, 0)
            endif
         endif
      endif
      call ga_sync()
99999 continue
c
      task_wall = task_wall + util_wallsec()
      task_cpu  = task_cpu  + util_cpusec()
      if (oprint_task_times) then
        if (ga_nodeid().eq.0) then
          write(6,00001) task_cpu, task_wall
00001     format(/' Task  times  cpu:',f11.1,'s     wall:',f11.1,'s')
          call util_flush(6)
        endif
      endif
      call ecce_print1('single task cpu time',mt_dbl,task_cpu,1)
      call ecce_print1('single task wall time',mt_dbl,task_wall,1)
      call ga_sync()
c
      return
c
      end
      logical function task_version(rtdb)
      implicit none
      integer rtdb
      call util_version()
      task_version = .true.
      end
      logical function task_shell(rtdb)
      implicit none
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb
c
      integer i
      character*1023 command
      integer range(3)
      integer util_system
      external util_system
c
      if (.not.rtdb_get(rtdb,'task:shell range',mt_int,3,range)) then
         range(1) = 0
         range(2) = 0
         range(3) = 1
      endif
      if (.not. rtdb_cget(rtdb,'task:shell command',1,command))
     $     call errquit('task shell: shell command missing',0)
c
      task_shell = .true.
c
      do i = range(1), range(2), range(3)
         if (ga_nodeid() .eq. i)
     $        task_shell = util_system(command) .eq. 0
      enddo
c
      end



