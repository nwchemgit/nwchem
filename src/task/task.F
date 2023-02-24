c
c $Id$
c
C> \ingroup task
C> @{
c
C> \brief The main task driver called for any task.
C>
C> This routine is called by ALL processes.  Task_input has read the
C> task directive and put stuff into the database.  Get the
C> data out and invoke the desired action.
C>
C> If the operation is in the list of those supported by generic 
C> routines (currently energy, gradient, hessian, optimize), then 
C> the generic routine is called.  Otherwise, a match is attempted 
C> for a specialized routine.  If no operation is specified 
C> and no specialized routine located, then it is assumed that
C> a generic energy calculation is required.
c
      subroutine task(rtdb)
c
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "util.fh"
#include "global.fh"
#include "mafdecls.fh"
#include "stdio.fh"
      integer rtdb              !< [Input] the RTDB handle
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
      character*32 key
      logical status, ignore
      character*80 errmsg
      double precision task_cpu, task_wall
      logical oprint_task_times, odebug
      character*32 bqname
c
      logical ocons
      logical rjhtest, raktest, atwtest, task_check_sum, cphf,
     $     jantest, uccsdtest
      logical property, dplot, gapss, task_shell, moints_full
      logical argos,analyz,diana,esp,prepar,nwmd
      logical kgdtest, jvltest, jvltest2
      logical task_energy, task_gradient, task_optimize, task_freq
      logical task_mepgs, task_string
      logical task_hessian, task_saddle, task_lst, task_neb
      logical task_dynamics, task_rfcf
      logical task_qmd
      logical orism
      logical omm
      logical oqmmm,opes, task_python, task_thermo,task_bq,task_nbo
      logical osmd, ofragment
c     logical oembed
      logical task_version, task_vscf
      logical task_raman 
      logical task_paw 
      logical task_pspw,task_band,task_ofpw,task_smd
      logical task_drdy
      logical task_et
      logical qmmm_init
      logical qmmm_end
c     logical embed_init
      logical cons_active
      logical task_property
      logical task_ccca
      logical lucia
c
c      logical task_pes
c      external task_pes
c
      logical task_smd_fragment
      external task_smd_fragment
      external qmmm_init
c
c     external embed_init
c
      external cons_active
      external qmmm_end
      external jantest, uccsdtest
      external rjhtest, raktest, atwtest, task_check_sum, cphf
      external property, dplot, gapss, task_shell, moints_full
      external analyz,diana,esp,prepar,nwmd,argos
      external task_energy, task_gradient, task_optimize, task_freq
      external task_mepgs, task_string
      external task_hessian, task_saddle, task_lst, task_neb
      external task_dynamics, task_rfcf
      external task_qmd
      external task_version, task_python, task_thermo,task_nbo
      external kgdtest, task_vscf
      external task_raman 
      external task_paw 
      external task_pspw,task_band,task_ofpw,task_smd
      external task_drdy
      external task_et
      external task_property
      external task_ccca
      external lucia
c
      external mmi_end
      logical mmi_end
      external mmi_init
      logical mmi_init
c
      logical bq_init,bq_end
      external bq_init,bq_end
c
      logical  task_qmmm_esp
      external task_qmmm_esp
c
      logical task_qmmm_fep
      external task_qmmm_fep
c
      logical task_qmmm_abc
      external task_qmmm_abc
c
      logical task_qmmm_pes
      external task_qmmm_pes
c
      logical task_qmmm_optimize
      external task_qmmm_optimize
c
      logical task_smd_dynamics
      external task_smd_dynamics
c
      logical task_qmmm_dynamics
      external task_qmmm_dynamics
c
      logical task_dntmc
      external task_dntmc

      logical task_rt_tddft
      external task_rt_tddft

      logical task_fcidump
      external task_fcidump

      logical task_gwmol
      external task_gwmol

      logical task_aoints
      external task_aoints

      logical task_grdtrc
      external task_grdtrc

      logical task_nadtrc
      external task_nadtrc

c      logical embed_optimize
c      external embed_optimize
c
c      logical embed_prepare
c      external embed_prepare
c
c      logical embed_frequency
c      external embed_frequency
c
C#ifdef USE_SUBGROUPS
C      logical util_sgstart
C      external util_sgstart
C#endif
      logical do_cosmo, flag
c
      oprint_task_times = util_print('task time',print_default)
      task_wall = -util_wallsec()
      task_cpu  = -util_cpusec()
*
C#ifdef USE_SUBGROUPS
C      if(util_sgstart(rtdb)) then
C#endif
c
      if (rtdb_get(rtdb, 'task:debug', mt_log, 1, odebug))
     $     call util_debug(rtdb)
C
      if (.not.rtdb_get(rtdb,'task:mm',mt_log,1,omm))
     &    omm = .false.
C
      if (.not.rtdb_get(rtdb,'task:QMMM',mt_log,1,oqmmm))
     &    oqmmm = .false.
c      
      if (.not.rtdb_get(rtdb,'task:smd',mt_log,1,osmd))
     &    osmd = .false.
c      
      if (.not.rtdb_get(rtdb,'task:fragment',mt_log,1,ofragment))
     &    ofragment = .false.
c
      if (.not.rtdb_get(rtdb,'task:pes',mt_log,1,opes))
     &    opes = .false.
c
c      if (.not.rtdb_get(rtdb,'task:embed',mt_log,1,oembed))
c     &    oembed = .false.
c
      if (rtdb_cget(rtdb,'bq',1,bqname)) then
         task_bq = .true.
      else
         task_bq = .false. 
      end if
c      
      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = ' '
c
      if(.not.rtdb_cget(rtdb,'task:theory',1,theory))
     + call errquit('task: no task input for theory?',0, RTDB_ERR)
c
      if (.not. rtdb_get(rtdb, 'task:ignore', mt_log, 1, ignore))
     $     ignore = .false.
c
c     If this task is going to use cosmo then we need to (re)run the
c     wavefunction calculation. To force this we set the 
c     <theory>:converged flag to .false. if it exists. Otherwise the
c     code will report that the wavefunction is already converged
c     and exit the optimizer rather than calculating the solvation
c     energy.
c
      if (.not. rtdb_get(rtdb, 'slv:cosmo', mt_log, 1, do_cosmo))
     $     do_cosmo = .false.
      if (do_cosmo) then
        key = theory(1:inp_strlen(theory))//":converged"
        if (rtdb_get(rtdb, key, mt_log, 1, flag)) then
          if (.not. rtdb_put(rtdb, key, mt_log, 1, .false.))
     $      call errquit("task: could not write to rtdb",0,RTDB_ERR)
        endif
      endif
c
c     if needed initialize qmmm for the duration of the task
c
      if (oqmmm) then
        status = qmmm_init(rtdb)
        if ((ga_nodeid().eq.0) .and. (.not. status)) then
              call errquit("failed to initialize qmmm", 0, CALC_ERR)
        end if
      endif
c
      if (omm) then
        status = mmi_init(rtdb)
        if ((ga_nodeid().eq.0) .and. (.not. status)) then
              call errquit("failed to initialize mm", 0, CALC_ERR)
        end if
      endif
cc
c      if (oembed) then
c        status = embed_init(rtdb)
c        if ((ga_nodeid().eq.0) .and. (.not. status)) then
c              call errquit("failed to initialize embed", 0, CALC_ERR)
c        end if
c      endif
c
      ocons = cons_active(rtdb)
      if(ocons) call cons_init(rtdb)
c
      if(task_bq) then
        status = bq_init(rtdb)
        if ((ga_nodeid().eq.0) .and. (.not. status)) then
              call errquit("failed to initialize bq", 0, CALC_ERR)
        end if
      end if
c
c 
c
      if(osmd.and.ofragment) then
         status = task_smd_fragment(rtdb)
         call ga_sync()
         goto 99999
      end if
      if(oqmmm.and.opes) then
         status = task_qmmm_pes(rtdb)
         call ga_sync()
         goto 99999
      endif
c
c      if(opes.and.(.not.oqmmm)) then
c         status = task_pes(rtdb)
c         call ga_sync()
c         goto 99999
c      endif
c
      if (inp_compare(.false.,'ccca',theory)) then
         status = task_ccca(rtdb)
         call ga_sync()
         goto 99999
      endif
c
c Look for speciallized routine tasks here 
c
      if (operation .eq. ' ') then ! Look for specialized routine
         if (inp_compare(.false.,'rjhtest', theory)) then
            status = rjhtest(rtdb)
         else if (inp_compare(.false.,'jantest', theory)) then
            status = jantest(rtdb)
         else if (inp_compare(.false.,'test', theory)) then
            call util_test(rtdb)
            status = .true.
         else if (inp_compare(.false.,'uccsdtest', theory)) then
            status = uccsdtest(rtdb)
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
         else if (inp_compare(.false.,'md', theory)) then
            status = nwmd(rtdb)
         else if (inp_compare(.false.,'argos', theory)) then
            status = argos(rtdb)
         else if (inp_compare(.false.,'cphf', theory)) then
            status = cphf(rtdb)
         else if (inp_compare(.false.,'nbo', theory)) then
            status = task_nbo(rtdb)
         else if (inp_compare(.false.,'drdy', theory).or.
     &            inp_compare(.false.,'dirdyvtst', theory).or.
     &            inp_compare(.false.,'drdy_nwchem', theory)) then
            status = task_drdy(rtdb)
         else if (inp_compare(.false.,'shell', theory)) then
            status = task_shell(rtdb)
         else if (inp_compare(.false.,'moints', theory)) then
            status = moints_full(rtdb)
         else if (inp_compare(.false.,'dplot', theory)) then
            status = dplot(rtdb)
         elseif(inp_compare(.false.,'esp',theory)) then
            if(oqmmm) then
              status = task_qmmm_esp(rtdb)
            else
              status = esp(rtdb)
            end if
         elseif(inp_compare(.false.,'prepare',theory)) then
          status=prepar(rtdb)
         elseif(inp_compare(.false.,'analysis',theory)) then
          status=analyz(rtdb)
         elseif(inp_compare(.false.,'analyze',theory)) then
          status=analyz(rtdb)
         elseif(inp_compare(.false.,'diana',theory)) then
          status=diana(rtdb)
         elseif(inp_compare(.false.,'lucia',theory)) then
          status=lucia(rtdb)
         else
            operation = 'energy' ! No matching specialized operation
            goto 100            !  execute the default generic operation
         endif
         if ((ga_nodeid().eq.0) .and. (.not. status)) then
            errmsg = ' '
            write(errmsg,'(a,'' failed '')')theory(1:inp_strlen(theory))
            if (ga_nodeid().eq.0) then
               if (ignore) then
                  write(LuOut,*) ' !! warning: ',
     1                        errmsg(1:inp_strlen(errmsg))
                  call util_flush(6)
               else
                  call errquit(errmsg, 0, CALC_ERR)
               endif
            endif
         endif
         call ga_sync()
         goto 99999
      endif
c
 100  if (inp_compare(.false., 'energy', operation)) then
         status = task_energy(rtdb)
      else if (inp_compare(.false., 'gradient', operation)) then
         status = task_gradient(rtdb)
      else if (inp_compare(.false., 'aoints', operation)) then
         status = task_aoints(rtdb)
      else if (inp_compare(.false., 'grdtrc', operation)) then
         status = task_nadtrc(rtdb)
      else if (inp_compare(.false., 'nadtrc', operation)) then
         status = task_nadtrc(rtdb)
      else if (inp_compare(.false., 'property', operation)) then
         status = task_property(rtdb)
      else if (inp_compare(.false., 'hessian', operation)) then
         status = task_hessian(rtdb)
      else if (inp_compare(.false., 'optimize', operation)) then
         if(oqmmm) then
           status = task_qmmm_optimize(rtdb)
c        else if (oembed) then
c          status = embed_optimize(rtdb)
         else
           status = task_optimize(rtdb)
         endif
      else if (inp_compare(.false., 'saddle', operation)) then
         if(oqmmm) then
           status = task_qmmm_optimize(rtdb)
c        else if (oembed) then
c          status = embed_optimize(rtdb)
         else
           status = task_saddle(rtdb)
         endif
      else if(inp_compare(.false.,'neb',operation)) then
         status=task_neb(rtdb)
      else if (inp_compare(.false., 'lst', operation)) then
         status = task_lst(rtdb)
      else if (inp_compare(.false., 'et', operation)) then
         status = task_et(rtdb)
      else if (inp_compare(.false., 'rt_tddft', operation)) then
         status = task_rt_tddft(rtdb)
      else if (inp_compare(.false., 'fcidump', operation)) then
         status = task_fcidump(rtdb)
      else if (inp_compare(.false., 'freq', operation)) then
c        if (oembed) then
c          status = embed_frequency(rtdb)
c        else
           status = task_freq(rtdb)
c        end if
      else if (inp_compare(.false.,'mepgs', operation)) then
         status = task_mepgs(rtdb)
      else if (inp_compare(.false.,'string', operation)) then
         status = task_string(rtdb)
      else if (inp_compare(.false.,'dynamics',operation)) then
         if(oqmmm) then
           status = task_qmmm_dynamics(rtdb)
         else if(osmd) then
           status = task_smd_dynamics(rtdb)
         else
           status = task_dynamics(rtdb)
         endif
      else if (inp_compare(.false.,'qmd',operation)) then
           status = task_qmd(rtdb)
      else if (inp_compare(.false.,'rfcf',operation)) then
           status = task_rfcf(rtdb)
      else if(inp_compare(.false.,'fep',operation)) then
         if(oqmmm) then
           status = task_qmmm_fep(rtdb)
         else
           call errquit('task: fep needs QMMM',0, INPUT_ERR)
         end if
      else if(inp_compare(.false.,'abc',operation)) then
         if(oqmmm) then
           status = task_qmmm_abc(rtdb)
         else
           call errquit('task: abc needs QMMM',0, INPUT_ERR)
         end if
      else if(inp_compare(.false.,'thermodynamics',operation)) then
          status=task_thermo(rtdb)
c     else if(inp_compare(.false.,'prepare',operation)) then
c         status=embed_prepare(rtdb)
      else if (inp_compare(.false.,'raman',operation)) then
         status=task_raman(rtdb)
      else if(inp_compare(.false.,'vscf',operation)) then
         status=task_vscf(rtdb)
      else if(inp_compare(.false.,'dntmc',operation)) then
         status=task_dntmc(rtdb)
      else if (inp_compare(.false., 'gw', operation)) then
         status = task_gwmol(rtdb)
c
c     **** check for theories that have specific operations ****
c
c     *** task_pspw handles pspw specific "operations" ***
      else if (inp_compare(.false.,'pspw',theory)) then
        status = task_pspw(rtdb)
c
c     *** task_band handles band specific "operations" ***
      else if (inp_compare(.false.,'band',theory)) then
        status = task_band(rtdb)

c     *** task_paw handles paw specific "operations" ***
      else if (inp_compare(.false.,'paw',theory)) then
        status = task_paw(rtdb)

c     *** task_ofpw handles ofpw specific "operations" ***
      else if (inp_compare(.false.,'ofpw',theory)) then
        status = task_ofpw(rtdb)

c     *** task_smd handles smd specific "operations" ***
      else if (inp_compare(.false.,'smd',theory)) then
        write(*,*) "task_smd is disabled for now"
c        write(*,*) "going to task_smd"
c        status = task_smd(rtdb)

      else
         call errquit('task: unknown generic operation',0, INPUT_ERR)
      endif
      if ((ga_nodeid().eq.0) .and. (.not. status)) then
         errmsg = ' '
         write(errmsg,'(a,'' '',a,'' failed '')')
     $      theory(1:inp_strlen(theory)), 
     $      operation(1:inp_strlen(operation))
         if (ga_nodeid().eq.0) then
            if (ignore) then
               write(LuOut,*)' !! warning: ',
     1                 errmsg(1:inp_strlen(errmsg))
               call util_flush(6)
            else
               call errquit(errmsg, 0, CALC_ERR)
            endif
         endif
      endif
cc
      call ga_sync()
99999 continue
c
      if (oqmmm) then
        status = qmmm_end(rtdb)

        if ((ga_nodeid().eq.0) .and. (.not. status)) then
              call errquit("failed qmmm_end", 0, CALC_ERR)
        end if
        status=rtdb_delete(rtdb, 'task:QMMM')

      endif

      if(ocons) call cons_end(rtdb)
c
      if(task_bq) then
       status = bq_end(rtdb)
       if ((ga_nodeid().eq.0) .and. (.not. status)) then
              call errquit("failed to uninitialize bq", 0, CALC_ERR)
        end if
      end if

      if(omm) then
       status = mmi_end(rtdb)
       if ((ga_nodeid().eq.0) .and. (.not. status)) then
           call errquit("failed to end mm", 0, CALC_ERR)
        end if
      end if



      task_wall = task_wall + util_wallsec()
      task_cpu  = task_cpu  + util_cpusec()
      if (oprint_task_times) then
        if (ga_nodeid().eq.0) then
          write(LuOut,00001) task_cpu, task_wall
00001     format(/' Task  times  cpu:',f11.1,'s     wall:',f11.1,'s')
          call util_flush(6)
        endif
      endif
C#ifdef USE_SUBGROUPS
C      endif
C      call util_sgend(rtdb)
C#endif
      call ecce_print1('single task cpu time',mt_dbl,task_cpu,1)
      call ecce_print1('single task wall time',mt_dbl,task_wall,1)
      call ga_sync()
c
      return
c
      end
c
C> \brief Performs the "task version" operation
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function task_version(rtdb)
      implicit none
      integer rtdb !< [Input] the RTDB handle
      call util_version()
      task_version = .true.
      return
      end
c
C> \brief Performs the "task shell" operation
C>
C> \return Return .true. if successfull, and .false. otherwise
c
      logical function task_shell(rtdb)
      implicit none
#include "errquit.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer rtdb !< [Input] the RTDB handle
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
     $     call errquit('task shell: shell command missing',0, RTDB_ERR)
c
      task_shell = .true.
c
      do i = range(1), range(2), range(3)
         if (ga_nodeid() .eq. i)
     $        task_shell = util_system(command) .eq. 0
      enddo
c
      return
      end
C> @}
