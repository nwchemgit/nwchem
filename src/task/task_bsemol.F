c
c     Task handler for Molecular BSE calcualtions
c     Daniel Mejia-Rodriguez, 2022
c
      logical function task_bsemol(rtdb)
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
      character*32 theory
      character*255 title
      integer iproc, ipol
      logical oprint, rodft, dftmp2, status, converged, ignore
      double precision cputime, walltime
c
      logical bsemol, task_energy
      external bsemol, task_energy
c
c      call ecce_print_module_entry('task gwmol')
c
c     Initialization 
c
      task_bsemol = .false.
c
c     Time at the start
c
      cputime  = util_cpusec()
      walltime = util_wallsec()
c
c     Invalidate status
c
      if (.not. rtdb_put( rtdb, 'task:status', mt_log, 1, .false.))
     &  call errquit('task_bse: failed to invalidate status',0,RTDB_ERR)
c
c     Print BSE header
c
      oprint = util_print('information', print_low)
      if (oprint) then
        if (.not. rtdb_cget(rtdb, 'title', 1,  title)) title = ' '
        if (ga_nodeid().eq.0) then
          write(LuOut,*)
          call util_print_centered(
     &      LuOut, 'NWChem BSE Module', 40, .true.)
          write(LuOut,*)
          write(LuOut,*)
          if (title.ne.' ') then
            call util_print_centered(LuOut, title, 40, .false.)
            write(LuOut,*)
            write(LuOut,*)
          endif
        endif
      endif
c
c     Check if theory has been specified
c
      if (.not. rtdb_cget(rtdb, 'task:theory', 1, theory))
     &  call errquit('task_bsemol: theory not specified',0, RTDB_ERR)
c
c     Restrict to DFT for now (SCF could also be used)
c
      if (.not.inp_compare(.false.,'dft',theory))
     &  call errquit('task_bsemol: BSE is only implemented for DFT', 0,
     &                INPUT_ERR)
c
      if (rtdb_get(rtdb, 'dft:rodft', mt_log, 1, rodft))
     &  call errquit('task_bsemol: BSE is not implemented for RODFT', 0,
     &                INPUT_ERR)
c
      if (rtdb_get(rtdb, 'dft:dftmp2', mt_log, 1, dftmp2))
     &  call errquit('task_bsemol: BSE is not implemented for double'//
     &               ' hybrid functionals',0,INPUT_ERR)
c
      if (.not.rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &  call errquit('task_bsemol: failed to read ipol', 0, RTDB_ERR)

      if (.not. rtdb_get(rtdb, 'task:ignore', mt_log, 1, ignore))
     &   ignore = .true. 
c
c     Make sure a DFT solution is available
c
      status = rtdb_get(rtdb, 'dft:converged', mt_log, 1, converged)
      if (.not.status) then
        status = task_energy(rtdb)
        if (.not.rtdb_get(rtdb, 'dft:converged', mt_log, 1, converged))
     &    call errquit('task_gwmol: failed to read convergence info',0,
     &                  RTDB_ERR)   
      endif

      status = converged


      if (ga_nodeid().eq.0) then
        write(luout,*) " "
        write(luout,'(29X,A)') "***********************"
        write(luout,'(29X,A)') "***       BSE       ***"
        write(luout,'(29X,A)') "***********************"
        write(luout,*) " "

        if (.not.converged) then
          write(LuOut,*) 
          write(LuOut,'(27X,A)') '!!!       WARNING:       !!!'
          write(LuOut,'(27X,A)') '!!! DFT DID NOT CONVERGE !!!'
          write(LuOut,*) 
        endif

        call util_flush(luout)
      endif

      if ((.not.status) .and. (.not.ignore)) 
     &  call errquit('task_bsemol: Energy failure',0,CALC_ERR)
c
c     Actually compute GW 
c
      status = bsemol(rtdb)
c
      cputime  = util_cpusec() - cputime
      walltime = util_wallsec() - walltime
c
       if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     &     call errquit('task_bse: failed to set status',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cputime))
     &     call errquit('task_bse: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, walltime))
     &     call errquit('task_bse: failed storing walltime',0,
     &       RTDB_ERR)
c
      task_bsemol = status
c
      end function
