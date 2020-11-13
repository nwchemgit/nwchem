c
c     Task handler for GW calcualtions
c
      logical function task_gw(rtdb)
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
      logical oprint, rodft, dftmp2, status
      double precision cputime, walltime
c
      logical dft_energy_gw
      external dft_energy_gw
c
c     Initialization 
c
      task_gw = .false.
      iproc = ga_nodeid()
c
c     Time at the start
c
      cputime  = util_cpusec()
      walltime = util_wallsec()
c
      if (.not. rtdb_put( rtdb, 'task:status', mt_log, 1, .false.))
     &  call errquit('task_gw: failed to invalidate status',0,RTDB_ERR)
c
      oprint = util_print('information', print_low)
      if (oprint) then
        if (.not. rtdb_cget(rtdb, 'title', 1,  title)) title = ' '
        if (iproc.eq.0) then
          write(LuOut,*)
          call util_print_centered(
     &      LuOut, 'NWChem GW@DFT Module', 40, .true.)
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
     &  call errquit('task_gw: theory not specified',0, RTDB_ERR)
c
      if (.not.inp_compare(.false.,'dft',theory))
     &  call errquit('task_gw: GW is only implemented for DFT', 0,
     &                INPUT_ERR)
c
      if (rtdb_get(rtdb, 'dft:rodft', mt_log, 1, rodft))
     &  call errquit('task_gw: GW is not implemented for RODFT', 0,
     &                INPUT_ERR)
c
      if (rtdb_get(rtdb, 'dft:dftmp2', mt_log, 1, dftmp2))
     &  call errquit('task_gw: GW is not implemented for double'//
     &               ' hybrid functionals',0,INPUT_ERR)
c
      if (.not.rtdb_get(rtdb, 'dft:ipol', mt_int, 1, ipol))
     &  call errquit('task_gw: failed to read ipol', 0, RTDB_ERR)
c
      if (ipol.ne.1)
     &  call errquit('task_gw: spin-polarized version is not'//
     &               ' available', 0, INPUT_ERR)
c
c     Make sure a DFT solution is available
c
      status = dft_energy_gw(rtdb)
c
      cputime  = util_cpusec() - cputime
      walltime = util_wallsec() - walltime
c
       if (.not. rtdb_put(rtdb, 'task:status', mt_log, 1, .true.))
     &     call errquit('task_gw: failed to set status',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:cputime', mt_dbl, 1, cputime))
     &     call errquit('task_gw: failed storing cputime',0,
     &       RTDB_ERR)
      if (.not. rtdb_put(rtdb, 'task:walltime', mt_dbl, 1, walltime))
     &     call errquit('task_gw: failed storing walltime',0,
     &       RTDB_ERR)
c
      task_gw = status
c
      end function
