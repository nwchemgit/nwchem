c
c     == DFT energy ==
      logical function dft_energy(rtdb)
c
      implicit none
c
#include "rtdb.fh"
c
      integer rtdb
c
      logical nwdft
      external nwdft
      logical status
c
c     == Add theory flag ==
      status = rtdb_cput(rtdb,'dft:theory', 1, 'dft')
c
      dft_energy = nwdft(rtdb)
      call grid_cleanup(.true.)
c
      end
c
c     == Spin-Orbit DFT energy ==
      logical function sodft_energy(rtdb)
c
      implicit none
c
#include "rtdb.fh"
c
      integer rtdb
c
      logical nwdft
      external nwdft
      logical status
c
c     == Add theory flag ==
      status = rtdb_cput(rtdb,'dft:theory', 1, 'sodft')
c
      sodft_energy = nwdft(rtdb)
      call grid_cleanup(.true.)
c
      end
c
c     == DFT RT_TDDFT ==
      logical function dft_rt_tddft(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      logical nwdft
      external nwdft
      logical status
c
c     == Add theory flag ==
      status = rtdb_cput(rtdb,'dft:theory', 1, 'dft')
c
c     == Add real-time TDDFT dynamics flag ==
      if (.not. rtdb_put(rtdb,'dft:lrt_tddft',mt_log,1,.true.))
     &    call errquit('dft_rt_tddft: rtdb_put failed',0,RTDB_ERR)
c
      dft_rt_tddft = nwdft(rtdb)
      call grid_cleanup(.true.)
c
      end
c
c     == Spin-Orbit RT_TDDFT ==
      logical function sodft_rt_tddft(rtdb)
c
      implicit none
c
#include "errquit.fh"
#include "mafdecls.fh"
#include "rtdb.fh"
c
      integer rtdb
c
      logical nwdft
      external nwdft
      logical status
c
c     == Add theory flag ==
      status = rtdb_cput(rtdb,'dft:theory', 1, 'sodft')
c
c     == Add real-time TDDFT dynamics flag ==
      if (.not. rtdb_put(rtdb,'dft:lrt_tddft',mt_log,1,.true.))
     &    call errquit('sodft_rt_tddft: rtdb_put failed',0,RTDB_ERR)
c
      sodft_rt_tddft = nwdft(rtdb)
      call grid_cleanup(.true.)
c
      end
c
c     == Main DFT driver ==
      logical function nwdft(rtdb)
c     
C$Id: nwdft.F,v 1.76 2008-11-11 00:40:27 niri Exp $
c
      implicit none
c     
c     == Include standard input and memory checking common structures ==
c     
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "util.fh"
#include "cgridfile.fh"
#include "cosmo.fh"
c
      integer rtdb
c     
c     == Local declarations ==
      double precision energy
      integer iproc,itype
      integer bases(3), i3
      logical LResult, oprint, converged
      logical dyall_mod_dir  ! do modified Dirac calculation
c
      logical dft_main0d, movecs_converged,grid_reopen,xc_gotxc  
      external dft_main0d, movecs_converged,grid_reopen,xc_gotxc  
      logical grid_ok,l1ecache
      integer igok
      character*80 theory
      logical lrt_tddft
      character*32 pname
c     
      pname = "nwdft:"
      nwdft = .false.
      iproc  = ga_nodeid()
c
c     == Make sure database is parallel, push context down to DFT ==
      LResult = rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(rtdb, 'dft')
c
c     == Init pstat ==
      call dft_pstat_init(rtdb)
c
c     == Init fdist ==
      call dft_fdist_init(rtdb)
c
c     == If Dyall-mod-dir is set, exit with error ==
      if (rtdb_get(rtdb,'dyall_mod_dir',mt_log,1,dyall_mod_dir)) then
        if (dyall_mod_dir) call errquit(pname//
     &      'Cannot do DFT calculation with Dyall-modified-Dirac method'
     &      ,99, RTDB_ERR)
      end if
c
c     == Print header ==
      oprint = util_print('information', print_low)
      lrt_tddft = .false.
      if (.not. rtdb_get(rtdb,'dft:lrt_tddft',mt_log,1,lrt_tddft))
     &   lrt_tddft = .false.
c
      if (oprint) then
        if (.not. rtdb_cget(rtdb, 'title', 1, title))
     &       title = ' '
        if(iproc.eq.0) then
          if (.not.lrt_tddft) then
            write(LuOut,*)
            call util_print_centered(
     &         LuOut, 'NWChem DFT Module', 40, .true.)
          else
            write(LuOut,*)
            call util_print_centered(
     &         LuOut, 'NWChem Real-Time TDDFT Module', 40, .true.)
          end if
          write(LuOut,*)
          write(LuOut,*)
          if (title .ne. ' ') then
            call util_print_centered(LuOut, title, 40, .false.)
            write(LuOut,*)
            write(LuOut,*)
          endif
        endif
      endif
c
c     == Gather input from all available resources e.g. input unit, ==
c     == old checkpoint files, etc. ==
      call dft_rdinput(rtdb)
c
c     == Cosmo initialization ==
      cosmo_last = .false.
      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,cosmo_on)) then
         if(cosmo_on) then
            call cosmo_initialize(rtdb,geom,ao_bas_han,oprint)
c
c           Turn cosmo on, we want to run the calculation
c           Start with gas_phase run
c
            cosmo_last = .true.
            cosmo_on = .true.
            if(.not.rtdb_get(rtdb,'cosmo_phase',mt_int,1,cosmo_phase))
     >         cosmo_phase = 1

         endif
      endif
c
c     == Check to see if system is already converged? ==
      if (movecs_in .eq. movecs_out) then
         if (.not. rtdb_cget(rtdb, 'dft:theory', 1, theory))
     $        call errquit(pname//'Theory not specified',0, RTDB_ERR)
c
          if(iproc.eq.0) then
           if (theory .eq. 'sodft') then
             call util_print_centered(
     &         LuOut, 'Spin-Orbit DFT', 40, .true.)
             write(LuOut,*)
             write(LuOut,*)
           endif
          endif
c
         if(theory .eq. 'dft'.or.theory.eq.'sodft')then 
            if(.not.xc_gotxc()) then
               grid_ok=.true.
            else
               grid_ok=grid_reopen(geom)
            endif
            igok=0
            if(grid_ok) igok=1
            call ga_igop(999,igok,1,'+')
            grid_ok=igok.eq.ga_nnodes()
            if(.not.grid_ok) call grid_cleanup(.true.)
            if (movecs_converged(rtdb, ao_bas_han,theory, movecs_in).
     &      and.grid_ok)     then
c
c     check if we can get the grid from a file
c         
               if(theory.eq.'sodft') then
               if (.not. rtdb_get(rtdb,'sodft:energy',mt_dbl,1,energy)) 
     &           call errquit(
     &           pname//'Failed getting converged energy',0, RTDB_ERR)
               else
               if (.not. rtdb_get(rtdb,'dft:energy', mt_dbl,1,energy)) 
     &           call errquit(
     .           pname//'Failed getting converged energy',0, RTDB_ERR)
               endif
               if (iproc.eq.0 .and. oprint)then
                  write(LuOut,1011)energy
                  call util_flush(Luout)
               endif
               converged = .true.
               nwdft = .true.
               goto 1101
            endif
         endif
      endif
 1011 format(/'  The DFT is already converged '//,
     &     '         Total DFT energy =', f20.12/)
c     
c     == Initialize the integral object ==
      if(XCFIT.and.CDFIT) then
        i3 = 3
        bases(1) = AO_bas_han
        bases(2) = CD_bas_han
        bases(3) = XC_bas_han
      elseif((.not.XCFIT).and.CDFIT) then
        i3 = 2
        bases(1) = AO_bas_han
        bases(2) = CD_bas_han
      elseif((.not.CDFIT).and.XCFIT) then
        i3 = 2
        bases(1) = AO_bas_han
        bases(2) = XC_bas_han
      else
        i3 = 1
        bases(1) = AO_bas_han
      endif
      call int_init(rtdb, i3, bases)
c
c     == Get 1-el in cache ==
      call schwarz_init(geom, AO_bas_han)
c
      l1ecache=.false.
      if ((nbf_ao*nbf_ao*3*8)/ga_nnodes() .lt. 8000000) then ! Cache if < 8MB/node
         if (iproc.eq.0 .and. oprint)then
            write(LuOut,*) ' Caching 1-el integrals '
            call util_flush(Luout)
         endif
         l1ecache=.true.
         call int_1e_cache_ga(AO_bas_han, oskel)
      endif
c
c     == Build the grid ==
      grid_written=.false.
      call grid_setquad(rtdb)
c
c     == Analyze set parameters for consistency ==
      call dft_inpana(rtdb)
c     
c     == Determine system type (molecule, polymer,surface,crystal) ==
c     == if no system type exists in db --> default = molecular case ==    
      if(.not. geom_systype_get(geom,itype)) itype = 0
c
      if(itype.eq.0) then
#ifdef GA_TRACE
      call trace_init(100000)                                 ! initialize trace
#endif
        nwdft = dft_main0d(rtdb)
#ifdef GA_TRACE
      call trace_end(ga_nodeid())                            ! end trace
#endif
      else
        call errquit(pname//'Problem with system type',0, INPUT_ERR)
      endif
c     
c     == Terminate integrals ==
      if(l1ecache) call int_1e_uncache_ga()
      call schwarz_tidy
      call int_terminate
c
c     == Go here if already converged ==
c
 1101 continue
c
c     == Clean up cosmo stuff, if necessary ==
      if (cosmo_on) then
         call cosmo_tidy()
         cosmo_on = .false.
         cosmo_phase = 1
      endif
c     
c     == Destroy basis set handles and geom object ==
      XCFIT    = .TRUE.
      if( XC_bas_han.eq.-99999 ) XCFIT = .FALSE.
      CDFIT    = .TRUE.
      if( CD_bas_han.eq.-99999 ) CDFIT = .FALSE.

      if(XCFIT)then
       if(.not.bas_destroy(XC_bas_han))
     &   call errquit(pname//'Destroying XC basis handle',0, BASIS_ERR)
      endif
      if(CDFIT)then
       if(.not.bas_destroy(CD_bas_han))
     &  call errquit(pname//'Destroying CD basis handle',0, BASIS_ERR)
      endif
      if (.not.(
     &  bas_destroy(AO_bas_han)
     &  .and. geom_destroy(geom)))
     &  call errquit (pname//'Destroying geom and basis handles',0,
     &       BASIS_ERR)
c
c     == Flush output buffer ==
      if(iproc.eq.0) call util_flush(LuOut)
c
      call dft_pstat_print
      call util_print_pop
c
      return
      end
