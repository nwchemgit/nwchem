*      logical function sodft_aoints(rtdb)
*      implicit none
*#include "rtdb.fh"
*      integer rtdb
*c
*      logical nwdft
*      external nwdft
*      logical status
*c
*cgk deub
*      write(*,*)'gk: in sodft_aoints in (nwdft.F) calling nwdft'
*cgk end
*      status = rtdb_cput(rtdb,'dft:theory', 1, 'sodft')
*      sodft_energy = nwdft(rtdb)
*cgk deub
*      write(*,*)'gk: in sodft_aoints back from nwdft'
*cgk end
*      call grid_cleanup(.true.)
*c
*      end
      logical function sodft_aoints(rtdb)
c     
c>>>  driver 
c     
C$Id: nwdft.F,v 1.76 2008/11/11 00:40:27 niri Exp $
      implicit none
#include "errquit.fh"
c****
c**** nwchem handles
c****
      integer rtdb
c     
c-->  Include standard input and memory checking common structures.
c     
c     
c>>>  original NWdft includes ...
c     
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "stdio.fh"
c     
c>>>  includes added on addition of NWdft to nwchem
c     
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "inp.fh"
#include "cdft.fh"
#include "util.fh"
#include "cgridfile.fh"
#include "cosmo.fh"
c     
c     local declarations
c
      double precision energy
      integer iproc,itype
      integer bases(3), i3
      logical LResult, oprint, converged
      logical dyall_mod_dir  ! do modified Dirac calculation
c
      logical wrt_dft_aoints
      external wrt_dft_aoints
      logical nadct_trace
      external nadct_trace
      logical dft_main0d, movecs_converged,grid_reopen,xc_gotxc  
      external dft_main0d, movecs_converged,grid_reopen,xc_gotxc  
      logical grid_ok,l1ecache
      integer igok
      character*80 theory
      character*30 operation
c     
cgk debug
*     write(*,*)'gk: entered sodft_aoints...'
*     write(*,*)'gk:    do some hand shaking so we can transfer some of'
*    & ,' this to aoints'
cgk end
      sodft_aoints = .false.
      iproc  = ga_nodeid()
c
c     Make sure database is parallel, push context down to DFT.
c
      LResult = rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(rtdb, 'dft')

      if (.not. rtdb_cget(rtdb, 'task:operation', 1, operation))
     $     operation = ' '

      if(.not.rtdb_cput(rtdb,'dft:theory', 1, 'sodft')) 
     &  call errquit('dft_aoints:  can not put theory on rtdb',
     &   0,RTDB_ERR)
c
c     init pstat
c	
cgk debug
*      write(*,*)'gk: calling dft_pstat_init from sodft_aoints'
cgk end
      call dft_pstat_init(rtdb)
cgk debug
*     write(*,*)'gk: back from dft_pstat_init from sodft_aoints'
cgk end
c
c     init fdist
c	
cgk debug
*     write(*,*)'gk: calling dft_fdist_init from sodft_aoints'
cgk end
      call dft_fdist_init(rtdb)
cgk debug
*     write(*,*)'gk: back from dft_fdist_init from sodft_aoints'
cgk end
c
c     Check for aoints...probably not neccessary now.
c
*     if (.not. rtdb_cget(rtdb, 'task:operation',1,operation))
*    $       call errquit('sodft_aoints: task:operation is not in rtdb',
*    &        0,RTDB_ERR)
c
      oprint = util_print('information', print_low)
      if (oprint) then
        if (.not. rtdb_cget(rtdb, 'title', 1, title))
     &       title = ' '
        if(iproc.eq.0) then
          call util_print_centered(
     &         LuOut, 'NWChem DFT Module', 40, .true.)
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
c     Gather input from all available resources e.g. input unit,
c     old checkpoint files, etc.
c
cgk debug
*     write(*,*)'gk: calling dft_rdinput from sodft_aoints'
cgk end
      call dft_rdinput(rtdb)
cgk debug
*     write(*,*)'gk: back from dft_rdinput from sodft_aoints'
cgk end
c
c     ----- cosmo initialization ----
c
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
c     Can we check here to see if system if already converged?
c
*     if (movecs_in .eq. movecs_out) then
*        if (.not. rtdb_cget(rtdb, 'dft:theory', 1, theory))
*    $        call errquit('dft: theory not specified',0, RTDB_ERR)
*        if(theory .eq. 'dft'.or.theory.eq.'sodft')then 
*           if(.not.xc_gotxc()) then
*              grid_ok=.true.
*           else
*              grid_ok=grid_reopen(geom)
*           endif
*           igok=0
*           if(grid_ok) igok=1
*           call ga_igop(999,igok,1,'+')
*           grid_ok=igok.eq.ga_nnodes()
*           if(.not.grid_ok) call grid_cleanup(.true.)
*           if (movecs_converged(rtdb, ao_bas_han,theory, movecs_in).
*    &      and.grid_ok)     then
c
c     check if we can get the grid from a file
c         
*              if(theory.eq.'sodft') then
*              if (.not. rtdb_get(rtdb, 'sodft:energy', mt_dbl, 1, 
*    &              energy)) call errquit(
*    .             'sodft_aoints: failed getting converged energy',0, 
*    .              RTDB_ERR)
*              else
*              if (.not. rtdb_get(rtdb, 'dft:energy', mt_dbl, 1, 
*    &              energy)) call errquit(
*    .              'sodft_aoints: failed getting converged energy',0,
*    &       RTDB_ERR)
*              endif
*              if (iproc.eq.0 .and. oprint)then
*                 write(LuOut,1011)energy
*                 call util_flush(Luout)
*              endif
*              converged = .true.
*              sodft_aoints = .true.
*              goto 1101
*           endif
*        endif
*     endif
*1011 format(/'  The DFT is already converged '//,
*    &     '         Total DFT energy =', f20.12/)
c     
c     Initialize the integral object
c    
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
c
cgk debug
*     write(*,*)'gk: !!!!!!!!!!!!!!calling int_init from sodft_aoints'
*     write(*,*)'gk:  i3 = ', i3
cgk end
      call int_init(rtdb, i3, bases)
cgk debug
*     write(*,*)'gk: !!!!!!!!!!!!!!back from int_init from sodft_aoints'
cgk end
c
c     get 1-el in cache
c
cgk debug
*     write(*,*)'gk: !!!!!!!!!!!calling schwartz_init from sodft_aoints'
cgk end
c     do we want Schwarz screening?
      call schwarz_init(geom, AO_bas_han)
cgk debug
*     write(*,*)'gk: !!!!!!!!!back from schwartz_init from sodft_aoints'
cgk end
c
      l1ecache=.false.
*     if ((nbf_ao*nbf_ao*3*8)/ga_nnodes() .lt. 8000000) then! Cache if < 8MB/node
*        if (iproc.eq.0 .and. oprint)then
*           write(LuOut,*) ' Caching 1-el integrals '
*           call util_flush(Luout)
*        endif
*        l1ecache=.true.
*        call int_1e_cache_ga(AO_bas_han, oskel)
*     endif
c
c     Build the grid.
c
      grid_written=.false.
cgk debug
*     write(*,*)'gk: calling get_setquad in sodft_aoints'
*     write(*,*)'gk:    building the grid'
cgk end
      call grid_setquad(rtdb)
cgk debug
*     write(*,*)'gk: back from get_setquad in sodft_aoints'
cgk end
c
c     Analyze set parameters for consistency
c
cgk debug
*     write(*,*)'gk: calling dft_inpana in sodft_aoints'
*     write(*,*)'gk:    analysing dft parameters '
cgk end
      call dft_inpana(rtdb)
cgk debug
*     write(*,*)'gk: back from dft_inpana in sodft_aoints'
cgk end
c
c     
c--> determine system type (molecule, polymer,surface,crystal)
c    if no system type exists in db --> default = molecular case     
c
      if(.not. geom_systype_get(geom,itype))
     &    itype = 0
c
      if(itype.eq.0) then
#ifdef GA_TRACE
      call trace_init(100000)                                 ! initialize trace
#endif
        if (inp_compare(.false., 'aoints', operation)) then
cgk debug
*     write(*,*)'gk: calling wrt_dft_aoints from sodft_aoints'
cgk end
         sodft_aoints = wrt_dft_aoints(rtdb)
cgk debug
*     write(*,*)'gk: back from from sodft_aoints'
cgk end
        else if(inp_compare(.false., 'grdtrc', operation)) then
           call int_terminate
           call intd_init(rtdb, i3, bases)
           sodft_aoints = nadct_trace(rtdb,.false.,.false.)
        else if(inp_compare(.false., 'nadtrc', operation)) then
           call int_terminate
           call intd_init(rtdb, i3, bases)
           sodft_aoints = nadct_trace(rtdb,.true.,.true.)
#ifdef GA_TRACE
      call trace_end(ga_nodeid())                            ! end trace
#endif
        else
          call errquit('sodft_aoints: task not defined',0, 
     &     INPUT_ERR)
        endif
      else
        call errquit('sodft_aoints: problem with system type',0, 
     &   INPUT_ERR)
      endif
c     
c     terminate integrals 
c     
cgk debug
*      write(*,*)'gk: starting dft cleanup in sodft_aoints'
cgk end
      if(l1ecache) call int_1e_uncache_ga()
cgk doing this in nadctrc, etc.
c      call schwarz_tidy
c      call int_terminate
cgk end
c
c     Go here if already converged
c
 1101 continue
c
c     Clean up cosmo stuff if necessary
c
      if (cosmo_on) then
         call cosmo_tidy()
         cosmo_on = .false.
         cosmo_phase = 1
      endif
c     
c     destroy BS handles and geom object
c
      XCFIT    = .TRUE.
      if( XC_bas_han.eq.-99999 ) XCFIT = .FALSE.
      CDFIT    = .TRUE.
      if( CD_bas_han.eq.-99999 ) CDFIT = .FALSE.

      if(XCFIT)then
       if(.not.bas_destroy(XC_bas_han))
     &     call errquit ('dft:destroying XC basis handle',0, BASIS_ERR)
      endif
      if(CDFIT)then
       if(.not.bas_destroy(CD_bas_han))
     &     call errquit ('dft:destroying CD basis handle',0, BASIS_ERR)
      endif
      if (.not.(
     &     bas_destroy(AO_bas_han)
     &     .and. geom_destroy(geom)))
     &     call errquit ('dft:destroying geom and basis handles',0,
     &       BASIS_ERR)
c
c--> Flush ouput buffer.
c
      if(iproc.eq.0) call util_flush(LuOut)

      call dft_pstat_print
      call util_print_pop
cgk debug
*     write(*,*)'gk: leaving sodft_aoints'
cgk end
c
      return
c
      end
 
