      logical function dft_energy(rtdb)
      implicit none
#include "rtdb.fh"
      integer rtdb
c
      logical nwdft
      external nwdft
      logical status
c
      status = rtdb_cput(rtdb,'dft:theory', 1, 'dft')
      dft_energy = nwdft(rtdb)
      call grid_cleanup(rtdb,.true.)
c
      end
      logical function sodft_energy(rtdb)
      implicit none
#include "rtdb.fh"
      integer rtdb
c
      logical nwdft
      external nwdft
      logical status
c
      status = rtdb_cput(rtdb,'dft:theory', 1, 'sodft')
      sodft_energy = nwdft(rtdb)
      call grid_cleanup(rtdb,.true.)
c
      end
      logical function nwdft(rtdb)
c     
c>>>  driver 
c     
C$Id: nwdft.F,v 1.52 2001-01-12 00:18:00 edo Exp $
      implicit none
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
#include "cdft.fh"
#include "util.fh"
#include "cgridfile.fh"
c     
c     local declarations
c
#ifdef DFT_TIMER
      character*(nw_max_path_len) timings
#endif
      double precision energy
      integer iproc,itype
      integer bases(3), i3
      logical LResult, oprint, converged
      logical dyall_mod_dir  ! do modified Dirac calculation
c
      logical dft_main0d, movecs_converged
      external dft_main0d, movecs_converged
      logical ocosmo
      character*80 theory
c     
      nwdft = .false.
      iproc  = ga_nodeid()
c
c     If Dyall-mod-dir is set, exit with error
c
      if (rtdb_get(rtdb,'dyall_mod_dir',mt_log,1,dyall_mod_dir)) then
        if (dyall_mod_dir) call errquit(
     &      'Cannot do DFT calculation with Dyall-modified-Dirac method'
     &      ,99)
      end if
c
c     Make sure database is parallel, push context down to DFT.
c
      LResult = rtdb_parallel(.true.)
      call util_print_push
      call util_print_rtdb_load(rtdb, 'dft')
#ifdef DFT_TIMER
c     open timings file as unit 7 in the perm directory
c     only process 0 will write to it
      if(iproc.eq.0) then
c      call util_file_name('dftim',.true.,.false.,timings)
      call util_file_name('dftim',.false.,.false.,timings)
      open(7,file=timings,form='formatted',status='unknown',err=1010)
      endif
#endif
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
      call dft_rdinput(rtdb)
c
c     bail out if cosmo active and spherical bfns
c
      if ( rtdb_get(rtdb,'slv:cosmo',mt_log,1,ocosmo)) then
         if(ocosmo.and.bas_is_spherical(ao_bas_han)) call errquit(
     .        ' cosmo does not handle spherical basis yet',0)
      endif
c
c     Build the grid.
c
      grid_written=.false.
      call grid_setquad(rtdb)
c
c     Analyze set parameters for consistency
c
      call dft_inpana(rtdb)
c
c     Can we check here to see if system if already converged?
c
      if (movecs_in .eq. movecs_out) then
         if (.not. rtdb_cget(rtdb, 'dft:theory', 1, theory))
     $        call errquit('dft: theory not specified',0)
         if(theory .eq. 'dft')then 
            if (movecs_converged(rtdb, ao_bas_han, 'dft', movecs_in))
     &           then
               if (.not. rtdb_get(rtdb, 'dft:energy', mt_dbl, 1, 
     &              energy)) call errquit(
     .              'dft_scf: failed getting converged energy',0)
               if (ga_nodeid().eq.0 .and. oprint)then
                  write(LuOut,1011)energy
                  call util_flush(6)
               endif
               converged = .true.
               nwdft = .true.
               goto 1101
            endif
         else if(theory .eq. 'sodft')then 
            if (movecs_converged(rtdb, ao_bas_han, 'sodft', movecs_in))
     &           then
               if (.not. rtdb_get(rtdb, 'sodft:energy', mt_dbl, 1, 
     &              energy)) call errquit(
     .             'dft_scf: failed getting converged energy',0)
               if (ga_nodeid().eq.0 .and. oprint)then
                  write(LuOut,1011)energy
                  call util_flush(6)
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
      if (i3 .gt. 1) call int_app_set_no_texas(rtdb)
      call int_init(rtdb, i3, bases)
c     
c--> determine system type (molecule, polymer,surface,crystal)
c    if no system type exists in db --> default = molecular case     
c
      if(.not. geom_systype_get(geom,itype))
     &    itype = 0
c
      if(itype.eq.0) then
        nwdft = dft_main0d(rtdb)
      else
        call errquit('nwdft: problem with system type',0)
      endif
c     
c     terminate integrals 
c     
      call int_terminate
      if (i3 .gt. 1) call int_app_unset_no_texas(rtdb)
c
c     Go here if already converged
c
 1101 continue
c     
c     destroy BS handles and geom object
c
      XCFIT    = .TRUE.
      if( XC_bas_han.eq.-99999 ) XCFIT = .FALSE.
      CDFIT    = .TRUE.
      if( CD_bas_han.eq.-99999 ) CDFIT = .FALSE.

      if(XCFIT)then
       if(.not.bas_destroy(XC_bas_han))
     &     call errquit ('dft:destroying XC basis handle',0)
      endif
      if(CDFIT)then
       if(.not.bas_destroy(CD_bas_han))
     &     call errquit ('dft:destroying CD basis handle',0)
      endif
      if (.not.(
     &     bas_destroy(AO_bas_han)
     &     .and. geom_destroy(geom)))
     &     call errquit ('dft:destroying geom and basis handles',0)
c
c--> Flush ouput buffer.
c
      if(iproc.eq.0) call util_flush(LuOut)

      call util_print_pop
c
#ifdef DFT_TIMER
      close(7)
#endif
      return
c
#ifdef DFT_TIMER
 1010 call errquit('dft: failed to open timings file',0)
#endif
      end
 
