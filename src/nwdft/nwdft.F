      logical function nwdft(rtdb)
c     
c>>>  driver 
c     
C$Id: nwdft.F,v 1.25 1996-08-19 15:19:53 d3g681 Exp $
      implicit none
      logical CDFIT,XCFIT
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
c     
c     local declarations
c
#ifdef DFT_TIMER
      character*(nw_max_path_len) timings
#endif
      character*255 title
      integer iproc,itype,ncent,ickcon,lckcon
      logical LResult,oprint
c>>>  
c     

      iproc  = ga_nodeid()
c
      call util_print_push
      call util_print_rtdb_load(rtdb, 'dft')
#ifdef DFT_TIMER
c     open timings file as unit 7 in the scratch directory
c     only process 0 will write to it
      call util_file_name('dftim',.true.,.false.,timings)
      open(7,file=timings,form='formatted',status='unknown',err=1010)
#endif
      oprint = util_print('information', print_low)
      if (oprint) then
        if (.not. rtdb_cget(rtdb, 'title', 1, title))
     $       title = ' '
        if(iproc.eq.0) then
          call util_print_centered(
     .         LuOut, 'NWCHEM DFT Module', 40, .true.)
          write(LuOut,*)
          write(LuOut,*)
          if (title .ne. ' ') then
            call util_print_centered(6, title, 40, .false.)
            write(LuOut,*)
            write(LuOut,*)
          endif
        endif
      endif
c
c-->  Allocate memory for checkpoint file constants.
c
      if(.not.MA_Push_Get(MT_Int,200,'pdft10 const',lckcon,ickcon))
     &  call errquit('nwdft: cannot allocate pdft10 const',0)
c
c-->  Gather input from all available resources e.g. input unit,
c     old checkpoint files, etc.
c
      call dft_rdinput(rtdb, 
     &             Int_MB(ickcon), 200)
c     
c--> determine system type (molecule, polymer,surface,crystal)
c    if no system type exists in db --> default = molecular case     
c
      if(.not. geom_systype_get(geom,itype))
     &    itype = 0
c
      if(itype.eq.0) then
        call dft_main0d(rtdb, 
     &       Int_MB(ickcon), 200)
cMG      elseif(itype.eq.3) then
cMG        if(.not. geom_ncent(geom,ncent))
cMG     &     call errquit('nwdft: ncent get failed',0)
cMG        call main3d(geom,rtdb,ncent,AO_bas_han,CD_bas_han,XC_bas_han)
      else
        call errquit('nwdft: problem with system type',0)
      endif
c
c--> Clean-up memory.
c
      if(.not.ma_pop_stack(lckcon))
     &  call errquit('nwdft: cannot pop stack',0)
c     
c     terminate integrals 
c     
      call int_terminate
c     
c     destroy BS handles and geom object
c

      XCFIT    = .TRUE.
      if( XC_bas_han.eq.-99999 ) XCFIT = .FALSE.
      CDFIT    = .TRUE.
      if( CD_bas_han.eq.-99999 ) CDFIT = .FALSE.

      if(XCFIT)then
       if(.not.bas_destroy(XC_bas_han))
     $     call errquit ('dft:destroying XC basis handle',0)
      endif
      if(CDFIT)then
       if(.not.bas_destroy(CD_bas_han))
     $     call errquit ('dft:destroying CD basis handle',0)
      endif
      if (.not.(
     &     bas_destroy(AO_bas_han)
     &     .and. geom_destroy(geom)))
     $     call errquit ('dft:destroying geom and basis handles',0)


c
c--> Flush ouput buffer.
c
      if(iproc.eq.0) call flush(LuOut)

      nwdft = .true.
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
 
