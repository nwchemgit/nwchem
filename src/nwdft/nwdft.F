      logical function nwdft(rtdb)
c     
c>>>  driver 
c     
C$Id: nwdft.F,v 1.13 1995-03-01 23:00:57 og845 Exp $
      implicit real*8 (a-h,o-z)
      logical LResult,converged
      logical ok,CDFIT,XCFIT
      logical status
c****
c**** nwchem handles
c****
      integer rtdb, geom, AO_bas_han, CD_bas_han, XC_bas_han
c     
c-->  Include standard input and memory checking common structures.
c     
c     
c>>>  original NWdft includes ...
c     
#include "mafdecls.h"
#include "global.fh"
#include "tcgmsg.fh"
#include "io.fh"
c     
c>>>  includes added on addition of NWdft to nwchem
c     
#include "bas.fh"
#include "geom.fh"
#include "rtdb.fh"
#include "context.fh"
c     
c>>>  
c     

      iproc  = ga_nodeid()
C
C	assign stdout to unit 6
C
      lwr=6
c
c     Push context down to DFT (copied from scf.F)
c
      if (.not. context_push('dft'))
     $     call errquit('dft: context_push failed', 0)
c
c**** status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c-->  Allocate memory for checkpoint file constants.
c
      LResult = MA_Push_Get(MT_Int,200,'pdft10 const',lckcon,ickcon)
c
c-->  Gather input from all available resources e.g. input unit,
c     old checkpoint files, etc.
c
      call rdinput(rtdb, geom, AO_bas_han, 
     &             CD_bas_han, XC_bas_han, Int_MB(ickcon), 200)
c     
c--> determine system type (molecule, polymer,surface,crystal)
c    if no system type exists in db --> default = molecular case     
c
      if(.not. geom_systype_get(geom,itype))
     &    itype = 0
c
      if(itype.eq.0) then
        call scfdriver(rtdb,  geom, AO_bas_han, 
     &       CD_bas_han, XC_bas_han, Int_MB(ickcon), 200)
      elseif(itype.eq.3) then
        if(.not. geom_ncent(geom,ncent))
     &     call errquit('nwdft: ncent get failed',0)
        call main3d(geom,rtdt,ncent)
      else
        call errquit('nwdft: problem with system type',0)
      endif
c
c--> Clean-up memory.
c
      LResult = MA_Pop_Stack(lckcon) ! moved here (HLT)
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
      if(iproc.eq.0) call flush(lwr)

      nwdft = .true.
      if (.not. context_pop('dft')) 
     &     call errquit('context_pop failed',0)
      close(7)
      return
      end
 
