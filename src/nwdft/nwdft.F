      logical function nwdft(rtdb)
c     
c>>>  driver 
c     
      implicit real*8 (a-h,o-z)
      logical LResult
      logical ok,valo
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
#include "handles.fh"
#include "mafdecls.h"
#include "global.fh"
#include "tcgmsg.fh"
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
c     
c-->  Initialize I/O logical units.
c     
      if( iproc.eq.0 ) call stdio_init
c
c     Push context down to DFT (copied from scf.F)
c
      if (.not. context_push('dft'))
     $     call errquit('dft: context_push failed', 0)
c
c**** status = rtdb_parallel(.true.) ! Broadcast reads to all processes
c
c-->  Read input deck to determine if there are any flags requiring
c     the checkpoint file to be opened as old.
c
      if( iproc.eq.0 ) call makechk
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
c-->
c     
      call scfdriver(rtdb,  geom, AO_bas_han, 
     &               CD_bas_han, XC_bas_han, Int_MB(ickcon), 200)
c
c--> Clean-up memory.
c
      LResult = MA_Pop_Stack(lckcon) ! moved here (HLT)
c
c--> Flush ouput buffer.
c
      call flush(lwr)
c
c--> Print memory usage summary.
c
      if( iproc.eq.0 ) call MA_Summarize_Allocated_Blocks()
      if( iproc.eq.0 ) call MA_Print_stats()
c     
c-->  Clean-up parallel environment.
c     
      call ga_terminate()
      call pend

      nwdft = .true.
      return
      end
 
