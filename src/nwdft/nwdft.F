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
c>>>  begin include "handles.fh"
#include "handles.fh"
c>>>  end include "handles.fh"
c>>>  begin include "mafdecls.h"
#include "mafdecls.h"
c>>>  end include "mafdecls.h"
c>>>  begin include "global.fh"
#include "global.fh"
c>>>  end include "global.fh"
c>>>  begin include "tcgmsg.fh"
#include "tcgmsg.fh"
c>>>  end include "tcgmsg.fh"
c     
c>>>  includes added on addition of NWdft to nwchem
c     
c>>>  begin include "bas_dft.fh"
cremoved "bas.fh"
#include "bas.fh"
c>>>  end include "bas.fh"
c>>>  begin include "geom.fh"
#include "geom.fh"
c>>>  end include "geom.fh"
c>>>  begin include "rtdb.fh"
#include "rtdb.fh"
c>>>  end include "rtdb.fh"
c>>>  begin include "context.fh"
#include "context.fh"
c>>>  end include "context.fh"
c     
c>>>  
c     
      integer ienv(4)

      iproc  = nodeid()
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
c-->  Initialize memory.  
c
c**** if( iproc.eq.0 ) call rdmem(Mstack,Mheap)
c****
c****      msg = 10000
c****      call brdcst(msg,Mstack,mitob(1),0)
c****      msg = msg + 1
c****      call brdcst(msg,Mheap,mitob(1),0)

c
c>>>  initialized memory allocation in stand-alone nwdft
c>>>  not needed in nwchem version
c>>>  call getmem(Mstack,Mheap)
c
c**** valo=.true.
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
      write (6,*) 'nwdft: calling scfdriver'
      call scfdriver(rtdb,  geom, AO_bas_han, 
     &               CD_bas_han, XC_bas_han, Int_MB(ickcon), 200)
c
c--> Clean-up memory.
c
c****
c**** these seem to be dead code  (HT)
c****
c**** if(Int_MB(ickcon+43).eq.3)then
c**** LResult = MA_Pop_Stack(lfxyz)
c**** LResult = MA_Pop_Stack(lgmat)
c**** LResult = MA_Pop_Stack(lamat)
c**** endif
c**** LResult = MA_Pop_Stack(lcharge) ! now deallocated in rdinput
c**** if(iproc.eq.0) then
c****   LResult = MA_Pop_Stack(lcoord) ! now deallocated in rdinput
c       LResult = MA_Pop_Stack(latmtyp) ! old array removed from sym
c****   LResult = MA_Pop_Stack(lckcon)
c**** endif
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
 
