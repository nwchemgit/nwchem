      subroutine util_allocga()
c $Id$
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer ncrap,g_crapa,g_crapb,
     A     mappa(8),ii
      integer global_b,global
      logical status
chack
chack one big ga_create to bypass shmem/armci bug 
chack for asymmetric armci memory allocation
chack 
c
      global_b = ga_memory_avail()
      If ( Global_B .ge. 0) then
         global = MA_SizeOf(MT_Byte, global_b, MT_Dbl)
      Else
         global = 0
      EndIf
      ncrap=sqrt(0.85d0*(global*ga_nnodes()))
#if 0
       if(ga_nodeid().eq.0)
     W      write(0,*) ga_nodeid(),' ncrap ',ncrap,global
       call ga_sync()
#endif
       mappa(1)=1
       do ii=2,ga_nnodes()
          mappa(ii)=(ii-1)*ncrap/ga_nnodes()+1
       enddo
      status =  ga_create_irreg(MT_DBL, ncrap, ncrap,'da',
     $     mappa, ga_nnodes(), 
     $     mappa, 1, g_crapa)
       if (.not. ga_destroy(g_crapa)) call errquit
     E      (' gadestroy failed ',0,0)
#if 0
       if(ga_nodeid().eq.0) write(0,*) ' gadestr ok'
#endif
       return
       end
