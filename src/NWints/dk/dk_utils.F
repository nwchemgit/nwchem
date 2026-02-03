      subroutine ga_serialcopy(g_inpg,g_inworld)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
      integer g_inpg,g_inworld
c
      integer itype,dim
      integer k_a,l_a,i,j
c
      if(ga_nodeid().eq.0) then
         call ga_inquire(g_inpg, itype, dim, dim)
         if (.not. ma_push_get(MT_DBL,dim*dim,'MA temp',l_a,k_a))
     $        call errquit('serialc: ma_push_get failed MA temp',
     Z        0       ,     MA_ERR)
         call ga_get(g_inpg,1,dim,1,dim,dbl_mb(k_a),dim)
         call ga_put(g_inworld,1,dim,1,dim,dbl_mb(k_a),dim)
         if (.not. ma_pop_stack(l_a))
     $        call errquit('serialc: ma_pop_stack l_a',0, MA_ERR)
      endif
      call ga_sync()
      return
      end
