      block data nxtval_ga_data
*
* $Id$
*
      implicit none
      integer g_x
      common/c_nxtval_ga/g_x
      data g_x/0/
      end
      subroutine nxtval_ga_initialize()
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
      integer g_x
      common/c_nxtval_ga/g_x
c
*ga:1:0
      if (.not. ga_create(MT_INT, 1, 1, 'nxtval_ga', 0, 0, g_x))
     $     call errquit('nxtval_ga: ga_create failed', 1, GA_ERR)
      call ga_zero(g_x)
c
      end
      integer function nxtval_ga(nnodes)
      implicit none
#include "errquit.fh"
#include "global.fh"
      integer nnodes
      integer g_x
      common/c_nxtval_ga/g_x
c
      if (g_x .eq. 0) call errquit('nxtval_ga: not initialized',0,
     &       GA_ERR)
c
      if (nnodes .gt. 0) then
         nxtval_ga = ga_read_inc(g_x, 1, 1, 1)
      else
         call ga_sync
         call ga_zero(g_x)
      endif
c
      end
      subroutine nxtval_ga_terminate()
#include "global.fh"
      integer g_x
      common/c_nxtval_ga/g_x
c
      if (.not. ga_destroy(g_x)) call errquit
     $     ('nxtval_ga: ga destroy failed', 0, GA_ERR)
      g_x = 0
c
      end

