      Subroutine grid_cleanup(okeep)
c
c$Id: grid_cleanup.F,v 1.12 2004-04-17 01:45:05 edo Exp $
c
      implicit none
#include "errquit.fh"
c      
#include "cdft.fh"
#include "mafdecls.fh"
#include "geom.fh"
#include "cgridfile.fh"
c
      logical grid_file_close
      logical okeep
      logical xc_gotxc
      external xc_gotxc
c
c     case of pure HF X
c
      if(.not.xc_gotxc()) return
c
c
      call grid_closecache
      if(.not.largenode) then
         if(.not.grid_file_close(okeep)) call
     .        errquit(' grclup: file_close failed',0, DISK_ERR)
         if(.not.okeep) then
            open(unit=gr_unit,file=grinfo)
            close(unit=gr_unit,status='delete')
         endif
      endif

c
      return
      end
      subroutine grid_closecache()
      implicit none
#include "mafdecls.fh"
#include "cgridfile.fh"
#include "errquit.fh"
      if (l_cache .ne. -1) then
         if (.not. ma_free_heap(l_cache)) then
            call MA_summarize_allocated_blocks()
            call errquit('grid_file_close: free heap?',0, MA_ERR)
          endif
      end if
      l_cache=-1
      return
      end
      logical function grid_opencache()
      implicit none
#include "mafdecls.fh"
#include "cgridfile.fh"
#include "errquit.fh"
#include "stdio.fh"
c
      grid_opencache=.true.
      if (l_cache .eq. -1) then
         if (.not. ma_alloc_get(MT_DBL, max_cache_rec*rec_len_in_dbl, 
     $        'grid', l_cache, k_cache)) then
            write(LuOut,*) ' grid_file: Failed to allocate cache '//
     $           '... halving'
            max_cache_rec = max_cache_rec/2
            grid_opencache=.false.
         endif
      endif
      return
      end
      logical function grid_clinit(rtdb)
      implicit none
#include "cdft.fh"
#include "mafdecls.fh"
#include "cgridfile.fh"
#include "rtdb.fh"
      integer rtdb
      grid_clinit=rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc,xfac)
      call util_file_name('grinfo',.true.,.true.,grinfo)
      return
      end
