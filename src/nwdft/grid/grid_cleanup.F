      Subroutine grid_cleanup(okeep)
c
c$Id$
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
      call grid_closecache
c
      if(.not.xc_gotxc().or.(.not.file_opened)) return
c
      if(.not.largenode) then
         if(.not.grid_file_close(okeep)) call
     .        errquit(' grclup: file_close failed',0, DISK_ERR)
      endif
      file_opened=.false.

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
#include "util.fh"
      integer rtdb
      logical use_nwxc, lhashfxc
      grid_clinit=rtdb_get(rtdb, 'dft:xfac', mt_dbl, numfunc,xfac)
      use_nwxc = util_module_avail("nwxc")
      if (use_nwxc) then
        call nwxc_rtdb_load(rtdb,"dft",use_nwxc)
      endif
      if (use_nwxc) then
        call nwxc_getvals("nwxc_is_on",use_nwxc)
      endif
      if (use_nwxc) then
        call nwxc_getvals("nwxc_has_hfx",lhashfxc)
        if (lhashfxc) then
          call nwxc_getwght("nwxc_wght_hfx",xfac(1))
        endif
      endif
      grid_clinit = grid_clinit.or.use_nwxc
      call util_file_name('grinfo',.true.,.true.,grinfo)
      return
      end
