      Subroutine grid_cleanup(rtdb)
c
c$Id: grid_cleanup.F,v 1.4 2000-07-29 01:11:57 edo Exp $
c
      implicit none
c      
#include "rtdb.fh"
#include "cdft.fh"
#include "mafdecls.fh"
#include "geom.fh"
c
      integer rtdb
      logical grid_on_file,grid_file_close
c
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &   grid_on_file))then
         grid_on_file = .false.
      endif
      if (grid_on_file)then
c
         if(.not.grid_file_close(.true.)) call
     .        errquit(' grclup: file_close failed',0)
c
!         if (.not.MA_Free_Heap(lrchi))
!     &      call errquit('grid_cleanup: cannot free heap?',0)
!         if (.not.MA_Free_Heap(lrdelchi))
!     &      call errquit('grid_cleanup: cannot free heap?',0)
c
      endif
      return
      end
