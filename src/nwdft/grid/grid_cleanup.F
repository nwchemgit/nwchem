      Subroutine grid_cleanup(rtdb)
c
c$Id: grid_cleanup.F,v 1.2 2000-03-27 19:15:31 edo Exp $
c
      implicit none
c      
#include "rtdb.fh"
#include "cdft.fh"
#include "mafdecls.fh"
c
      integer rtdb
      logical grid_on_file
c
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &   grid_on_file))then
         grid_on_file = .false.
      endif
      if (grid_on_file)then
c
         call grid_file_close(.true.)
c
!         if (.not.MA_Free_Heap(lrchi))
!     &      call errquit('grid_cleanup: cannot free heap?',0)
!         if (.not.MA_Free_Heap(lrdelchi))
!     &      call errquit('grid_cleanup: cannot free heap?',0)
c
      endif
      return
      end
