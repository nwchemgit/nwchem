      subroutine grid_read_ck(rtdb)
c
c$Id: grid_read_ck.F,v 1.1 1999-04-29 16:36:57 d3h449 Exp $
c
      implicit none
      integer rtdb
c
#include "cdft.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "cgridfile.fh"
c
      integer iqt, ictr
      integer nqpts
      integer lrqbuf, irqbuf
      integer nqpts_per_buf
      logical grid_file_rewind
      external grid_file_rewind
c
      if (.not.MA_Push_get(MT_dbl, buffer_size, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 0)
c
c     rewind grid pts file
c
      if (.not. grid_file_rewind())
     $   call errquit('grid_setspac: rewinding gridpts?', 0)
c
      nqpts_per_buf = n_per_rec
c
      do 100 iqt = 1, n_tot_tasks
c
         call grid_file_read(nqpts_per_buf, nqpts, ictr, dbl_mb(irqbuf))
         write(6,*)' ictr, nqpts: ', ictr, nqpts
c
  100 continue
c
      if (.not. MA_Pop_Stack(lrqbuf))
     &   call errquit('grid_setspac: pop stack failed.',0)
      return
      end 
