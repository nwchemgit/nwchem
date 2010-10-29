      subroutine grid_read_ck(rtdb)
c
c$Id$
c
      implicit none
#include "errquit.fh"
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
      double precision rad
c
      if (.not.MA_Push_get(MT_dbl, buffer_size, 
     &                     'quad pts buffer', lrqbuf, irqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 0,
     &       MA_ERR)
c
c     rewind grid pts file
c
      if (.not. grid_file_rewind())
     $   call errquit('grid_setspac: rewinding gridpts?', 0,
     &       UNKNOWN_ERR)
c
      nqpts_per_buf = n_per_rec
c
      do 100 iqt = 1, n_tot_tasks
c
         call grid_file_read(nqpts_per_buf, nqpts, ictr, rad,
     ,        dbl_mb(irqbuf),nsubb)
          if(nqpts.eq.0) goto 100
c         write(6,*)' ictr, nqpts: ', ictr, nqpts
c
  100 continue
c
      if (.not. MA_Pop_Stack(lrqbuf))
     &   call errquit('grid_setspac: pop stack failed.',0, MA_ERR)
      return
      end 
