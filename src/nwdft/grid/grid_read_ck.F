      subroutine grid_read_ck(rtdb)
c
c$Id: grid_read_ck.F,v 1.4 2000-03-20 21:11:49 edo Exp $
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
      double precision rad
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
         call grid_file_read(nqpts_per_buf, nqpts, ictr, rad,
     ,        dbl_mb(irqbuf),nsubb)
          if(nqpts.eq.0) goto 100
c         write(6,*)' ictr, nqpts: ', ictr, nqpts
c
  100 continue
c
      if (.not. MA_Pop_Stack(lrqbuf))
     &   call errquit('grid_setspac: pop stack failed.',0)
      return
      end 
