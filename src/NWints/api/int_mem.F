      subroutine int_mem(max1e, maxg, mscratch_1e, mscratch_2e)
      implicit none
#include "apiP.fh"
c::passed
      integer max1e, maxg, mscratch_1e, mscratch_2e
c
      if (init_int.ne.0) then
        max1e       = isz_1e
        maxg        = isz_2e4c
        mscratch_1e = mem_1e
        mscratch_2e = mem_2e4c
      else
        write(6,*)'int_mem: int_init not called '
        stop ' error '
      endif
      end
      subroutine int_mem_1e(max1e, mscratch_1e)
      implicit none
#include "apiP.fh"
c::passed
      integer max1e, mscratch_1e
c
      if (init_int.ne.0) then
        max1e       = isz_1e
        mscratch_1e = mem_1e
      else
        write(6,*)'int_mem_1e: int_init not called '
        stop ' error '
      endif
      end
      subroutine int_mem_2e4c(maxg, mscratch_2e)
      implicit none
#include "apiP.fh"
c::passed
      integer maxg, mscratch_2e
c
      if (init_int.ne.0) then
        maxg        = isz_2e4c
        mscratch_2e = mem_2e4c
      else
        write(6,*)'int_mem_2e4c: int_init not called '
        stop ' error '
      endif
      end
