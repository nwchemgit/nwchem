      Subroutine grid_update_rtdb(rtdb)
c
c$Id$
c
      implicit none
#include "errquit.fh"
c      
#include "mafdecls.fh"
#include "rtdb.fh"
#include "cdft.fh"
#include "util.fh"
c
      integer rtdb
      logical oprint
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
c
c     Read/set/put the numerical quadrature input parameters.
c
c     reset gridtype; set nagrid, nradpts, nrad, nang
c
      if (.not. rtdb_cput(rtdb, 'dft:gridtype', 1, gridtype))
     &   call errquit('grid_update_rtdb: rtdb_put failed', 1, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb, 'dft:nagrid', mt_int, 1, nagrid))
     &   call errquit('grid_update_rtdb: rtdb_put failed', 2, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'dft:nradpts', mt_int, 1, nradpts))
     &   call errquit('grid_update_rtdb: rtdb_put failed', 3, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'dft:nrad', mt_int, dft_ntags_bsmx, nrad))
     &   call errquit('grid_update_rtdb: rtdb_put failed', 4, RTDB_ERR)
c
      if (.not. rtdb_put(rtdb,'dft:nang', mt_int, dft_ntags_bsmx, nang))
     &   call errquit('grid_update_rtdb: rtdb_put failed', 5, RTDB_ERR)
      return
      end
