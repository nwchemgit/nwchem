      Subroutine grid_gen_pts(rtdb)
c
c$Id$
c
      implicit none
#include "errquit.fh"
c      
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "cdft.fh"
#include "geom.fh"
#include "cgridfile.fh"
c
      integer rtdb
      integer lcoord, icoord, lcharge, icharge, ltags, itags
      integer mxnrad, nqlen, lqs, iqs
      integer ang_quad_buf_size, nq_tot
      logical grid_file_balance,grid_getnsubb,oldgrid
      external  grid_file_balance,grid_getnsubb
c
c     Need coordinates of all atoms.     
c
c     allocate space for atomic coordinates and charges
c
      if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &   icoord))call errquit(
     (     'grid_gen_pts: failed to alloc coordinates',0, MA_ERR)
      if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &   icharge))call errquit(
     (     'grid_gen_pts: failed to alloc charges',0, MA_ERR)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit(
     (     'grid_gen_pts: failed to alloc center tags',0, MA_ERR)
c
c     Get ncenter tags, coordinates, and charges from the geometry object.
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &     call errquit('grid_gen_pts: geom_cart_get failed',74,
     &       GEOM_ERR)
      oldgrid=.true.
c
c     A list of pruned (or unpruned) r, theta and phi is generated here.
c
      if (.not. rtdb_get(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &   call errquit('grid_gen: rtdb_get failed', 119, RTDB_ERR)
c
      if(oldgrid) then
      nqlen = ncenters*mxnrad
      if (.not.MA_Push_get(MT_int,4*nqlen,'shell list',
     &                     lqs,iqs))
     &   call errquit('grid_gen: cannot allocate shell list',0, MA_ERR)
      call grid_list(rtdb, int_mb(iqs), nqlen, nq_tot)
      endif
c
c     The points (rotationally invariant if desired) are generated and written to disk
c     with fixed number of points (buffer size).
c
c     buffer size needed to hold largest angular quadrature is nqmax
c
      ang_quad_buf_size = nqmax
c
      if (.not. MA_Pop_Stack(lqs))
     &   call errquit('grid_gen_pts: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(ltags))
     &   call errquit('grid_gen_pts: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(lcharge))
     &   call errquit('grid_gen_pts: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(lcoord))
     &   call errquit('grid_gen_pts: pop stack failed.',0, MA_ERR)
      return
      end
