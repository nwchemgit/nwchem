      Subroutine grid_gen_pts(rtdb)
c
c$Id: grid_gen_pts.F,v 1.11 2000-03-17 02:58:56 edo Exp $
c
      implicit none
c      
#include "inp.fh"
#include "rtdb.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "tcgmsg.fh"
#include "util.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "geom.fh"
#include "sym.fh"
c
      integer rtdb
      logical  grid_on_file
      integer lcoord, icoord, lcharge, icharge, ltags, itags
      integer mxnrad, nqlen, lqs, iqs
      integer laqbuf, iaqbuf
      integer ang_quad_buf_size, nq_tot
c
c     Need coordinates of all atoms.     
c
c     allocate space for atomic coordinates and charges
c
      if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &   icoord))call errquit(
     (     'grid_gen_pts: failed to alloc coordinates',0)
      if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &   icharge))call errquit(
     (     'grid_gen_pts: failed to alloc charges',0)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit(
     (     'grid_gen_pts: failed to alloc center tags',0)
c
c     Get ncenter tags, coordinates, and charges from the geometry object.
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &     call errquit('grid_gen_pts: geom_cart_get failed',74)

c
c     A list of pruned (or unpruned) r, theta and phi is generated here.
c
      if (.not. rtdb_get(rtdb, 'dft:mxnrad', mt_int, 1, mxnrad))
     &   call errquit('grid_gen: rtdb_get failed', 119)
c
      nqlen = ncenters*mxnrad
      if (.not.MA_Push_get(MT_int,4*nqlen,'shell list',
     &                     lqs,iqs))
     &   call errquit('grid_gen: cannot allocate shell list',0)
      call grid_list(rtdb, int_mb(iqs), nqlen, nq_tot)
c
c     The points (rotationally invariant if desired) are generated and written to disk
c     with fixed number of points (buffer size).
c
c     buffer size needed to hold largest angular quadrature is nqmax
c
!      ang_quad_buf_size = nqmax
c
c     hardwire!!
c
      ang_quad_buf_size = 5810
c
      if (.not.MA_Push_get(MT_dbl, 4*ang_quad_buf_size,
     &                     'ang quad pts buffer', laqbuf, iaqbuf))
     &   call errquit('grid_gen: cannot allocate quad pt buffer', 0)
c
      if (.not. rtdb_get(rtdb, 'dft:grid_on_file', mt_log, 1,
     &   grid_on_file))then
         grid_on_file = .false.
      endif
      if (grid_on_file)then
         call grid_tasks(rtdb, Dbl_MB(icoord), int_mb(iqs), nqlen,
     &                   dbl_mb(iaqbuf), ang_quad_buf_size,
     ,        nq_tot)
c
c        try to balance the number of tasks on each node
c
cedo         call grid_file_balance()
c
      endif
      if (.not. MA_Pop_Stack(laqbuf))
     &   call errquit('grid_gen_pts: pop stack failed.',0)
      if (.not. MA_Pop_Stack(lqs))
     &   call errquit('grid_gen_pts: pop stack failed.',0)
      if (.not. MA_Pop_Stack(ltags))
     &   call errquit('grid_gen_pts: pop stack failed.',0)
      if (.not. MA_Pop_Stack(lcharge))
     &   call errquit('grid_gen_pts: pop stack failed.',0)
      if (.not. MA_Pop_Stack(lcoord))
     &   call errquit('grid_gen_pts: pop stack failed.',0)
      return
      end
