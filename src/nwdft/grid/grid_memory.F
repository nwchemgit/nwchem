      Subroutine grid_memory(rtdb)
c
c$Id$
c
      implicit none
#include "errquit.fh"
c      
#include "rtdb.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "geom.fh"
c
      integer rtdb
      logical oprint
      integer lcoord, icoord, lcharge, icharge, ltags, itags
      integer lqxyz, iqxyz, lrq0, irq0, lexpo, iexpo
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
c
c     allocate space for atomic coordinates and charges
c
      if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &   icoord))call errquit(
     .     'grid_memory: failed to alloc coordinates',0, MA_ERR)
      if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &   icharge))call errquit(
     .     'grid_memory: failed to alloc charges',0, MA_ERR)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit(
     .     'grid_memory: failed to alloc center tags',0, MA_ERR)
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &   call errquit('rdinput: geom_cart_get failed',74, GEOM_ERR)
c
      if (.not.MA_Push_get(mt_dbl,3*nqmax,'qxyz',lqxyz,iqxyz))
     &   call errquit('grid_memory: cannot allocate qxyz',0, MA_ERR)
      if (.not.MA_Push_get(mt_dbl,ncenters,'rq0',lrq0,irq0))
     &   call errquit('grid_memory: cannot allocate rq0',0, MA_ERR)
      if (.not.MA_Push_Get(MT_Dbl,nbf_ao_mxprim,'expo',lexpo,iexpo))
     &   call errquit('grid_memory: cannot allocate expo',0, MA_ERR)
c
      call mbf_ao_max(rtdb,  dbl_mb(irq0), 
     &                dbl_mb(iexpo), dbl_mb(icoord))
c
      if (.not.ma_pop_stack(lexpo))
     &   call errquit('grid_memory: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lrq0))
     &   call errquit('grid_memory: cannot pop stack',0, MA_ERR)
      if (.not.ma_pop_stack(lqxyz))
     &   call errquit('grid_memory: cannot pop stack',0, MA_ERR)
      if (.not. MA_Pop_Stack(ltags))
     &   call errquit('grid_memory: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(lcharge))
     &   call errquit('grid_memory: pop stack failed.',0, MA_ERR)
      if (.not. MA_Pop_Stack(lcoord))
     &   call errquit('grid_memory: pop stack failed.',0, MA_ERR)
      return
      end
