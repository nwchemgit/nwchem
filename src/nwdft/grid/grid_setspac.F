      subroutine grid_setspac(rtdb)
c
c$Id: grid_setspac.F,v 1.1 1999-04-23 16:30:36 d3h449 Exp $
c
      implicit none
      integer rtdb
c
#include "cdft.fh"
#include "rtdb.fh"
#include "geom.fh"
#include "mafdecls.fh"
#include "cgridfile.fh"
c
      integer iqt, ictr
      integer nqpts
      integer lcoord, icoord, lcharge, icharge, ltags, itags
      integer lrqbuf, irqbuf
      integer lrij, irij, lrq, irq, lp, ip
      integer nqpts_per_buf, ndbl
      logical grid_file_rewind
      external grid_file_rewind
c
c     Need coordinates of all atoms.     
c
c     allocate space for atomic coordinates and charges
c
      if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &   icoord))call errquit('grid_gen_pts: failed to alloc 
     &   coordinates',0)
      if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &   icharge))call errquit('grid_setspac: failed to alloc 
     &   charges',0)
      if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &   ltags, itags))call errquit('grid_setspac: failed to 
     &   alloc center tags',0)
c
c     Get ncenter tags, coordinates, and charges from the geometry object.
c     
      if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                        Dbl_MB(icoord), Dbl_MB(icharge)))
     &     call errquit('grid_setspac: geom_cart_get failed',74)
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
      ndbl = (ncenters*(ncenters+1))/2
      if (.not.MA_Push_Get(MT_Dbl,ndbl,'Rij',lRij,iRij))
     &   call errquit('  in grid_setspac: cannot allocate Rij',0)
      call a_dist(Dbl_MB(icoord), Dbl_MB(iRij), ncenters)
      if (.not.MA_Push_get(mt_dbl,ncenters*nqpts_per_buf,'rq',lrq,irq))
     &   call errquit('grid_setspac: cannot allocate rq',0)
      if (.not.MA_Push_get(mt_dbl,ncenters*nqpts_per_buf,'p',lp,ip))
     &   call errquit('grid_setspac: cannot allocate rq',0)
c
      do 100 iqt = 1, n_tot_tasks
c
         call grid_file_read(nqpts_per_buf, nqpts, ictr, dbl_mb(irqbuf))
         write(6,*)' ictr, nqpts: ', ictr, nqpts
         call grid_beckew(dbl_mb(irqbuf), dbl_mb(icoord), dbl_mb(irij), 
     &                    dbl_mb(irq), dbl_mb(ip), ictr, nqpts)
c
  100 continue
c
      if (.not. MA_Pop_Stack(lp))
     &   call errquit('grid_setspac: pop stack failed.',0)
      if (.not. MA_Pop_Stack(lrq))
     &   call errquit('grid_setspac: pop stack failed.',0)
      if (.not. MA_Pop_Stack(lrij))
     &   call errquit('grid_setspac: pop stack failed.',0)
      if (.not. MA_Pop_Stack(lrqbuf))
     &   call errquit('grid_setspac: pop stack failed.',0)
      if (.not. MA_Pop_Stack(ltags))
     &   call errquit('grid_setspac: pop stack failed.',0)
      if (.not. MA_Pop_Stack(lcharge))
     &   call errquit('grid_setspac: pop stack failed.',0)
      if (.not. MA_Pop_Stack(lcoord))
     &   call errquit('grid_setspac: pop stack failed.',0)
      return
      end 
