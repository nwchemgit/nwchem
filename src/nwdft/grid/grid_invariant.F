      Subroutine grid_invariant
c
c$Id: grid_invariant.F,v 1.5 1999-11-02 02:40:43 edo Exp $
c
      implicit none
c      
#include "mafdecls.fh"
#include "util.fh"
#include "global.fh"
#include "stdio.fh"
#include "cdft.fh"
#include "geom.fh"
c
      logical oprint, oprint_grid
      integer lcoord, icoord, lcharge, icharge, ltags, itags
c
c     Set print options.
c
      oprint = util_print('quadrature', print_high)
      oprint_grid = util_print('griddebug', print_debug)
c
c     Define the rotation matrix that will be used to ensure that the quadrature
c     grid is rotationally invariant.
c
      if (iqrot.eq.0)then
c
c        Turn off the rotational invariance by defining a unit matrix.
c
         qrot(1,1) = 1.D0
         qrot(2,1) = 0.D0
         qrot(3,1) = 0.D0
         qrot(1,2) = 0.D0
         qrot(2,2) = 1.D0
         qrot(3,2) = 0.D0
         qrot(1,3) = 0.D0
         qrot(2,3) = 0.D0
         qrot(3,3) = 1.D0
c
      else
c
c        allocate space for atomic coordinates and charges
c
         if (.not. Ma_Push_Get(MT_Dbl,ncenters*3,'coordinates',lcoord,
     &      icoord))call errquit('grid_invariant: failed to alloc 
     &      coordinates',0)
         if (.not. Ma_Push_Get(MT_Dbl,ncenters,'charges',lcharge,
     &      icharge))call errquit('grid_invariant: failed to alloc 
     &      charges',0)
         if (.not. Ma_Push_Get(MT_Byte, ncenters*16, 'center tags',
     &      ltags, itags))call errquit('grid_invariant: failed to alloc 
     &      center tags',0)
c     
         if (.not. geom_cart_get(geom, ncenters, Byte_MB(itags),
     &                           Dbl_MB(icoord), Dbl_MB(icharge)))
     &        call errquit('grid_invariant: geom_cart_get failed',74)
c
         call grid_invar_mat(dbl_mb(icoord), dbl_mb(icharge))
         if (.not. MA_Pop_Stack(ltags))
     &      call errquit('grid_invariant: pop stack failed.',0)
         if (.not. MA_Pop_Stack(lcharge))
     &      call errquit('grid_invariant: pop stack failed.',0)
         if (.not. MA_Pop_Stack(lcoord))
     &      call errquit('grid_invariant: pop stack failed.',0)
      endif
c
c     debug writes
c
      if (ga_nodeid().eq.0.and.oprint_grid)then
         write(6,*)' qrot: '
         call output(qrot, 1, 3, 1, 3, 3, 3, 1)
      endif
      return
      end
      Subroutine grid_invar_mat(coord, charge)
c
c$Id: grid_invariant.F,v 1.5 1999-11-02 02:40:43 edo Exp $
c
      implicit none
c
#include "stdio.fh"
#include "cdft.fh"
c      
      double precision eps
      parameter (eps = 1.d-20)
      double precision coord(3,ncenters), charge(ncenters)
c
c     allocate space for M and T (rotational invariance; 
c     see B.G. Johnson et al, CPL, 220, 377 (1994), and auxiliary 
c     species R_a - T and (R_a -T)(R_a - T)^T
c
      double precision tvec(3), ramt(3), ramtm(3,3), eigs(3), work(12)
      double precision ztot, lramt
c
      integer ictr, icart, jcart, info
c
      call dzero(qrot, 9)
      call dzero(tvec, 3)
c
c     atomic case 
c
      if(ncenters.eq.1) then
         qrot(1,1)=1d0
         qrot(2,2)=1d0
         qrot(3,3)=1d0
         return
      endif
c
c     calculate the T vector
c
      ztot = 0.d0
      do ictr = 1, ncenters
         ztot = ztot+charge(ictr)
         do icart = 1, 3
            tvec(icart) = tvec(icart) + 
     &                    charge(ictr)*coord(icart,ictr)
         enddo
      enddo
c
c     normalize T
c
      if (abs(ztot).gt.EPS) then
         do icart = 1, 3
            tvec(icart) = tvec(icart)/ztot
         enddo
      else
         write(LuOut,*)
         write(LuOut,*) 'grid_invariant:  ',
     &                  'The total nuclear charge is zero.'
         call errquit('Exiting in grid_invariant.',4)
      endif
c
c     calculate the M matrx
c
      do ictr = 1, ncenters
c
c        Form R_a - T
c
         call dzero(ramt, 3)
         lramt = 0.d0
         do icart = 1, 3
            ramt(icart) = coord(icart,ictr) - tvec(icart)
            lramt = lramt + ramt(icart)*ramt(icart)
         enddo
c
c        form (R_a -T)(R_a - T)^T
c
         do icart = 1, 3
            do jcart = 1, 3
               ramtm(icart,jcart) = ramt(icart)*ramt(jcart)
            enddo
         enddo
c
c        form M in qrot
c
         do icart = 1, 3
            qrot(icart,icart) = qrot(icart,icart) + 
     &                          charge(ictr)*lramt
            do jcart = 1, 3
               qrot(icart,jcart) = qrot(icart,jcart) - 
     &                             ramtm(icart,jcart)*charge(ictr)
            enddo
         enddo
      enddo
c
c     diagonalize M
c
      call dzero(eigs, 3)
      call dzero(work, 12)
      call dsyev('V', 'L', 3, qrot, 3, eigs, work, 12, info)
      if (info.ne.0) then
         write(LuOut,*)
         write(LuOut,*) 'grid_invariant:  failure in dsyev.'
         call errquit('Exiting in grid_invariant.',5)
      endif
      return
      end
