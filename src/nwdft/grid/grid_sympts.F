      subroutine grid_sympts0(geom,nops,k_symops,l_symops,
     ,     k_gridscr,l_gridscr,k_iscr,l_iscr,mxang)
c
c$Id: grid_sympts.F,v 1.3 2000-03-30 03:06:00 edo Exp $
c
      implicit none
#include "mafdecls.fh"
      integer geom
      integer nops ! [out]
      integer mxang ! [in] max no. angular pts
      integer k_symops,l_symops ! [out]  symops handles
      integer k_gridscr,l_gridscr ! [out]  grid handles
      integer k_iscr,l_iscr ! [out] iscr handles
c      
      integer sym_number_ops
      external sym_number_ops
c     
      nops = sym_number_ops(geom)
      if(nops.ne.0) then
         if (.not.MA_push_Get(MT_Dbl,nops*3*4,'symops',
     ,        l_symops,k_symops)) call errquit(
     '        '     grid_sympt: malloc get failed',1)
         call sym_ops_get(geom,nops,dbl_mb(k_symops))
      else
         l_symops=-999
      endif
      if (.not.MA_push_Get(MT_Dbl,mxang*4,'gridscr',
     ,     l_gridscr,k_gridscr)) call errquit(
     '     '     grid_sympt: malloc get failed',2)
      if (.not.MA_push_Get(MT_int,mxang,'iscr',
     ,     l_iscr,k_iscr)) call errquit(
     '     '     grid_sympt: malloc get failed',3)
      return
      end
      subroutine grid_sympts2(l_symops,l_gridscr,l_iscr)
      implicit none
#include "mafdecls.fh"
      integer l_symops ! [in]  symops handle
      integer l_gridscr,l_iscr ! [in]  grid handle
c      
c     
      if (.not.MA_pop_stack(l_iscr)) call errquit(
     '     '     grid_sympt: free stack failed',0)
      if (.not.MA_pop_stack(l_gridscr)) call errquit(
     '     '     grid_sympt: free stack failed',1)
      if(l_symops.ne.-999) then
      if (.not.MA_pop_stack(l_symops)) call errquit(
     '     '     grid_sympt: free stack failed',2)
      endif
      return
      end
      Subroutine grid_sympts(dosympts,r,xyz, ictr, irsh, 
     ,     iang, iang_done,nq_left,
     ,     aqbuf, nq, nx,ny,nz,nxyz,
     ,     gridscr,nops,ops,iscr)
      implicit none
#include "cdft.fh"
#include "util.fh"
#include "global.fh"
c      
      logical dosympts ! [in]
      integer iang ! [input]
      integer ictr ! [input]
      integer irsh ! [input]
      integer nq   ! [input]
      integer iscr(*)
c      
      double precision gridscr(4,*)
      integer nops
      double precision ops(3,nops,4)
      double precision xyz(3) ! [input] quadr center
      logical  oprint_grid
      double precision aqbuf(4,*), qwght_nq ! quadr pts and wght
      integer nx,ny,nz      ! [in] cube def
      double precision r ! [in]
      integer nxyz
      logical grid_incube
      external grid_incube
      integer iang_done ! [in/out] angular quadra analyzed
c
c     Compute the quadrature points & weights for a given 
c     set of radial shells.
c     Keep only point belonging to cube
c
      integer k, l, ia_ictr,nq_left,k1,h,n_equiv
      double precision  tx, ty, tz, xor, yor, zor
      double precision x,y,z,x1,y1,z1,ddd,wsum
c
      oprint_grid = util_print('griddebug', print_debug)
c
c
      ia_ictr = iatype(ictr)
      if(iang_done.ne.iang) then
         nq_left = 0
         do k = 1, ntheta(iang)
c     
            tx=za_leb(1,k,iang)
            ty=za_leb(2,k,iang)
            tz=za_leb(3,k,iang)
            xor=(qrot(1,1)*tx + qrot(1,2)*ty + qrot(1,3)*tz)
            yor=(qrot(2,1)*tx + qrot(2,2)*ty + qrot(2,3)*tz)
            zor=(qrot(3,1)*tx + qrot(3,2)*ty + qrot(3,3)*tz)
            nq_left = nq_left + 1
            
            gridscr(1,nq_left) =  xor
            gridscr(2,nq_left) =  yor
            gridscr(3,nq_left) =  zor
            gridscr(4,nq_left) = wa_leb(k,iang)
!            gridscr(4,nq_left) = wa_leb(k,iang)*rwghts(irsh,ia_ictr)
         enddo
cedo         write(6,'(4(A,I4))') ' iang ',iang,' initial nq ',nq_left,
cedo     ,        ' nops ',nops
cedo         do h=1,nops
cedo            write(6,*) ' h ',h
cedo            write(6,'(3f12.6)')
cedo     .           ops(1,h,1),ops(1,h,2),ops(1,h,3),
cedo     .           ops(2,h,1),ops(2,h,2),ops(2,h,3),
cedo     .           ops(3,h,1),ops(3,h,2),ops(3,h,3)
cedo         enddo
c
c     apply symmetry operations
c
         if(dosympts) then
            do k=1,nq_left
               if(gridscr(4,k).ne.0d0) then
cedo            write(6,'(" k ",i3,4(F12.6))') k,
cedo     .           gridscr(1,k),
cedo     .           gridscr(2,k),
cedo     .           gridscr(3,k),
cedo     .           gridscr(4,k)
                  n_equiv=1
                  x=xyz(1) +gridscr(1,k)
                  y=xyz(2) +gridscr(2,k)
                  z=xyz(3) +gridscr(3,k)
cedo            write(6,'(" k1 ",i3,4(F12.6))') k,x,y,z
                  do h=1,nops
                     x1=ops(1,h,1)*x+ops(1,h,2)*y+ops(1,h,3)*z
                     y1=ops(2,h,1)*x+ops(2,h,2)*y+ops(2,h,3)*z
                     z1=ops(3,h,1)*x+ops(3,h,2)*y+ops(3,h,3)*z
cedo            write(6,'("  h0",i3,4(F12.6))') h,x1,y1,z1
                     x1=x1-xyz(1) 
                     y1=y1-xyz(2) 
                     z1=z1-xyz(3) 
cedo            write(6,'("  h1",i3,4(F12.6))') h,x1,y1,z1
                     do k1=k+1,nq_left
                        if(k.ne.k1) then
                           ddd=(x1-gridscr(1,k1))**2+
     +                          (y1-gridscr(2,k1))**2+
     +                          (z1-gridscr(3,k1))**2
                           if(ddd.lt.1d-10) then
                              if(gridscr(4,k1).ne.0d0) then
                                 if(abs(gridscr(4,k1)-gridscr(4,k)).
     .                                gt.1d-5) then
cedo                              write(6,'(" w ",2(I4,1pe20.12))') 
cedo     .                                k,gridscr(4,k),k1,gridscr(4,k1)
                                 endif
                                 gridscr(4,k1)=0d0
                                 n_equiv=n_equiv+1
                              endif
                              goto 111
                           endif
                        endif
                     enddo
cedo                  write(6,*) ' no equiv found ???',222
 111                 continue
                  enddo
                  gridscr(4,k)=n_equiv*gridscr(4,k)

cedo            write(6,'(2(A,I4),f20.12)') 
cedo     ,           ' k ',k,' equiv = ',n_equiv,gridscr(4,k)
               endif
            enddo
         endif
c     remove zeroes
         call grid_compress(nq_left,gridscr,iscr,1d-19)
c
c     check weight normaliz
c
         if(dosympts) then
            wsum=0d0
            do k=1,nq_left
               wsum=wsum+gridscr(4,k)
            enddo
cedo         write(6,'(A,I4),A,1PE20.12)') '         final nq ',nq_left
            wsum=wsum/dacos(-1d0)*.25d0
            if(abs(wsum-1d0).gt.1d-6) call errquit(
     '           'grid_sympts: lost weights normaliz ',111)
         endif
         iang_done=iang
      endif

c     
c
      nq=0
      do k = 1, nq_left
cedo         xor=(gridscr(1,k)-xyz(1))*r+xyz(1)
cedo         yor=(gridscr(2,k)-xyz(2))*r+xyz(2)
cedo         zor=(gridscr(3,k)-xyz(3))*r+xyz(3)
         xor=gridscr(1,k)*r
         yor=gridscr(2,k)*r
         zor=gridscr(3,k)*r
         if(grid_incube(xor,yor,zor,nx,ny,nz,nxyz,
     ,        lcube(ia_ictr))) then
            nq=nq+1
            aqbuf(1,nq) = xyz(1) +xor
            aqbuf(2,nq) = xyz(2) +yor
            aqbuf(3,nq) = xyz(3) +zor
            aqbuf(4,nq) = gridscr(4,k)*rwghts(irsh,ia_ictr)

         endif
      enddo
      return
      end
