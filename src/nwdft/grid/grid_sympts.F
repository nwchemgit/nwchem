      subroutine grid_sympts0(geom,nops,
     ,     k_symops,l_symops,
     ,     k_gridscr,l_gridscr,k_iscr,l_iscr,mxang)
c
c$Id: grid_sympts.F,v 1.13 2001-10-25 02:57:36 edo Exp $
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
     '        '     grid_sympt: mpush get failed',1)
         call sym_ops_get(geom,nops,dbl_mb(k_symops))
         call grid_c5trash(nops,dbl_mb(k_symops))
      else
         l_symops=-999
      endif
      if (.not.MA_push_Get(MT_Dbl,mxang*4,'gridscr',
     ,     l_gridscr,k_gridscr)) call errquit(
     '     '     grid_sympt: mpush get failed',2)
      if (.not.MA_push_Get(MT_int,mxang,'iscr',
     ,     l_iscr,k_iscr)) call errquit(
     '     '     grid_sympt: mpush get failed',3)
      return
      end
      subroutine grid_sympts2(l_symops,l_gridscr)
      implicit none
#include "mafdecls.fh"
      integer l_symops ! [in]  symops handle
      integer l_gridscr ! [in]  grid handle
c      
c     
      if (.not.MA_chop_stack(l_gridscr)) call errquit(
     '     '     grid_sympt: free stack failed',1)
      if(l_symops.ne.-999) then
      if (.not.MA_pop_stack(l_symops)) call errquit(
     '     '     grid_sympt: free stack failed',2)
      endif
      return
      end
      Subroutine grid_sympts(dosympts,r,xyz, ictr, irsh, 
     ,     iang, iang_done,nq_left,
     ,     aqbuf, nq, fnxyz,nxyz,
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
      double precision aqbuf(4,*) ! quadr pts and wght
      integer fnxyz      ! [in] cube def
      integer nx,ny,nz      !  cube def
      double precision r ! [in]
      integer nxyz
      double precision nxyz2
      logical grid_incube
      external grid_incube
      integer iang_done ! [in/out] angular quadra analyzed
c
c     Compute the quadrature points & weights for a given 
c     set of radial shells.
c     Keep only point belonging to cube
c
      integer k, ia_ictr,nq_left,k1,h,n_equiv
      double precision  tx, ty, tz, xor, yor, zor
      double precision x,y,z,x1,y1,z1,ddd,wsum
c
      oprint_grid = util_print('griddebug', print_debug)
c
c
      call grid_n2nxyz(nxyz,fnxyz,nx,ny,nz)
      nxyz2=nxyz*0.5d0
      ia_ictr = iatype(ictr)
      if(iang_done.ne.iang) then
         nq_left = 0
         do k = 1, ntheta(iang)
c     
            tx=za_leb(1,k,iang)
            ty=za_leb(2,k,iang)
            tz=za_leb(3,k,iang)
            nq_left = nq_left + 1
            
            gridscr(1,nq_left) =
     .           qrot(1,1)*tx + qrot(1,2)*ty + qrot(1,3)*tz
            gridscr(2,nq_left) =  
     .           qrot(2,1)*tx + qrot(2,2)*ty + qrot(2,3)*tz
            gridscr(3,nq_left) =  
     .           qrot(3,1)*tx + qrot(3,2)*ty + qrot(3,3)*tz
            gridscr(4,nq_left) = wa_leb(k,iang)
         enddo
c
c     apply symmetry operations
c
         if(dosympts) then
            do k=1,nq_left
               if(gridscr(4,k).ne.0d0) then
                  n_equiv=1
                  x=xyz(1) +gridscr(1,k)
                  y=xyz(2) +gridscr(2,k)
                  z=xyz(3) +gridscr(3,k)
                  do h=1,nops
                     x1=ops(1,h,1)*x+ops(1,h,2)*y+ops(1,h,3)*z-xyz(1)
                     y1=ops(2,h,1)*x+ops(2,h,2)*y+ops(2,h,3)*z-xyz(2)
                     z1=ops(3,h,1)*x+ops(3,h,2)*y+ops(3,h,3)*z-xyz(3)
                     do k1=k+1,nq_left
                        ddd=(x1-gridscr(1,k1))*(x1-gridscr(1,k1))+
     +                       (y1-gridscr(2,k1))*(y1-gridscr(2,k1))+
     +                       (z1-gridscr(3,k1))*(z1-gridscr(3,k1))
                        if(ddd.lt.1d-10) then
                           if(gridscr(4,k1).ne.0d0) then
                              gridscr(4,k1)=0d0
                              n_equiv=n_equiv+1
                           endif
                           goto 111
                        endif
                     enddo
cedo                  write(6,*) ' no equiv found ???',222
 111                 continue
                  enddo
                  gridscr(4,k)=n_equiv*gridscr(4,k)

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
     '           'grid_sympts: lost weights normaliz ',
     ,           log10((abs(wsum-1d0))))
         endif
         iang_done=iang
      endif

c     
c
      nq=0
      do k = 1, nq_left
         xor=gridscr(1,k)*r
         yor=gridscr(2,k)*r
         zor=gridscr(3,k)*r
         if(fnxyz.eq.-99.or.
     .        grid_incube(xor,yor,zor,nx,ny,nz,fnxyz,nxyz2,
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
      subroutine grid_c5trash(nops,symops)
      implicit none
      integer nops,nops_red
      double precision symops(3,nops,4)
      double precision deter3
      external deter3
c
      integer ii
      double precision dd,matt(3,3),trtr
      integer angle
      double precision pi
      double precision aa(3,120,4)
      integer iscr(120)
c     
c     remove any 5-fold operation
c     since it does not belong to Oh
c
      pi=4d0*atan(1d0)
      nops_red=0
      do ii=1,nops
         matt(1,1)=symops(1,ii,1)
         matt(2,1)=symops(2,ii,1)
         matt(3,1)=symops(3,ii,1)
         matt(1,2)=symops(1,ii,2)
         matt(2,2)=symops(2,ii,2)
         matt(3,2)=symops(3,ii,2)
         matt(1,3)=symops(1,ii,3)
         matt(2,3)=symops(2,ii,3)
         matt(3,3)=symops(3,ii,3)
         dd=deter3(matt)
         trtr=matt(1,1)+matt(2,2)+matt(3,3)
         angle=abs(nint(acos(mod(0.5d0*(trtr-dd),1d0))*180d0/pi))
         if(mod(angle,36).eq.0.and.(angle/36).lt.5) then
c     do nothing
         else
            nops_red=nops_red+1
            iscr(nops_red)=ii
         endif
      enddo
      call dcopy(3*4*nops,symops,1,aa,1)
      call grid_symcomp(nops_red,symops,nops,aa,iscr)
      nops=nops_red
      return
      end
      subroutine grid_symcomp(nops_red,symops_red,nops,symops,iscr)
      implicit none
      integer nops,nops_red
      double precision symops_red(3,nops_red,4)
      double precision symops(3,nops,4)
      integer kk,iscr(72)
c
c     compact array to get nops_red stride
c
      do kk=1,nops_red
         symops_red(1,kk,1)=symops(1,iscr(kk),1)
         symops_red(2,kk,1)=symops(2,iscr(kk),1)
         symops_red(3,kk,1)=symops(3,iscr(kk),1)
         symops_red(1,kk,2)=symops(1,iscr(kk),2)
         symops_red(2,kk,2)=symops(2,iscr(kk),2)
         symops_red(3,kk,2)=symops(3,iscr(kk),2)
         symops_red(1,kk,3)=symops(1,iscr(kk),3)
         symops_red(2,kk,3)=symops(2,iscr(kk),3)
         symops_red(3,kk,3)=symops(3,iscr(kk),3)
      enddo
      return
      end
