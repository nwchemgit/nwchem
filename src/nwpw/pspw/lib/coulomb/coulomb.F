
*
* $Id: coulomb.F,v 1.9 2009-02-07 03:50:47 bylaska Exp $
*

      subroutine coulomb_init()
      implicit none 
#include "errquit.fh"

#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
c     real*8 vc(nfft3d)
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl

*     **** local variables ****
      integer npack0,nfft3d,G(3)
      integer i,j,k
      integer zero,qzero,pzero,taskid
      integer nx,ny,nxh,nyh
      real*8  fourpi,gg
      logical value
      integer tmp1(2)

*     **** external functions ****
*     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx

      call nwpw_timing_start(7)
      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

*     **** allocate vc memory ****
      value = MA_alloc_get(mt_dbl,npack0,'vc',vc_hndl,vc_indx)
      if (.not. value) call errquit('out of heap memory',0, MA_ERR)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      call Parallel2d_taskid_i(taskid)
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      nxh=nx/2
      nyh=ny/2

*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
c     call D3dB_ktoqp(1,k+1,qzero,pzero)
c     zero = (qzero-1)*(nxh+1)*ny
c    >     + j*(nxh+1)
c    >     + i+1
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)

*     ***** form Vc = 4*pi/G**2  *****
      fourpi = 4.0d0*(4.0d0*datan(1.0d0))
      do i = 1,nfft3d
         
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

         if ((pzero.eq.taskid) .and. (i.eq.zero)) then
            dbl_mb(tmp1(1)+i-1) = 0.0d0
         else
            dbl_mb(tmp1(1)+i-1) = fourpi/gg
         end if
        
      end do
      call Pack_t_pack(0,dbl_mb(tmp1(1)))
      call Pack_t_Copy(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx))
      value = MA_pop_stack(tmp1(2))

      call nwpw_timing_end(7)


      return
      end

      subroutine coulomb_end()
      implicit none
#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl
      logical value

      value = MA_free_heap(vc_hndl)
      return
      end


      subroutine coulomb_v(dng,vc_out)
      implicit none
      complex*16 dng(*)
      complex*16 vc_out(*)

#include "mafdecls.fh"


*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl

      call nwpw_timing_start(7)
      call Pack_tc_Mul(0,dbl_mb(vc_indx),dng,vc_out)
      call nwpw_timing_end(7)

      return 
      end


      real*8 function coulomb_e(dng)
      implicit none
#include "errquit.fh"
      complex*16 dng(*)

#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
*     real*8 vc(nfft3d)
*     common / vc_block / vc
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl


*     **** local variables ****
      integer npack0
      real*8 ec

c     real*8  tmp1(*)
      integer tmp1(2)
      logical value

*     **** external functions ****
      real*8 lattice_omega
      external lattice_omega

      call nwpw_timing_start(7)
      call Pack_npack(0,npack0)
      value = MA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      
      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),ec)
       
      ec = 0.5d0*ec*lattice_omega()

      value = MA_pop_stack(tmp1(2))
      call nwpw_timing_end(7)

      coulomb_e = ec
      return 
      end





      subroutine coulomb_euv(dng,euv)
*
* $Id: coulomb.F,v 1.9 2009-02-07 03:50:47 bylaska Exp $
*
      implicit none 
#include "errquit.fh"
      complex*16 dng(*)
      real*8 euv(3,3)

#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      integer vc_indx,vc_hndl
      common / vc_block / vc_indx,vc_hndl


*     **** local variables ****
      integer npack0,nfft3d,G(2,3)
      integer i,j
      integer u,v,s
      logical value

      real*8 pi,fourpi,scal,ss,sum
      real*8 hm(3,3),Bus(3,3),ecoul 
      integer tmp1(2),tmp2(2)

*     **** external functions ****
c     real*8 G(nfft3d,3)
      integer  G_indx
      external G_indx

      real*8   lattice_unitg,lattice_omega,coulomb_e
      external lattice_unitg,lattice_omega,coulomb_e



      pi     = 4.0d0*datan(1.0d0)
      fourpi = 4.0d0*pi
      scal   = 1.0d0/(2.0d0*pi)

*     *** define hm ****
      do j=1,3
      do i=1,3
         hm(i,j) = scal*lattice_unitg(i,j)
      end do
      end do      


      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)

      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G1',G(2,1),G(1,1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G2',G(2,2),G(1,2))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G3',G(2,3),G(1,3))
      if (.not. value) call errquit('out of stack  memory',0, MA_ERR)

      value = MA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      value = MA_push_get(mt_dbl,npack0,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      call dcopy(nfft3d,dbl_mb(G_indx(1)),1,dbl_mb(G(1,1)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(2)),1,dbl_mb(G(1,2)),1)
      call dcopy(nfft3d,dbl_mb(G_indx(3)),1,dbl_mb(G(1,3)),1)
      call Pack_t_pack(0,dbl_mb(G(1,1)))
      call Pack_t_pack(0,dbl_mb(G(1,2)))
      call Pack_t_pack(0,dbl_mb(G(1,3)))

*     **** tmp2(G) = (n(G)**2) * (4*pi/G**2)**2  ****
      call Pack_ct_Sqr(0,dng,dbl_mb(tmp1(1)))
      call Pack_tt_Mul(0,dbl_mb(tmp1(1)),dbl_mb(vc_indx),
     >                                   dbl_mb(tmp2(1)))
c      call Pack_tt_Mul(0,dbl_mb(tmp2(1)),dbl_mb(vc_indx),
c     >                                   dbl_mb(tmp2(1)))
      call Pack_tt_Mul2(0,dbl_mb(vc_indx),dbl_mb(tmp2(1)))


*     **** Bus = Sum(G) (omega/4*pi)*tmp2(G)*Gu*Gs ****
      call dcopy(9,0.0d0,0,Bus,1)
      ss     = lattice_omega()/fourpi
      do u=1,3
      do s=u,3
        call Pack_tt_Mul(0,dbl_mb(G(1,u)),
     >                     dbl_mb(G(1,s)),
     >                     dbl_mb(tmp1(1)))
        call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(tmp2(1)),sum)

        Bus(u,s) = ss*sum
      end do
      end do
      do u=1,3
      do s=u+1,3
         Bus(s,u) = Bus(u,s)
      end do
      end do

      ecoul = coulomb_e(dng)
      do v=1,3
      do u=1,3
         euv(u,v) = -ecoul*hm(u,v)
         do s=1,3
            euv(u,v) = euv(u,v) + Bus(u,s)*hm(s,v)
         end do
      end do
      end do
 
      value = MA_pop_stack(tmp2(2))
      value = MA_pop_stack(tmp1(2))
      value = MA_pop_stack(G(2,3)) 
      value = MA_pop_stack(G(2,2)) 
      value = MA_pop_stack(G(2,1))

      return
      end
