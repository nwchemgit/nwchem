*
* $Id: wgc.F 24638 2013-10-11 23:57:35Z bylaska $
*

*     ****************************************
*     *                                      *
*     *              wgc_init                *
*     *                                      *
*     ****************************************
      subroutine wgc_init(rho0)
      implicit none 
      real*8 rho0

#include "errquit.fh"
#include "mafdecls.fh"

*     **** common block used for coulomb.f ****
      real*8  wgc_rho0,wgc_tune1,wgc_tune2
      integer wgc_indx,wgc_hndl
      logical wgc_eq_tune
      common / wgc_block / wgc_rho0,wgc_tune1,wgc_tune2,
     >                     wgc_indx,wgc_hndl,wgc_eq_tune

*     **** local variables ****
      real*8 one3rd
      parameter (one3rd=1.0d0/3.0d0)
      double precision toll
      parameter (toll=1.0d-16)

      integer npack0,nfft3d,G(3)
      integer i,j,k
      integer zero,qzero,pzero,taskid
      integer nx,ny
      real*8  pi,gg,scal1,ss
      logical value
      integer tmp1(2)

*     **** external functions ****
      integer  G_indx
      external G_indx
      real*8   control_wgc_alphabeta
      external control_wgc_alphabeta
 

      call nwpw_timing_start(7)
      call Parallel2d_taskid_i(taskid)

      call D3dB_nfft3d(1,nfft3d)
      call Pack_npack(0,npack0)
      G(1) = G_indx(1)
      G(2) = G_indx(2)
      G(3) = G_indx(3)

      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      scal1 = 1.0d0/(nx*ny*nz)

      wgc_rho0  = rho0
      wgc_tune1 = control_wgc_alphabeta(1)
      wgc_tune2 = control_wgc_alphabeta(2)
      wgc_eq_tune = (dabs(wgc_tune1-wgc_tune2).lt.1.0d-6)

      ss = scal1/(2.0d0*wgc_tune1*wgc_tune2
     >                 *rho0**(wgc_tune1+wgc_tune2-2.0d0))

*     **** allocate vc memory ****
      value = MA_alloc_get(mt_dbl,npack0,'wgc',wgc_hndl,wgc_indx)
      if (.not.value) call errquit('wgc_init:out of heap',0,MA_ERR)

      value = MA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('wgc_init:out of stack',1,MA_ERR)


*     ***** find the G==0 point in the lattice *****
      i=0
      j=0
      k=0
      call D3dB_ijktoindexp(1,i+1,j+1,k+1,zero,pzero)


*     ***** form Vc = 4*pi/G**2  *****
      pi = (4.0d0*datan(1.0d0))
      kf = (3.0d0*pi*pi*wgc_rho0)**one3rd

      do i = 1,nfft3d
         
         gg  = ( dbl_mb(G(1)+i-1)*dbl_mb(G(1)+i-1)
     >         + dbl_mb(G(2)+i-1)*dbl_mb(G(2)+i-1)
     >         + dbl_mb(G(3)+i-1)*dbl_mb(G(3)+i-1) )

         if (((pzero.eq.taskid) .and. (i.eq.zero)).or.
     E    (abs(gg) .lt.toll))  then
            dbl_mb(tmp1(1)+i-1) = 0.0d0
         else
            eta = dsqrt(gg)/(2.0d0*kf)
            eta2 = eta*eta
            tf = 3.0d0*rho0/(kf*kf)
            vw = 4.0d0*rho0/gg
            aa1 = (1.0d0-eta2)/(4.0d0*eta)
            aa2 = dlog(dabs((1.0d0+eta)/(1.0d0-eta)))
            lind = tf*(0.5d0 + aa1*aa2)
            dbl_mb(tmp1(1)+i-1) = ss*(1.0d0/lind-1.0d0/vw-1.0d0/tf)
         end if
        
      end do
      call Pack_t_pack(0,dbl_mb(tmp1(1)))
      call Pack_t_Copy(0,dbl_mb(tmp1(1)),dbl_mb(wgc_indx))

      value = MA_pop_stack(tmp1(2))
      if (.not. value) call errquit('wgc_init:popping stack',3,MA_ERR)

      call nwpw_timing_end(7)


      return
      end

*     ****************************************
*     *                                      *
*     *              wgc_end                 *
*     *                                      *
*     ****************************************
      subroutine wgc_end()
      implicit none

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      real*8  wgc_rho0,wgc_tune1,wgc_tune2
      integer wgc_indx,wgc_hndl
      logical wgc_eq_tune
      common / wgc_block / wgc_rho0,wgc_tune1,wgc_tune2,
     >                     wgc_indx,wgc_hndl,wgc_eq_tune

      if (.not.MA_free_heap(wgc_hndl))
         call errquit('wgc_end:freeing heap',0,MA_ERR)

      return
      end

*     ****************************************
*     *                                      *
*     *              wgc_v                   *
*     *                                      *
*     ****************************************
      subroutine wgc_v(ispin,dn,v_out)
      implicit none
      integer ispin
      real*8  dn(*)
      real*8  v_out(*)

#include "mafdecls.fh"
#include "errquit.fh"


*     **** common block used for coulomb.f ****
      real*8  wgc_rho0,wgc_tune1,wgc_tune2
      integer wgc_indx,wgc_hndl
      logical wgc_eq_tune
      common / wgc_block / wgc_rho0,wgc_tune1,wgc_tune2,
     >                     wgc_indx,wgc_hndl,wgc_eq_tune

      integer nfft3d,n2ft3d
      

      call nwpw_timing_start(7)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)

      if (.not. MA_push_get(mt_dcpl,nfft3d,'rho',rho(2),rho(1))) then
         call errquit('wgc_v:out of stack',0,MA_ERR)
      end if
      if (.not. MA_push_get(mt_dcpl,nfft3d,'tmp',tmp(2),tmp(1))) then
         call errquit('wgc_v:out of stack',1,MA_ERR)
      end if

      !**** alpha==beta ****
      if (wgc_eq_tune) then
         call D3dB_rr_Sum(1,dn,dn(n2ft3d*(ispin-1),dcpl_mb(rho(1)))
         call D3dB_r_Power1(1,wgc_tune1,dcpl_mb(rho(1)))
         call Pack_c_SMul1(0,dcpl_mb(rho(1)))
         call D3dB_rc_fft3f(1,dcpl_mb(rho(1)))
         call Pack_c_pack(0,dcpl_mb(rho(1)))

         call Pack_tc_Mul(0,dbl_mb(wgc_indx),dcpl_mb(rho(1)),v_out)
         call Pack_c_SMul1(0,2.0d0*wgc_tune1,v_out)
         call Pack_c_unpack(0,v_out)
         call D3dB_cr_fft3f(1,v_out)

         call D3dB_rr_Sum(1,dn,dn(n2ft3d*(ispin-1),dcpl_mb(rho(1)))
         call D3dB_r_Power1(1,(wgc_tune1-1.0d0),dcpl_mb(rho(1)))
         call D3dB_rr_Multiply2(1,dcpl_mb(rho(1)),v_out)

      !**** alpha!=beta ****
      else
         call D3dB_rr_Sum(1,dn,dn(n2ft3d*(ispin-1),dcpl_mb(rho(1)))
         call D3dB_r_Power1(1,wgc_tune1,dcpl_mb(rho(1)))
         call Pack_c_SMul1(0,dcpl_mb(rho(1)))
         call D3dB_rc_fft3f(1,dcpl_mb(rho(1)))
         call Pack_c_pack(0,dcpl_mb(rho(1)))
         call Pack_tc_Mul(0,dbl_mb(wgc_indx),
     >                      dcpl_mb(rho(1)),
     >                      v_out)
         call Pack_c_SMul1(0,wgc_tune2,v_out)
         call Pack_c_unpack(0,v_out)
         call D3dB_cr_fft3f(1,v_out)
         call D3dB_rr_Sum(1,dn,dn(n2ft3d*(ispin-1),dcpl_mb(rho(1)))
         call D3dB_r_Power1(1,(wgc_tune2-1.0d0),dcpl_mb(rho(1)))
         call D3dB_rr_Multiply1(1,dcpl_mb(rho(1)),v_out)


         call D3dB_rr_Sum(1,dn,dn(n2ft3d*(ispin-1),dcpl_mb(rho(1)))
         call D3dB_r_Power1(1,wgc_tune2,dcpl_mb(rho(1)))
         call Pack_c_SMul1(0,dcpl_mb(rho(1)))
         call D3dB_rc_fft3f(1,dcpl_mb(rho(1)))
         call Pack_c_pack(0,dcpl_mb(rho(1)))
         call Pack_tc_Mul(0,dbl_mb(wgc_indx),
     >                      dcpl_mb(rho(1)),
     >                      dcpl_mb(tmp(1)))
         call Pack_c_SMul1(0,wgc_tune1,dcpl_mb(tmp(1)))
         call Pack_c_unpack(0,dcpl_mb(tmp(1)))
         call D3dB_cr_fft3f(1,dcpl_mb(tmp(1)))
         call D3dB_rr_Sum(1,dn,dn(n2ft3d*(ispin-1),dcpl_mb(rho(1)))
         call D3dB_r_Power1(1,(wgc_tune1-1.0d0),dcpl_mb(rho(1)))
         call D3dB_rrr_MultiplyAdd(1,dcpl_mb(rho(1)),dcpl_mb(tmp(1)),
     >                             v_out)
      end if

      call nwpw_timing_end(7)

      return 
      end

*     ****************************************
*     *                                      *
*     *              wgc_e                   *
*     *                                      *
*     ****************************************
      real*8 function wgc_e(dng)
      implicit none
      complex*16 dng(*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      real*8  wgc_rho0,wgc_tune1,wgc_tune2
      integer wgc_indx,wgc_hndl
      logical wgc_eq_tune
      common / wgc_block / wgc_rho0,wgc_tune1,wgc_tune2,
     >                     wgc_indx,wgc_hndl,wgc_eq_tune

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

      wgc_e = ec
      return 
      end



*     ****************************************
*     *                                      *
*     *              wgc_euv                 *
*     *                                      *
*     ****************************************
      subroutine wgc_euv(dng,euv)
      implicit none 
      complex*16 dng(*)
      real*8 euv(3,3)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** common block used for coulomb.f ****
      real*8  wgc_rho0,wgc_tune1,wgc_tune2
      integer wgc_indx,wgc_hndl
      logical wgc_eq_tune
      common / wgc_block / wgc_rho0,wgc_tune1,wgc_tune2,
     >                     wgc_indx,wgc_hndl,wgc_eq_tune

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
      if (.not. value) call errquit('wgc_euv:out of stack',0, MA_ERR)
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G2',G(2,2),G(1,2))
      if (.not. value) call errquit('wgc_euv:out of stack',0, MA_ERR)
      value = MA_push_get(mt_dbl,nfft3d,
     >                     'G3',G(2,3),G(1,3))
      if (.not. value) call errquit('wgc_euv:out of stack',0, MA_ERR)

      value = MA_push_get(mt_dbl,npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not. value) call errquit('wgc_euv:out of stack',0, MA_ERR)

      value = MA_push_get(mt_dbl,npack0,'tmp2',tmp2(2),tmp2(1))
      if (.not. value) call errquit('wgc_euv:out of stack',0, MA_ERR)

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
 
      value =           MA_pop_stack(tmp2(2))
      value = value.and.MA_pop_stack(tmp1(2))
      value = value.and.MA_pop_stack(G(2,3)) 
      value = value.and.MA_pop_stack(G(2,2)) 
      value = value.and.MA_pop_stack(G(2,1))
      if (.not. value) call errquit('wgc_euv:popping stack',0,MA_ERR)

      return
      end

