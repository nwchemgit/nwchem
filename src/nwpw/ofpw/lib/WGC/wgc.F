*
* $Id$
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
#include "bafdecls.fh"
#include "wgc.fh"


*     **** local variables ****
      real*8 one3rd
      parameter (one3rd=1.0d0/3.0d0)
      double precision toll
      parameter (toll=1.0d-16)

      integer npack0,nfft3d,G(3)
      integer i,j,k
      integer zero,qzero,pzero,taskid
      integer nx,ny,nz
      real*8  pi,gg,ss,aa1,aa2,tf,kf,lind,eta1,eta2,eta,vw
      logical value
      integer tmp1(2)

*     **** external functions ****
      integer  G_indx
      external G_indx
      real*8   control_wgc_alphabetalambda
      external control_wgc_alphabetalambda
 

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
      wgc_scal1 = 1.0d0/(nx*ny*nz)

      wgc_rho0    = rho0
      wgc_tune1   = control_wgc_alphabetalambda(1)
      wgc_tune2   = control_wgc_alphabetalambda(2)
      wgc_lmd     = control_wgc_alphabetalambda(3)
      wgc_eq_tune = (dabs(wgc_tune1-wgc_tune2).lt.1.0d-6)

      !*******************************************************************************
      !**** An extra scaling of 1/(N1*N2*N3) is done to the WGG_kernel(G) so that ****
      !**** the rho(G) does not need to be scaled by 1/(N1*N2*N3)                 ****
      !**** Might want to scal rho(G) properly to make the code more transparent  ****
      !*******************************************************************************
      ss = wgc_lmd*wgc_scal1/(2.0d0*wgc_tune1*wgc_tune2
     >                       *rho0**(wgc_tune1+wgc_tune2-2.0d0))

*     **** allocate vc memory ****
      value = BA_alloc_get(mt_dbl,npack0,'wgc',wgc_hndl,wgc_indx)
      if (.not.value) call errquit('wgc_init:out of heap',0,MA_ERR)

      value = BA_push_get(mt_dbl,nfft3d,'tmp1',tmp1(2),tmp1(1))
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
      call Pack_tt_dot(0,dbl_mb(wgc_indx),dbl_mb(wgc_indx),aa1)

      value = BA_pop_stack(tmp1(2))
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "wgc.fh"


      if (.not.BA_free_heap(wgc_hndl))
     >   call errquit('wgc_end:freeing heap',0,MA_ERR)

      return
      end


*     ****************************************
*     *                                      *
*     *              wgc_rho                 *
*     *                                      *
*     ****************************************
      real*8 function wgc_rho()
      implicit none

#include "wgc.fh"

      wgc_rho = wgc_rho0
      return
      end

*     ****************************************
*     *                                      *
*     *              wgc_alpha               *
*     *                                      *
*     ****************************************
      real*8 function wgc_alpha()
      implicit none

#include "wgc.fh"

      wgc_alpha = wgc_tune1
      return
      end

*     ****************************************
*     *                                      *
*     *              wgc_beta                *
*     *                                      *
*     ****************************************
      real*8 function wgc_beta()
      implicit none

#include "wgc.fh"

      wgc_beta = wgc_tune2
      return
      end


*     ****************************************
*     *                                      *
*     *              wgc_lambda              *
*     *                                      *
*     ****************************************
      real*8 function wgc_lambda()
      implicit none

#include "wgc.fh"

      wgc_lambda = wgc_lmd
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

#include "bafdecls.fh"
#include "errquit.fh"
#include "wgc.fh"

*     **** local variables ****
      logical value
      integer nfft3d,n2ft3d
      integer rho1(2),tmp1(2)

      call nwpw_timing_start(7)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)

      value =      BA_push_get(mt_dbl,2*n2ft3d,'rho1',rho1(2),rho1(1))
     >        .and.BA_push_get(mt_dbl,2*n2ft3d,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('wgc_v:out of stack',0,MA_ERR)

      !**** alpha==beta ****
      if (wgc_eq_tune) then

         !**** the rho^alpha(G) is not beinge scaled by 1/(N1*N2*N3) ****
         call D3dB_rr_Sum(1,dn,dn(1+n2ft3d*(ispin-1)),dbl_mb(rho1(1)))
         call D3dB_r_notZero_Ends(1,dbl_mb(rho1(1)))
         call D3dB_r_Power1(1,wgc_tune1,dbl_mb(rho1(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rho1(1)))
         !call D3dB_rc_fft3f(1,dbl_mb(rho1(1)))
         call D3dB_rc_pfft3f(1,0,dbl_mb(rho1(1)))
         call Pack_c_pack(0,dbl_mb(rho1(1)))

         call Pack_tc_Mul(0,dbl_mb(wgc_indx),dbl_mb(rho1(1)),v_out)
         call Pack_c_SMul1(0,2.0d0*wgc_tune1,v_out)
         !call Pack_c_SMul1(0,wgc_tune1,v_out)
         call Pack_c_unpack(0,v_out)
         !call D3dB_cr_fft3b(1,v_out)
         call D3dB_cr_pfft3b(1,0,v_out)

         call D3dB_rr_Sum(1,dn,dn(1+n2ft3d*(ispin-1)),dbl_mb(rho1(1)))
         call D3dB_r_notZero_Ends(1,dbl_mb(rho1(1)))
         call D3dB_r_Power1(1,(wgc_tune1-1.0d0),dbl_mb(rho1(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rho1(1)))
         call D3dB_rr_Multiply2(1,dbl_mb(rho1(1)),v_out)

      !**** alpha!=beta ****
      else
         !**** the rho^alpha(G) is not being scaled by 1/(N1*N2*N3) ****
         call D3dB_rr_Sum(1,dn,dn(1+n2ft3d*(ispin-1)),dbl_mb(rho1(1)))
         call D3dB_r_notZero_Ends(1,dbl_mb(rho1(1)))
         call D3dB_r_Power1(1,wgc_tune1,dbl_mb(rho1(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rho1(1)))
         !call D3dB_rc_fft3f(1,dbl_mb(rho1(1)))
         call D3dB_rc_pfft3f(1,0,dbl_mb(rho1(1)))
         call Pack_c_pack(0,dbl_mb(rho1(1)))
         call Pack_tc_Mul(0,dbl_mb(wgc_indx),
     >                      dbl_mb(rho1(1)),
     >                      v_out)
         call Pack_c_SMul1(0,wgc_tune2,v_out)
         call Pack_c_unpack(0,v_out)
         !call D3dB_cr_fft3b(1,v_out)
         call D3dB_cr_pfft3b(1,0,v_out)
         call D3dB_rr_Sum(1,dn,dn(1+n2ft3d*(ispin-1)),dbl_mb(rho1(1)))
         call D3dB_r_notZero_Ends(1,dbl_mb(rho1(1)))
         call D3dB_r_Power1(1,(wgc_tune2-1.0d0),dbl_mb(rho1(1)))
         call D3dB_rr_Multiply2(1,dbl_mb(rho1(1)),v_out)


         !**** the rho^beta(G) is not being scaled by 1/(N1*N2*N3) ****
         call D3dB_rr_Sum(1,dn,dn(1+n2ft3d*(ispin-1)),dbl_mb(rho1(1)))
         call D3dB_r_notZero_Ends(1,dbl_mb(rho1(1)))
         call D3dB_r_Power1(1,wgc_tune2,dbl_mb(rho1(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rho1(1)))
         !call D3dB_rc_fft3f(1,dbl_mb(rho1(1)))
         call D3dB_rc_pfft3f(1,0,dbl_mb(rho1(1)))
         call Pack_c_pack(0,dbl_mb(rho1(1)))
         call Pack_tc_Mul(0,dbl_mb(wgc_indx),
     >                      dbl_mb(rho1(1)),
     >                      dbl_mb(tmp1(1)))
         call Pack_c_SMul1(0,wgc_tune1,dbl_mb(tmp1(1)))
         call Pack_c_unpack(0,dbl_mb(tmp1(1)))
         !call D3dB_cr_fft3b(1,dbl_mb(tmp1(1)))
         call D3dB_cr_pfft3b(1,0,dbl_mb(tmp1(1)))
         call D3dB_rr_Sum(1,dn,dn(1+n2ft3d*(ispin-1)),dbl_mb(rho1(1)))
         call D3dB_r_notZero_Ends(1,dbl_mb(rho1(1)))
         call D3dB_r_Power1(1,(wgc_tune1-1.0d0),dbl_mb(rho1(1)))
         call D3dB_r_Zero_Ends(1,dbl_mb(rho1(1)))
         call D3dB_rrr_MultiplyAdd(1,dbl_mb(rho1(1)),dbl_mb(tmp1(1)),
     >                             v_out)
      end if
      call D3dB_r_Zero_Ends(1,v_out)

      value =      BA_pop_stack(tmp1(2))
     >        .and.BA_pop_stack(rho1(2))
      if (.not.value) call errquit('wgc_v:popping stack',1,MA_ERR)

      call nwpw_timing_end(7)

      return 
      end


*     ****************************************
*     *                                      *
*     *              wgc_e                   *
*     *                                      *
*     ****************************************
      real*8 function wgc_e(ispin,dn)
      implicit none
      integer ispin
      real*8  dn(*)

#include "bafdecls.fh"
#include "errquit.fh"
#include "wgc.fh"


*     **** local variables ****
      logical value
      integer npack0,nfft3d,n2ft3d
      real*8 ec

      integer tmp1(2),rho1(2),rho2(2)

*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

      call nwpw_timing_start(7)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
      call Pack_npack(0,npack0)

      value =      BA_push_get(mt_dcpl,nfft3d,'rho1',rho1(2),rho1(1))
     >        .and.BA_push_get(mt_dcpl,nfft3d,'rho2',rho2(2),rho2(1))
     >        .and.BA_push_get(mt_dbl, npack0,'tmp1',tmp1(2),tmp1(1))
      if (.not.value) call errquit('wgc_e:out of stack',0,MA_ERR)

      !**** the rho^alpha(G) is not being scaled by 1/(N1*N2*N3) ****
      call D3dB_rr_Sum(1,dn,dn(1+n2ft3d*(ispin-1)),dcpl_mb(rho1(1)))
      call D3dB_r_notZero_Ends(1,dcpl_mb(rho1(1)))
      call D3dB_r_Power1(1,wgc_tune1,dcpl_mb(rho1(1)))
      call D3dB_r_Zero_Ends(1,dcpl_mb(rho1(1)))
      !call D3dB_rc_fft3f(1,dcpl_mb(rho1(1)))
      call D3dB_rc_pfft3f(1,0,dcpl_mb(rho1(1)))
      call Pack_c_pack(0,dcpl_mb(rho1(1)))
      if (wgc_eq_tune) then
         call Pack_ct_Sqr(0,dcpl_mb(rho1(1)),dbl_mb(tmp1(1)))
      else
         call D3dB_rr_Sum(1,dn,dn(1+n2ft3d*(ispin-1)),dcpl_mb(rho2(1)))
         call D3dB_r_notZero_Ends(1,dcpl_mb(rho2(1)))
         call D3dB_r_Power1(1,wgc_tune2,dcpl_mb(rho2(1)))
         call D3dB_r_Zero_Ends(1,dcpl_mb(rho2(1)))
         !call D3dB_rc_fft3f(1,dcpl_mb(rho2(1)))
         call D3dB_rc_pfft3f(1,0,dcpl_mb(rho2(1)))
         call Pack_c_pack(0,dcpl_mb(rho2(1)))
         call Pack_cct_conjgMul(0,dcpl_mb(rho1(1)),
     >                            dcpl_mb(rho2(1)),
     >                            dbl_mb(tmp1(1)))
      end if
      call Pack_tt_dot(0,dbl_mb(tmp1(1)),dbl_mb(wgc_indx),ec)

      !**************************************************************************************************
      !**** Note the energy is really omega*Sum(G) conj(rho^alpha(G))*rho^beta(G)|^2 * WGC_kernel(G) ****
      !**** This funny form is because rho(G) is not scaled 1/(N1*N2*N3) and an                      ****
      !**** extra scaling of 1/(N1*N2*N3) is done to the WGG_kernel(G)                               ****
      !**** Might want to scal rho^alpha/beta(G) properly to make the code more transparent          ****
      !**************************************************************************************************
      ec = ec*lattice_omega()*wgc_scal1
       

      value =      BA_pop_stack(tmp1(2))
     >        .and.BA_pop_stack(rho2(2))
     >        .and.BA_pop_stack(rho1(2))
      if (.not.value) call errquit('wgc_e:popping stack',1,MA_ERR)

      call nwpw_timing_end(7)

      wgc_e = ec
      return 
      end

