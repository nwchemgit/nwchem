      subroutine psi_Hv4(ispin,ne,
     >                 psi_k,psi_r,vsr_l,vlr_l,vc,xcp,Hpsi_k,
     >                 move,fion)

*
* $Id: psi_Hv4.F,v 1.2 2001-09-29 00:28:07 bylaska Exp $
*
      implicit none
      integer    ispin,ne(2)
      complex*16 psi_k(*)
      real*8     psi_r(*)
      complex*16 vsr_l(*)
      real*8     vlr_l(*)
      real*8     vc(*)
      real*8     xcp(*)
      complex*16 Hpsi_k(*)
      logical    move
      real*8     fion(3,*)

#include "mafdecls.fh"



*     **** local variables ****
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      real*8  scal1,scal2

*     **** tempory space ****
c     complex*16 vall(nfft3d)
c     complex*16 vpsi(nfft3d)
      integer nfft3d,n2ft3d,npack1
      integer vall(2),vpsi(2)
      logical value


*     **** external functions ****
      real*8   lattice_omega
      external lattice_omega

*     ***** timing variables *****
*     1 - fft_time
*     2 - dot_time
*     3 - ABC_time
*     4 - vxc_time
*     5 - vlocal_time
*     6 - vnl_time
*     7 - vh_time
*     8 - strfac_time
*     9 - maskC_time
*    10 - energy_time
*    11 - density_time
*    12 - allocate/deallocate time
      real*8 times(50),tim1,tim2
      common / timing / times


      call Pack_npack(1,npack1)
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call current_second(tim1)
      value = MA_push_get(mt_dcpl,nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not. value) call errquit('out of stack memory',0)
      call current_second(tim2)
      times(12) = times(12) + (tim2-tim1)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)


*     ****************************
*     **** add kinetic energy ****
*     ****************************
      call ke(ispin,ne,psi_k,Hpsi_k)


*     **************************
*     **** add non-local psp ***
*     **************************
      call v_nonlocal(ispin,ne,psi_k,Hpsi_k,move,fion)


*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vsr_l    ****
*     **********************************************************
      call Pack_c_SMul(0,scal2,vsr_l,dcpl_mb(vall(1)))
      

*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vall(1)))
      call D3dB_cr_fft3b(1,dcpl_mb(vall(1)))


*     **********************************
*     **** add vlr_l and vc to vall ****
*     **********************************
      call D3dB_rr_Sum(1,dcpl_mb(vall(1)),vlr_l,dcpl_mb(vall(1)))
      call D3dB_rr_Sum(1,dcpl_mb(vall(1)),vc,   dcpl_mb(vall(1)))


         
      do ms=1,ispin
         do n=n1(ms),n2(ms)

*           *************************************
*           **** vpsi_r = (vall_r+xcp)*psi_r ****
*           *************************************
            call D3dB_rr_Sum(1,dcpl_mb(vall(1)),
     >                       xcp(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))
            call D3dB_rr_Mul(1,dcpl_mb(vpsi(1)),
     >                       psi_r(1+(n-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))

*           *******************************
*           **** get vpsi_r in k-space ****
*           *******************************
            call D3dB_rc_fft3f(1,dcpl_mb(vpsi(1)))
            call Pack_c_pack(1,dcpl_mb(vpsi(1)))

*           ***********************
*           **** add to Hpsi_k ****
*           ***********************
            call Pack_cc_daxpy(1,(-scal1),
     >                           dcpl_mb(vpsi(1)),
     >                           Hpsi_k(1+(n-1)*npack1))


         end do
      end do


      call current_second(tim1)
      value = MA_pop_stack(vpsi(2))
      value = MA_pop_stack(vall(2))
      call current_second(tim2)
      times(12) = times(12) + (tim2-tim1)
      return
      end
