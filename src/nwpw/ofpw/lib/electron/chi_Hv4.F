*
* $Id$
*

*     *************************************
*     *                                   *
*     *             chi_Hv4               *
*     *                                   *
*     *************************************
*
*  This routine calculates 
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k +VSic*psi_k + FFT[(vall+xcp)*psi_r]
*
*  where vall = iFFT[Vsr_l]+ vc + vlr_l + Vfield
*
*  Note - the difference between psi_Hv4 and psi_H is that psi_Hv4 is designed
*         for computing vc*psi using aperiodic convolution.
*
*   Entry - ispin,ne            - spin and number of electron in each spin
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vsr_l               - short-range local psp in k-space
*           vlr_l               - long-range local psp in r-space
*           v_field,field_exist - external fields (e.g. QM/MM)
*           vc                  - coulomb potential in k-space
*           xcp                 - xc potential in r-space
*           move                - flag to compute ionic forces
*   Exit - Hpsi_k - gradient in k-space
*          fion   - ionic forces
*

      subroutine chi_Hv4(ispin,ne,
     >                 psi_k,psi_r,
     >                 vsr_l,vlr_l,
     >                 v_field,field_exist,
     >                 vc,xcp,v_tf,v_wgc,Hpsi_k,
     >                 move,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi_k(*)
      real*8     psi_r(*)
      complex*16 vsr_l(*)
      real*8     vlr_l(*)
      real*8     v_field(*)
      logical    field_exist
      real*8     vc(*)
      real*8     xcp(*)
      real*8     v_tf(*)
      real*8     v_wgc(*)
      complex*16 Hpsi_k(*)
      logical    move
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer ms,n,n1(2),n2(2),nne(2)
      integer nx,ny,nz
      real*8  scal1,scal2

*     **** tempory space ****
      integer nfft3d,n2ft3d,npack1
      integer vall(2),vpsi(2),Hpsi_r(2)
      logical value


*     **** external functions ****
      logical  pspw_SIC,pspw_SIC_relaxed
      logical  pspw_HFX,pspw_HFX_relaxed
      real*8   lattice_omega 
      external pspw_SIC,pspw_SIC_relaxed
      external pspw_HFX,pspw_HFX_relaxed
      external lattice_omega


      call Pack_npack(1,npack1)
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nfft3d(1,nfft3d)
      n2ft3d = 2*nfft3d
      call nwpw_timing_start(12)
      value = BA_push_get(mt_dcpl,nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not. value) call errquit('psi_H:out of stack memory',0,
     &       MA_ERR)
      call nwpw_timing_end(12)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)


*     ****************************
*     **** add kinetic energy ****
*     ****************************
      !call ke(ispin,ne,psi_k,Hpsi_k)
      nne(1) = 1
      nne(2) = ispin-1
      call ke(ispin,nne,psi_k,Hpsi_k)



*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vsr_l    ****
*     **********************************************************
      call Pack_c_SMul(0,scal2,vsr_l,dcpl_mb(vall(1)))
      
*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vall(1)))
      call D3dB_cr_fft3b(1,dcpl_mb(vall(1)))

*     *******************************************
*     **** add vlr_l, vc and v_field to vall ****
*     *******************************************
      call D3dB_rr_Sum(1,dcpl_mb(vall(1)),vlr_l,dcpl_mb(vall(1)))
      call D3dB_rr_Sum(1,dcpl_mb(vall(1)),vc,   dcpl_mb(vall(1)))
      call D3dB_rr_Sum(1,dcpl_mb(vall(1)),v_wgc,dcpl_mb(vall(1)))
      if (field_exist)
     >  call D3dB_rr_Sum(1,dcpl_mb(vall(1)),v_field,dcpl_mb(vall(1)))


*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*psi_r] **** 
*     **************************************************       
      do ms=1,ispin

*           **** vpsi_r = (vall_r+xcp)*psi_r ****
            call D3dB_rr_Sum(1,dcpl_mb(vall(1)),
     >                       xcp(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))
            call D3dB_rr_Sum2(1,v_tf(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))
            call D3dB_rr_Mul2(1,psi_r(1+(ms-1)*n2ft3d),
     >                          dcpl_mb(vpsi(1)))

*           **** get vpsi_r in k-space ****
            call D3dB_rc_pfft3f(1,1,dcpl_mb(vpsi(1)))
            call Pack_c_pack(1,dcpl_mb(vpsi(1)))

*           **** add to Hpsi_k ****
            call Pack_cc_daxpy(1,(-scal1)*dble(ne(ms)),
     >                           dcpl_mb(vpsi(1)),
     >                           Hpsi_k(1+(ms-1)*npack1))
      end do

      call nwpw_timing_start(12)

      value = .true.
      value = value.and.BA_pop_stack(vpsi(2))
      value = value.and.BA_pop_stack(vall(2))
      if (.not. value) call errquit('psi_H:popping stack memory',1,0)

      call nwpw_timing_end(12)
      return
      end

