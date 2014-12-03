*
* $Id$
*

*     *************************************
*     *                                   *
*     *          paw_psi_H                *
*     *                                   *
*     *************************************
*
*  This routine calculates 
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k +VSic*psi_k + FFT[(vall+xcp)*psi_r]
*
*  where vall = iFFT[Vl+Vc] + Vfield
*
*   Entry - ispin,ne            - spin and number of electron in each spin
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vl                  - local psp in k-space
*           vc                  - coulomb potential in k-space
*           xcp                 - xc potential in r-space
*           move                - flag to compute ionic forces
*   Exit - Hpsi_k - gradient in k-space
*          fion   - ionic forces
*
      subroutine paw_psi_H(ispin,ne,
     >                 psi_k,psi_r,
     >                 vl,
     >                 vc,xcp,Hpsi_k,
     >                 move,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi_k(*)
      real*8     psi_r(*)
      complex*16 vl(*)
      complex*16 vc(*)
      real*8     xcp(*)
      complex*16 Hpsi_k(*)
      logical    move
      real*8     fion(3,*)

#include "bafdecls.fh"

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
      if (.not. value) call errquit('out of stack memory',0,0)
      call nwpw_timing_end(12)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)


*     ************************************************
*     **** add kinetic energy - Hpsi_k = Ke*psi_k ****
*     ************************************************
      call ke(ispin,ne,psi_k,Hpsi_k)


*     *********************************************************
*     **** add non-local psp - Hpsi_k = Hpsi_k + Vnl*psi_k ****
*     *********************************************************
      call paw_nonlocal(ispin,ne,psi_k,Hpsi_k,move,fion)


*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vl + vc  ****
*     **********************************************************
      call Pack_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))
c      call Pack_cc_Sum(0,dcpl_mb(vall(1)),vc,dcpl_mb(vall(1)))
      call Pack_cc_Sum2(0,vc,dcpl_mb(vall(1)))
      
*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vall(1)))
      call D3dB_cr_fft3b(1,dcpl_mb(vall(1)))



*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*psi_r] **** 
*     **************************************************        
      do ms=1,ispin
         do n=n1(ms),n2(ms)

*           **** vpsi_r = (vall_r+xcp)*psi_r ****
            call D3dB_rr_Sum(1,dcpl_mb(vall(1)),
     >                       xcp(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))
c            call D3dB_rr_Mul(1,dcpl_mb(vpsi(1)),
c     >                       psi_r(1+(n-1)*n2ft3d),
c     >                       dcpl_mb(vpsi(1)))
            call D3dB_rr_Mul2(1,psi_r(1+(n-1)*n2ft3d),
     >                          dcpl_mb(vpsi(1)))

*           **** put vpsi_r in k-space ****
            call D3dB_rc_fft3f(1,dcpl_mb(vpsi(1)))
            call Pack_c_pack(1,dcpl_mb(vpsi(1)))

*           **** add to Hpsi_k ****
            call Pack_cc_daxpy(1,(-scal1),
     >                           dcpl_mb(vpsi(1)),
     >                           Hpsi_k(1+(n-1)*npack1))
         end do
      end do

      call nwpw_timing_start(12)
      value = BA_pop_stack(vpsi(2))
      value = BA_pop_stack(vall(2))
      call nwpw_timing_end(12)
      return
      end




*     *************************************
*     *                                   *
*     *         paw_psi_Hv4               *
*     *                                   *
*     *************************************
*
*  This routine calculates
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k + FFT[(vall+xcp)*psi_r]
*
*  where vall = iFFT[Vsr_l]+ vc 
*
*  Note - the difference between psi_Hv4 and psi_H is that psi_Hv4 is designed
*         for computing vc*psi using aperiodic convolution.
*
*   Entry - ispin,ne            - spin and number of electron in each spin
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vl               - short-range local psp in k-space
*           vc                  - coulomb potential in k-space
*           xcp                 - xc potential in r-space
*           move                - flag to compute ionic forces
*   Exit - Hpsi_k - gradient in k-space
*          fion   - ionic forces
*

      subroutine paw_psi_Hv4(ispin,ne,
     >                 psi_k,psi_r,
     >                 vl,
     >                 vc,xcp,Hpsi_k,
     >                 move,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi_k(*)
      real*8     psi_r(*)
      complex*16 vl(*)
      real*8     vc(*)
      real*8     xcp(*)
      complex*16 Hpsi_k(*)
      logical    move
      real*8     fion(3,*)

#include "bafdecls.fh"


*     **** local variables ****
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      real*8  scal1,scal2

*     **** tempory space ****
      integer nfft3d,n2ft3d,npack1
      integer vall(2),vpsi(2)
      logical value


*     **** external functions ****
      real*8   lattice_omega
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
      if (.not. value) call errquit('out of stack memory',0,0)
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
      call ke(ispin,ne,psi_k,Hpsi_k)

*     **************************
*     **** add non-local psp ***
*     **************************
c      call v_nonlocal(ispin,ne,psi_k,Hpsi_k,move,fion)


*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vsr_l    ****
*     **********************************************************
      call Pack_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))

*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vall(1)))
      call D3dB_cr_fft3b(1,dcpl_mb(vall(1)))

*     *******************************************
*     ****    add vc and v_field to vall     ****
*     *******************************************
c      call D3dB_rr_Sum(1,dcpl_mb(vall(1)),vc,dcpl_mb(vall(1)))
      call D3dB_rr_Sum2(1,vc,dcpl_mb(vall(1)))

*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*psi_r] ****
*     **************************************************
      do ms=1,ispin
         do n=n1(ms),n2(ms)

*           **** vpsi_r = (vall_r+xcp)*psi_r ****
            call D3dB_rr_Sum(1,dcpl_mb(vall(1)),
     >                       xcp(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))
c            call D3dB_rr_Mul(1,dcpl_mb(vpsi(1)),
c     >                       psi_r(1+(n-1)*n2ft3d),
c     >                       dcpl_mb(vpsi(1)))
            call D3dB_rr_Mul2(1,psi_r(1+(n-1)*n2ft3d),
     >                          dcpl_mb(vpsi(1)))

*           **** get vpsi_r in k-space ****
            call D3dB_rc_fft3f(1,dcpl_mb(vpsi(1)))
            call Pack_c_pack(1,dcpl_mb(vpsi(1)))

*           **** add to Hpsi_k ****
            call Pack_cc_daxpy(1,(-scal1),
     >                           dcpl_mb(vpsi(1)),
     >                           Hpsi_k(1+(n-1)*npack1))
         end do
      end do

      call nwpw_timing_start(12)
      value = BA_pop_stack(vpsi(2))
      value = BA_pop_stack(vall(2))
      call nwpw_timing_end(12)
      return
      end


