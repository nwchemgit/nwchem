*
* $Id$
*

*     *************************************
*     *                                   *
*     *             chi_H                 *
*     *                                   *
*     *************************************
*
*  This routine calculates 
*
*       Hpsi_k = vke2*psi_k +  + FFT[(vke+vall+xcp)*psi_r]
*
*  where vall = iFFT[Vl+Vc] + Vfield
*
*   Entry - ispin,ne            - spin and number of electron in each spin
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vl                  - local psp in k-space
*           v_field,field_exist - external fields (e.g. QM/MM)
*           vc                  - coulomb potential in k-space
*           xcp                 - xc potential in r-space
*           v_tf                - TF kinetic potential in r-space
*           v_wgc               - WGC kinetic potential in r-space
*           move                - flag to compute ionic forces
*   Exit - Hpsi_k - gradient in k-space
*          fion   - ionic forces
*
      subroutine chi_H(ispin,ne,
     >                 psi_k,psi_r,
     >                 vl,
     >                 v_field,field_exist,
     >                 vc,xcp,v_tf,v_wgc,Hpsi_k,
     >                 move,fion)
      implicit none
      integer    ispin,ne(2)
      complex*16 psi_k(*)
      real*8     psi_r(*)
      complex*16 vl(*)
      real*8     v_field(*)
      logical    field_exist
      complex*16 vc(*)
      real*8     xcp(*)
      real*8     v_tf(*)
      real*8     v_wgc(*)
      complex*16 Hpsi_k(*)
      logical    move
      real*8     fion(3,*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical dosic,dohfx
      integer ms,n,n1(2),n2(2),nne(2)
      integer nx,ny,nz
      real*8  scal1,scal2

*     **** tempory space ****
c     complex*16 vall(nfft3d)
c     complex*16 vpsi(nfft3d)
      integer nfft3d,n2ft3d,npack1,indx1,indx2
      integer vall(2),vpsi(2)
      logical value,done


*     **** external functions ****
      logical  D3dB_rc_pfft3_queue_filled
      real*8   lattice_omega
      external D3dB_rc_pfft3_queue_filled
      external lattice_omega


      call Pack_npack(1,npack1)
      
      call D3dB_nx(1,nx)
      call D3dB_ny(1,ny)
      call D3dB_nz(1,nz)
      call D3dB_nfft3d(1,nfft3d)
      call D3dB_n2ft3d(1,n2ft3d)
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



*     *********************************************************
*     **** add kinetic energy - Hpsi_k = Ke*psi_k ****
*     *********************************************************
      !call dcopy(2*ispin*npack1,0.0d0,0,Hpsi_k,1)
      nne(1) = 1
      nne(2) = ispin-1
      call ke(ispin,nne,psi_k,Hpsi_k)
      do ms=1,ispin
         call dscal(2*npack1,dble(ne(ms)),Hpsi_k(1+(ms-1)*npack1),1)
      end do


*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vl + vc  ****
*     **********************************************************
      call Pack_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))
      call Pack_cc_Sum2(0,vc,dcpl_mb(vall(1)))
      
*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vall(1)))
      call D3dB_cr_pfft3b(1,0,dcpl_mb(vall(1)))

*     ******************************
*     **** add v_field to vall  ****
*     ******************************
      if (field_exist)
     >  call D3dB_rr_Sum(1,dcpl_mb(vall(1)),v_field,dcpl_mb(vall(1)))


*     ****************************
*     **** add v_wgc to vall  ****
*     ****************************
      call D3dB_rr_Sum2(1,v_wgc,dcpl_mb(vall(1)))


*     **************************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*psi_r +Vsic*psi_r] **** 
*     ************************************************************** 

      do ms=1,ispin
         call D3dB_rr_Sum(1,dcpl_mb(vall(1)),
     >                      xcp(1+(ms-1)*n2ft3d),
     >                      dcpl_mb(vpsi(1)))
         call D3dB_rr_Sum2(1,v_tf(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))
         call D3dB_rr_Mul2(1,
     >                     psi_r(1+(ms-1)*n2ft3d),
     >                     dcpl_mb(vpsi(1)))
         call D3dB_rc_pfft3f(1,1,dcpl_mb(vpsi(1)))
         call Pack_c_pack(1,dcpl_mb(vpsi(1)))

*        **** add to Hpsi_k ****
         call Pack_cc_daxpy(1,(-scal1)*dble(ne(ms)),
     >                      dcpl_mb(vpsi(1)),
     >                      Hpsi_k(1+(ms-1)*npack1))
      end do
 
      call nwpw_timing_start(12)

      value = .true.
      value = value.and.BA_pop_stack(vpsi(2))
      value = value.and.BA_pop_stack(vall(2))
      if (.not. value) call errquit('chi_H:popping stack memory',1,0)

      call nwpw_timing_end(12)
      return
      end
