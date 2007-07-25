*
* $Id: chi_H.F,v 1.3 2007-07-25 01:29:54 bylaska Exp $
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
*           v_tf                - TF kinetic energy potential in r-space
*           move                - flag to compute ionic forces
*   Exit - Hpsi_k - gradient in k-space
*          fion   - ionic forces
*
      subroutine chi_H(ispin,ne,
     >                 psi_k,psi_r,
     >                 vl,
     >                 v_field,field_exist,
     >                 vc,xcp,v_tf,Hpsi_k,
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
      complex*16 Hpsi_k(*)
      logical    move
      real*8     fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical dosic,dohfx
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      real*8  scal1,scal2

*     **** tempory space ****
c     complex*16 vall(nfft3d)
c     complex*16 vpsi(nfft3d)
      integer nfft3d,n2ft3d,npack1,indx1,indx2
      integer vall(2),vpsi(2),Hpsi_r(2)
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
      n2ft3d = 2*nfft3d
      call nwpw_timing_start(12)
      value = MA_push_get(mt_dcpl,nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
     
      if (.not. value) call errquit('psi_H:out of stack memory',0,
     &       MA_ERR)
      call nwpw_timing_end(12)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)


      call dcopy(2*ispin*npack1,0.0d0,0,Hpsi_k,1)

*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vl + vc  ****
*     **********************************************************
      call Pack_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))
      call Pack_cc_Sum(0,dcpl_mb(vall(1)),vc,dcpl_mb(vall(1)))
      
*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call Pack_c_unpack(0,dcpl_mb(vall(1)))
      call D3dB_cr_fft3b(1,dcpl_mb(vall(1)))

*     ******************************
*     **** add v_field to vall  ****
*     ******************************
      if (field_exist)
     >  call D3dB_rr_Sum(1,dcpl_mb(vall(1)),v_field,dcpl_mb(vall(1)))


*     **************************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*psi_r +Vsic*psi_r] **** 
*     ************************************************************** 
      indx1 = 1
      indx2 = 1
      done  = .false.
      do while (.not.done)
         if (indx1.le.ispin) then
            ms = indx1

*           **** vpsi_r = (vall_r+xcp+v_tf)*psi_r + Hpsi_r****
            call D3dB_rr_Sum(1,dcpl_mb(vall(1)),
     >                       xcp(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))
            call D3dB_rr_Sum(1,dcpl_mb(vpsi(1)),
     >                       v_tf(1+(ms-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))
            call D3dB_rr_Mul(1,dcpl_mb(vpsi(1)),
     >                       psi_r(1+(indx1-1)*n2ft3d),
     >                       dcpl_mb(vpsi(1)))


*           **** put vpsi_r in k-space ****
            call D3dB_rc_pfft3f_queuein(1,dcpl_mb(vpsi(1)))
            indx1 = indx1 + 1

         end if

         if ((D3dB_rc_pfft3_queue_filled()).or.
     >       (indx1.gt.ispin)) then
            call D3dB_rc_pfft3f_queueout(1,dcpl_mb(vpsi(1)))

*           **** add to Hpsi_k ****
            call Pack_cc_daxpy(1,(-scal1),
     >                           dcpl_mb(vpsi(1)),
     >                           Hpsi_k(1+(indx2-1)*npack1))
            indx2 = indx2 + 1
         end if

         done = ((indx1.gt.ispin).and.(indx2.gt.ispin))
      end do
 

      call nwpw_timing_start(12)

      value = .true.
      value = value.and.MA_pop_stack(vpsi(2))
      value = value.and.MA_pop_stack(vall(2))
      if (.not. value) call errquit('chi_H:popping stack memory',1,0)

      call nwpw_timing_end(12)
      return
      end
