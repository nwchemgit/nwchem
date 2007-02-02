*
* $Id: cpsi_H.F,v 1.12 2007-02-02 01:10:09 d3p708 Exp $
*

      subroutine cpsi_H(ispin,ne,
     >                 psi_k,psi_r,vl,vc,xcp,Hpsi_k,
     >                 move,fion)
      implicit none
#include "errquit.fh"
      integer    ispin,ne(2)
      complex*16 psi_k(*)
      complex*16 psi_r(*)
      complex*16 vl(*)
      complex*16 vc(*)
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
      integer nfft3d,npack1,nb,neall,indx,indx2,nbrill
      integer vall(2),vpsi(2)
      logical value


*     **** external functions ****
      integer  brillioun_nbrillioun
      real*8   lattice_omega
      logical  cpsi_spin_orbit
      external cpsi_spin_orbit
      external brillioun_nbrillioun
      external lattice_omega


      
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      call nwpw_timing_start(12)
      value = MA_push_get(mt_dcpl,nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      call nwpw_timing_end(12)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)
      nbrill = brillioun_nbrillioun()
       

*     ****************************
*     **** add kinetic energy ****
*     ****************************
      call cke(ispin,ne,psi_k,Hpsi_k)

*     **************************
*     **** add non-local psp ***
*     **************************
      call cpsp_v_nonlocal(ispin,ne,psi_k,Hpsi_k,move,fion)
      call cpsp_v_spin_orbit(ispin,ne,psi_k,Hpsi_k,move,fion)

*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vl + vc  ****
*     **********************************************************
      call Cram_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))
      call Cram_cc_Sum(0,dcpl_mb(vall(1)),vc,dcpl_mb(vall(1)))
      call Cram_c_unpack(0,dcpl_mb(vall(1)))

*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call C3dB_cr_fft3b(1,dcpl_mb(vall(1)))

      do nb=1,nbrill
      do ms=1,ispin
         do n=n1(ms),n2(ms)

            indx  = 1 + (n-1)*nfft3d + (nb-1)*nfft3d*neall
            indx2 = 1 + (n-1)*npack1 + (nb-1)*npack1*neall
*           *************************************
*           **** vpsi_r = (vall_r+xcp)*psi_r ****
*           *************************************
            call C3dB_rc_Sum(1,
     >                       xcp(1+(ms-1)*nfft3d),
     >                       dcpl_mb(vall(1)),
     >                       dcpl_mb(vpsi(1)))
            call C3dB_cc_Mul(1,dcpl_mb(vpsi(1)),
     >                       psi_r(indx),
     >                       dcpl_mb(vpsi(1)))

*           *******************************
*           **** get vpsi_r in k-space ****
*           *******************************
            call C3dB_rc_fft3f(1,dcpl_mb(vpsi(1)))
            call Cram_c_pack(nb,dcpl_mb(vpsi(1)))

*           ***********************
*           **** add to Hpsi_k ****
*           ***********************

            call Cram_cc_daxpy(nb,(-scal1),
     >                         dcpl_mb(vpsi(1)),
     >                         Hpsi_k(indx2))
         end do
      end do
      end do


      call nwpw_timing_start(12)
      
      value = MA_pop_stack(vpsi(2))
      value = MA_pop_stack(vall(2))

      call nwpw_timing_end(12)
      return
      end




*     *************************************
*     *                                   *
*     *             cpsi_Horb              *
*     *                                   *
*     *************************************
*
*  This routine calculates
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k + + FFT[(vall+xcp)*psi_r]
*
*  where vall = iFFT[Vl+Vc]
*
*   Entry - ispin,ne            - spin and number of electrons in each spin
*           ms,n                - spin and electron number
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vall                  - local +external fields+coulomb in r-space
*           xcp                 - xc potential in r-space
*   Exit - Hpsi_k - gradient in k-space
*
      subroutine cpsi_Horb(nb,ispin,ms,
     >                 psi_r,
     >                 vall,
     >                 orb_k,orb_r,
     >                 Horb_k)
      implicit none
      integer    nb,ispin,ms
      complex*16 psi_r(*)
      complex*16 vall(*)
      complex*16 orb_k(*)
      complex*16 orb_r(*)
      complex*16 Horb_k(*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer nne(2)
      integer nx,ny,nz,npack1
      real*8  scal1,fion(3,1)

*     **** tempory space ****
      integer nfft3d,ne1,shiftk
      integer vpsi(2),spsi(2),tpsi(2)
      logical value
      complex*16 tmp


*     **** external functions ****
      integer cpsi_ne
      real*8   lattice_omega
      logical  control_spin_orbit
      complex*16 zdotc
      external lattice_omega
      external control_spin_orbit
      external zdotc,cpsi_ne

      if (control_spin_orbit()) then
        call cpsi_Horb2com(nb,ispin,ms,
     >                 psi_r,
     >                 vall,
     >                 orb_k,orb_r,
     >                 Horb_k)
        return
      end if

      ne1=cpsi_ne(1)
      call Cram_max_npack(npack1)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call C3dB_nfft3d(1,nfft3d)
      call nwpw_timing_start(12)
      value = MA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
c      value = value.and.MA_push_get(mt_dcpl,2*npack1*ne1,'Torb_k',
c     > tpsi(2),tpsi(1))
c      value = value.and.MA_push_get(mt_dcpl,2*nfft3d*ne1,'Sorb_k',
c     > spsi(2),spsi(1))
c      if (.not. value)
c     >  call errquit('cpsi_Horb:out of stack memory',0, MA_ERR)
      call nwpw_timing_end(12)

      shiftk=npack1*ne1
      scal1 = 1.0d0/dble(nx*ny*nz)
      
*     ************************************************
*     **** add kinetic energy - Horb_k = Ke*orb_k ****
*     ************************************************
      call cke_orb(nb,orb_k,Horb_k)

c      call dcopy(2*shiftk,0.0,0,dcpl_mb(tpsi(1)),1)
c      call cke_orb(nb,orb_k,dcpl_mb(tpsi(1)))
c      tmp=zdotc(npack1,orb_k,1,dcpl_mb(tpsi(1)),1)
c      write(*,*)"ke",tmp
c      call flush(6)

*     *********************************************************
*     **** add non-local psp - Horb_k = Horb_k + Vnl*orb_k ****
*     *********************************************************
      call cpsp_v_nonlocal_orb(nb,orb_k,Horb_k)

c      call dcopy(2*shiftk,0.0,0,dcpl_mb(tpsi(1)),1)
c      call cpsp_v_nonlocal_orb(nb,orb_k,dcpl_mb(tpsi(1)))
c      tmp=zdotc(npack1,orb_k,1,dcpl_mb(tpsi(1)),1)
c      write(*,*)"vnl",tmp
c      call flush(6)

*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*orb_r] ****
*     **************************************************

*     **** vpsi_r = (vall_r)*orb_r ****
      call C3dB_cc_Mul(1,vall(1+(ms-1)*nfft3d),
     >                   orb_r,
     >                   dcpl_mb(vpsi(1)))


*     **** put vpsi_r in k-space ****
      call C3dB_rc_fft3f(1,dcpl_mb(vpsi(1)))
      call Cram_c_pack(nb,dcpl_mb(vpsi(1)))

*     **** add to Horb_k ****
      call Cram_cc_daxpy(nb,(-scal1),
     >                   dcpl_mb(vpsi(1)),
     >                   Horb_k)

c      call dcopy(2*shiftk,0.0,0,dcpl_mb(tpsi(1)),1)
c      call Cram_cc_daxpy(nb,(-scal1),
c     >                   dcpl_mb(vpsi(1)),
c     >                   dcpl_mb(tpsi(1)))
c      tmp=zdotc(npack1,orb_k,1,dcpl_mb(tpsi(1)),1)
c      write(*,*)"vall",tmp
c      call flush(6)

      call nwpw_timing_start(12)
      value = .true.
c      value = value.and.MA_pop_stack(spsi(2))
c      value = value.and.MA_pop_stack(tpsi(2))
      value = value.and.MA_pop_stack(vpsi(2))
      if (.not.value)
     >  call errquit('psi_Horb: popping stack',1, MA_ERR)
      call nwpw_timing_end(12)

      return
      end

*     *************************************
*     *                                   *
*     *             cpsi_Horb2com         *
*     *                                   *
*     *************************************
*
*  This routine calculates
*
*       Hpsi_k = KE*psi_k + Vnl*psi_k + + FFT[(vall+xcp)*psi_r]
*
*  where vall = iFFT[Vl+Vc]
*
*   Entry - ispin,ne            - spin and number of electrons in each spin
*           ms,n                - spin and electron number
*           psi_k,psi_r         - orbitals in k-space and r-space
*           vall                  - local +external fields+coulomb in r-space
*           xcp                 - xc potential in r-space
*   Exit - Hpsi_k - gradient in k-space
*
      subroutine cpsi_Horb2com(nb,ispin,ms,
     >                 psi_r,
     >                 vall,
     >                 orb_k,orb_r,
     >                 Horb_k)
      implicit none
      integer    nb,ispin,ms
      complex*16 psi_r(*)
      complex*16 vall(*)
      complex*16 orb_k(*)
      complex*16 orb_r(*)
      complex*16 Horb_k(*)

#include "mafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      integer nne(2)
      integer nx,ny,nz
      real*8  scal1,fion(3,1)

*     **** tempory space ****
      integer nfft3d,npack1,ne1,shiftr,shiftk
      integer vpsi(2),tpsi(2),spsi(2)
      complex*16 tmp,tmp2
      logical value


*     **** external functions ****
      real*8   lattice_omega
      complex*16 zdotc
      integer cpsi_ne
      external lattice_omega
      external zdotc
      external cpsi_ne

      
c      value=MA_set_auto_verify(.true.)
      ne1=cpsi_ne(1)
      call Cram_max_npack(npack1)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call C3dB_nfft3d(1,nfft3d)
      shiftr=ne1*nfft3d
      shiftk=ne1*npack1
      call nwpw_timing_start(12)
      value = MA_push_get(mt_dcpl,2*nfft3d,'vpsi',vpsi(2),vpsi(1))
c      value = value.and.MA_push_get(mt_dcpl,2*shiftk,'Torb_k',
c     > tpsi(2),tpsi(1))
c      value = value.and.MA_push_get(mt_dcpl,2*shiftr,'Sorb_k',
c     > spsi(2),spsi(1))
      if (.not. value)
     >  call errquit('cpsi_Horb2com:out of stack memory',0, MA_ERR)
      call nwpw_timing_end(12)


      scal1 = 1.0d0/dble(nx*ny*nz)
      
*     ************************************************
*     **** add kinetic energy - Horb_k = Ke*orb_k ****
*     ************************************************
      call cke_orb(nb,orb_k,Horb_k)

c      call dcopy(4*shiftk,0.0,0,dcpl_mb(tpsi(1)),1)
c      call cke_orb(nb,orb_k,dcpl_mb(tpsi(1)))
c      tmp=zdotc(npack1,orb_k,1,dcpl_mb(tpsi(1)),1)
c      tmp=tmp+zdotc(npack1,orb_k(shiftk+1),
c     >           1,dcpl_mb(tpsi(1)+shiftk),1)   
c      value=MA_verify_allocator_stuff()
c      write(*,*)"ke",tmp
c      call flush(6)

*     *********************************************************
*     **** add non-local psp - Horb_k = Horb_k + Vnl*orb_k ****
*     *********************************************************
      call cpsp_v_nonlocal_orb(nb,orb_k,Horb_k)

c      call dcopy(4*shiftk,0.0,0,dcpl_mb(tpsi(1)),1)
c      call cpsp_v_nonlocal_orb(nb,orb_k,dcpl_mb(tpsi(1)))
c      tmp=zdotc(npack1,orb_k,1,dcpl_mb(tpsi(1)),1)
c      tmp=tmp+zdotc(npack1,orb_k(shiftk+1),
c     >           1,dcpl_mb(tpsi(1)+shiftk),1)   
c      value=MA_verify_allocator_stuff()
c      write(*,*)"vnl",tmp
c      call flush(6)
 
*     *********************************************************
*     **** add spin_orbit psp - Horb_k = Horb_k + Vso*orb_k ***
*     *********************************************************
      call cpsp_v_spin_orbit_orb(nb,orb_k,Horb_k)

c      call dcopy(4*shiftk,0.0,0,dcpl_mb(tpsi(1)),1)
c      call cpsp_v_spin_orbit_orb(nb,orb_k,dcpl_mb(tpsi(1)))
c      tmp=zdotc(npack1,orb_k,1,dcpl_mb(tpsi(1)),1)
c      tmp=tmp+zdotc(npack1,orb_k(shiftk+1),
c     >           1,dcpl_mb(tpsi(1)+shiftk),1)   
c      value=MA_verify_allocator_stuff()
c      write(*,*)"vnlso",tmp
c      call flush(6)

*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*orb_r] ****
*     **************************************************
*     **** vpsi_r = (vall_r)*orb_r ****
      call C3dB_cc_Mul(1,vall(1),
     >                   orb_r,
     >                   dcpl_mb(vpsi(1)))
      call C3dB_cc_Mul(1,vall(1+nfft3d),
     >                   orb_r(shiftr+1),
     >                   dcpl_mb(vpsi(1)+nfft3d))

      value=MA_verify_allocator_stuff()
*     **** put vpsi_r in k-space ****
      call C3dB_rc_fft3f(1,dcpl_mb(vpsi(1)))
      call Cram_c_pack(nb,dcpl_mb(vpsi(1)))

      call C3dB_rc_fft3f(1,dcpl_mb(vpsi(1)+nfft3d))
      call Cram_c_pack(nb,dcpl_mb(vpsi(1)+nfft3d))

*     **** add to Horb_k ****
      call Cram_cc_daxpy(nb,(-scal1),
     >                   dcpl_mb(vpsi(1)),
     >                   Horb_k)
      call Cram_cc_daxpy(nb,(-scal1),
     >                   dcpl_mb(vpsi(1)+nfft3d),
     >                   Horb_k(shiftk+1))


c      call dcopy(4*shiftk,0.0,0,dcpl_mb(tpsi(1)),1)
c      call Cram_cc_daxpy(nb,(-scal1),
c     >                   dcpl_mb(vpsi(1)),
c     >                   dcpl_mb(tpsi(1)))
c      call Cram_cc_daxpy(nb,(-scal1),
c     >                   dcpl_mb(vpsi(1)+nfft3d),
c     >                   dcpl_mb(tpsi(1)+shiftk))
c      tmp=zdotc(npack1,orb_k,1,dcpl_mb(tpsi(1)),1)
c      tmp=tmp+zdotc(npack1,orb_k(shiftk+1),
c     >           1,dcpl_mb(tpsi(1)+shiftk),1)   
c      value=MA_verify_allocator_stuff()
c      write(*,*)"vall",tmp
c      call flush(6)
      
c      value=MA_verify_allocator_stuff()
      call nwpw_timing_start(12)
      value = .true.
c      value = value.and.MA_pop_stack(spsi(2))
c      value = value.and.MA_pop_stack(tpsi(2))
      value = value.and.MA_pop_stack(vpsi(2))
      if (.not.value)
     >  call errquit('psi_Horb2com: popping stack',1, MA_ERR)
      call nwpw_timing_end(12)
      call nwpw_timing_start(12)

      return
      end


