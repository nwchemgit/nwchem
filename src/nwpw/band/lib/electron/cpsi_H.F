*
* $Id: cpsi_H.F,v 1.18 2007-11-17 00:05:44 bylaska Exp $
*

*     **********************************************
*     *                                            *
*     *              cpsi_H                        *
*     *                                            *
*     **********************************************

      subroutine cpsi_H(ispin,ne,
     >                 psi_k_tag,psi_r_tag,vl,vc,xcp,Hpsi_k_tag,
     >                 move,fion)
      implicit none
      integer    ispin,ne(2)
      integer    psi_k_tag
      integer    psi_r_tag
      complex*16 vl(*)
      complex*16 vc(*)
      real*8     xcp(*)
      integer    Hpsi_k_tag
      logical    move
      real*8     fion(3,*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      integer psi_shift,Hpsi_shift,nshift1,nshift2
      real*8  scal1,scal2

*     **** tempory space ****
      integer nfft3d,npack1,nb,nbq,neall,nbrill
      integer vall(2),vpsi(2)
      logical value


*     **** external functions ****
      logical  cpsi_spin_orbit
      integer  Pneb_nbrillq,Pneb_convert_nb
      integer  cpsi_data_get_chnk,cpsi_data_nsize
      real*8   lattice_omega
      external cpsi_spin_orbit
      external Pneb_nbrillq,Pneb_convert_nb
      external cpsi_data_get_chnk,cpsi_data_nsize
      external lattice_omega
      
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)

      value = MA_push_get(mt_dcpl,nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)
      nbrill = Pneb_nbrillq()
      nshift1 = cpsi_data_nsize(psi_r_tag)
      nshift2 = cpsi_data_nsize(Hpsi_k_tag)
       

*     ****************************
*     **** add kinetic energy ****
*     ****************************
      call cke(ispin,ne,psi_k_tag,Hpsi_k_tag)

*     **************************
*     **** add non-local psp ***
*     **************************
      call cpsp_v_nonlocal(  ispin,ne,psi_k_tag,Hpsi_k_tag,move,fion)
      call cpsp_v_spin_orbit(ispin,ne,psi_k_tag,Hpsi_k_tag,move,fion)

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

      do nbq=1,nbrill
         nb = Pneb_convert_nb(nbq)
         psi_shift   = cpsi_data_get_chnk(psi_r_tag,nbq)
         Hpsi_shift  = cpsi_data_get_chnk(Hpsi_k_tag,nbq)
         do ms=1,ispin
         do n=n1(ms),n2(ms)

c            indx  = 1 + (n-1)*nfft3d + (nb-1)*nfft3d*neall
c            indx2 = 1 + (n-1)*npack1 + (nb-1)*npack1*neall
*           *************************************
*           **** vpsi_r = (vall_r+xcp)*psi_r ****
*           *************************************
            call C3dB_rc_Sum(1,
     >                       xcp(1+(ms-1)*nfft3d),
     >                       dcpl_mb(vall(1)),
     >                       dcpl_mb(vpsi(1)))
            call C3dB_cc_Mul(1,dcpl_mb(vpsi(1)),
     >                       dbl_mb(psi_shift),
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
     >                         dbl_mb(Hpsi_shift))

            psi_shift  = psi_shift  + nshift1
            Hpsi_shift = Hpsi_shift + nshift2
         end do
         end do
      end do


      value =           MA_pop_stack(vpsi(2))
      value = value.and.MA_pop_stack(vall(2))
      if (.not.value)
     >  call errquit('cpsi_H: error popping stack',0,MA_ERR)

      return
      end


*     *************************************
*     *                                   *
*     *             cpsi_Horb             *
*     *                                   *
*     *************************************
*
*  This routine calculates
*
*       Horb_k = KE*orb_k + Vnl*orb_k + + FFT[(vall+xcp)*orb_r]
*
*  where vall = iFFT[Vl+Vc]
*
*   Entry - ispin,ne            - spin and number of electrons in each spin
*           ms,n                - spin and electron number
*           orb_k,orb_r         - orbitals in k-space and r-space
*           vall                  - local +external fields+coulomb in r-space
*           xcp                 - xc potential in r-space
*   Exit - Horb_k - gradient in k-space
*
      subroutine cpsi_Horb(nb,ispin,ms,
     >                 vall,
     >                 orb_k,orb_r,
     >                 Horb_k)
      implicit none
      integer    nb,ispin,ms
      complex*16 vall(*)
      complex*16 orb_k(*)
      complex*16 orb_r(*)
      complex*16 Horb_k(*)

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer nx,ny,nz,npack1
      real*8  scal1

*     **** tempory space ****
      integer nfft3d,ne1,shiftk
      integer vpsi(2)
      logical value


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

      value = MA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not. value)
     >  call errquit('cpsi_Horb:out of stack memory',0, MA_ERR)

      shiftk=npack1*ne1
      scal1 = 1.0d0/dble(nx*ny*nz)
      
*     ************************************************
*     **** add kinetic energy - Horb_k = Ke*orb_k ****
*     ************************************************
      call cke_orb(nb,orb_k,Horb_k)


*     *********************************************************
*     **** add non-local psp - Horb_k = Horb_k + Vnl*orb_k ****
*     *********************************************************
      call cpsp_v_nonlocal_orb(nb,orb_k,Horb_k)


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

      value = MA_pop_stack(vpsi(2))
      if (.not.value)
     >  call errquit('cpsi_Horb: popping stack',0,MA_ERR)

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
*       Horb_k = KE*orb_k + Vnl*orb_k + + FFT[(vall+xcp)*orb_r]
*
*  where vall = iFFT[Vl+Vc]
*
*   Entry - ispin,ne            - spin and number of electrons in each spin
*           ms,n                - spin and electron number
*           orb_k,orb_r         - orbitals in k-space and r-space
*           vall                  - local +external fields+coulomb in r-space
*           xcp                 - xc potential in r-space
*   Exit - Horb_k - gradient in k-space
*
      subroutine cpsi_Horb2com(nb,ispin,ms,
     >                 vall,
     >                 orb_k,orb_r,
     >                 Horb_k)
      implicit none
      integer    nb,ispin,ms
      complex*16 vall(*)
      complex*16 orb_k(*)
      complex*16 orb_r(*)
      complex*16 Horb_k(*)

#include "mafdecls.fh"
#include "errquit.fh"


*     **** local variables ****
      integer nx,ny,nz
      real*8  scal1

*     **** tempory space ****
      integer nfft3d,npack1,ne1,shiftr,shiftk
      integer vpsi(2)
      logical value


*     **** external functions ****
      real*8   lattice_omega
      integer  cpsi_ne
      external lattice_omega
      external cpsi_ne

      ne1=cpsi_ne(1)
      call Cram_max_npack(npack1)
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call C3dB_nfft3d(1,nfft3d)
      shiftr=ne1*nfft3d + 1
      shiftk=ne1*npack1 + 1

      value = MA_push_get(mt_dcpl,2*nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not. value)
     >  call errquit('cpsi_Horb2com:out of stack memory',0, MA_ERR)


      scal1 = 1.0d0/dble(nx*ny*nz)
      
*     ************************************************
*     **** add kinetic energy - Horb_k = Ke*orb_k ****
*     ************************************************
      call cke_orb(nb,orb_k,Horb_k)


*     *********************************************************
*     **** add non-local psp - Horb_k = Horb_k + Vnl*orb_k ****
*     *********************************************************
      call cpsp_v_nonlocal_orb(nb,orb_k,Horb_k)

 
*     *********************************************************
*     **** add spin_orbit psp - Horb_k = Horb_k + Vso*orb_k ***
*     *********************************************************
      call cpsp_v_spin_orbit_orb(nb,orb_k,Horb_k)


*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*orb_r] ****
*     **************************************************
*     **** vpsi_r = (vall_r)*orb_r ****
      call C3dB_cc_Mul(1,vall(1),
     >                   orb_r,
     >                   dcpl_mb(vpsi(1)))
      call C3dB_cc_Mul(1,vall(1+nfft3d),
     >                   orb_r(shiftr),
     >                   dcpl_mb(vpsi(1)+nfft3d))


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
     >                   Horb_k(shiftk))


      value = MA_pop_stack(vpsi(2))
      if (.not.value)
     >  call errquit('psi_Horb2com: popping stack',1, MA_ERR)

      return
      end



*     ***************************************************
*     *                                                 *
*     *                 cpsi_H_vall                     *
*     *                                                 *
*     ***************************************************
      subroutine cpsi_H_vall(ispin,ne,psi_k_tag,psi_r_tag,vall,
     >                       Hpsi_k_tag)
      implicit none
      integer    ispin,ne(2)
      integer    psi_k_tag
      integer    psi_r_tag
      complex*16 vall(*)
      integer    Hpsi_k_tag

#include "errquit.fh"
#include "mafdecls.fh"

*     **** local variables ****
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      integer psi_shift,Hpsi_shift,nshift1,nshift2
      real*8  scal1,scal2,fion(3)

*     **** tempory space ****
      integer nfft3d,nb,nbq,nbrill
      integer vpsi(2)
      logical value

*     **** external functions ****
      integer  Pneb_nbrillq,Pneb_convert_nb
      integer  cpsi_data_get_chnk,cpsi_data_nsize
      real*8   lattice_omega
      external Pneb_nbrillq,Pneb_convert_nb
      external cpsi_data_get_chnk,cpsi_data_nsize
      external lattice_omega
      !logical  cpsi_spin_orbit
      !external cpsi_spin_orbit
      
      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call C3dB_nfft3d(1,nfft3d)
c      call Cram_max_npack(npack1)

      value = MA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not.value) call errquit('cpsi_H_vall:out of stack',0,MA_ERR)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      nbrill = Pneb_nbrillq()
      nshift1 = cpsi_data_nsize(psi_r_tag)
      nshift2 = cpsi_data_nsize(Hpsi_k_tag)
       

*     ****************************
*     **** add kinetic energy ****
*     ****************************
      call cke(ispin,ne,psi_k_tag,Hpsi_k_tag)

*     **************************
*     **** add non-local psp ***
*     **************************
      call cpsp_v_nonlocal(  ispin,ne,psi_k_tag,Hpsi_k_tag,.false.,fion)
      call cpsp_v_spin_orbit(ispin,ne,psi_k_tag,Hpsi_k_tag,.false.,fion)

*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+)*psi_r]    ****
*     **************************************************
      do nbq=1,nbrill
         nb = Pneb_convert_nb(nbq)
         psi_shift  = cpsi_data_get_chnk( psi_r_tag,nbq)
         Hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nbq)
    
         do ms=1,ispin
         do n=n1(ms),n2(ms)

*           *************************************
*           **** vpsi_r = (vall_r)*psi_r ****
*           *************************************
            call C3dB_cc_Mul(1,vall(1+(ms-1)*nfft3d),
     >                         dbl_mb(psi_shift),
     >                         dcpl_mb(vpsi(1)))

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
     >                         dbl_mb(Hpsi_shift))

            psi_shift  = psi_shift  + nshift1
            Hpsi_shift = Hpsi_shift + nshift2
         end do
         end do
      end do

      value = MA_pop_stack(vpsi(2))
      if(.not.value) call errquit('cpsi_H_vall:popping stack',0,MA_ERR)

      return
      end

