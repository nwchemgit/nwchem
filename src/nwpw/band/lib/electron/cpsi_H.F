*
* $Id$
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

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical done,dohfx
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      integer ki,kj,nbi,nbj
      integer psi_shift,Hpsi_shift,Hpsir_shift,nshift1,nshift2
      real*8  scal1,scal2

*     **** tempory space ****
      integer nfft3d,npack0,nb,nbq,neall,nbrillq
      integer vall(2),vpsi(2),Hpsi_r_tag
      logical value


*     **** external functions ****
      logical  C3dB_rc_pfft3_queue_filled,band_HFX,band_HFX_relaxed
      integer  Pneb_nbrillq,Pneb_convert_nb
      integer  cpsi_data_get_chnk,cpsi_data_nsize,cpsi_data_alloc
      real*8   lattice_omega
      external C3dB_rc_pfft3_queue_filled,band_HFX,band_HFX_relaxed
      external Pneb_nbrillq,Pneb_convert_nb
      external cpsi_data_get_chnk,cpsi_data_nsize,cpsi_data_alloc
      external lattice_omega

      dohfx = (band_HFX().and.band_HFX_relaxed())

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call C3dB_nfft3d(1,nfft3d)
      call Cram_npack(0,npack0)

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      neall = ne(1) + ne(2)
      nbrillq = Pneb_nbrillq()
      nshift1 = cpsi_data_nsize(psi_r_tag)
      nshift2 = cpsi_data_nsize(Hpsi_k_tag)

      value = BA_push_get(mt_dcpl,nfft3d,'vall',vall(2),vall(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not. value) call errquit('out of stack memory',0, MA_ERR)

      if (dohfx) then
         Hpsi_r_tag = cpsi_data_alloc(nbrillq,ne(1)+ne(2),2*nfft3d)
      end if
       

*     ****************************
*     **** add kinetic energy ****
*     ****************************
      call cke(ispin,ne,psi_k_tag,Hpsi_k_tag)

*     **************************
*     **** add non-local psp ***
*     **************************
      call cpsp_v_nonlocal(ispin,ne,psi_k_tag,Hpsi_k_tag,move,fion)
      call cpsp_v_spin_orbit(ispin,ne,psi_k_tag,Hpsi_k_tag,
     >  move,fion)

*     *******************************************************
*     **** HF exchange corrections- Hpsi_r = Vhfx*psi_r ****
*     *******************************************************
      if (dohfx) then
         call band_potential_HFX(ispin,psi_r_tag,Hpsi_r_tag)
      end if



*     **********************************************************
*     **** add up k-space potentials, vall = scal2*vl + vc  ****
*     **********************************************************
c      call Cram_c_SMul(0,scal2,vl,dcpl_mb(vall(1)))
c      call Cram_cc_Sum(0,dcpl_mb(vall(1)),vc,dcpl_mb(vall(1)))
      call dcopy(2*npack0,vc,1,dcpl_mb(vall(1)),1)
      call daxpy(2*npack0,scal2,vl,1,dcpl_mb(vall(1)),1)
      call Cram_c_unpack(0,dcpl_mb(vall(1)))

*     **********************************************
*     **** fourier transform k-space potentials ****
*     **********************************************
      call C3dB_cr_pfft3b(1,0,dcpl_mb(vall(1)))

      ki  = 1
      nbi = 1 
      psi_shift   = cpsi_data_get_chnk(psi_r_tag,nbi)
      if (dohfx) Hpsir_shift = cpsi_data_get_chnk(Hpsi_r_tag,nbi)

      kj  = 1
      nbj = 1 
      Hpsi_shift  = cpsi_data_get_chnk(Hpsi_k_tag,nbj)
      done = .false.
      do while (.not.done)
         if (nbi.le.nbrillq) then
            ms=1
            if (ki.gt.n2(1)) ms = 2

*           **** vpsi_r = (vall_r+xcp)*psi_r ****
            call C3dB_rc_Sum(1,
     >                       xcp(1+(ms-1)*nfft3d),
     >                       dcpl_mb(vall(1)),
     >                       dcpl_mb(vpsi(1)))
c            call C3dB_cc_Mul(1,dcpl_mb(vpsi(1)),
c     >                       dbl_mb(psi_shift),
c     >                       dcpl_mb(vpsi(1)))
            call C3dB_bb_Mul2c(1,dbl_mb(psi_shift),
     >                           dcpl_mb(vpsi(1)))

            if (dohfx)
     >         call C3dB_bb_Sum2(1,dbl_mb(Hpsir_shift),dcpl_mb(vpsi(1)))

*           **** get vpsi_r in k-space ****
            call C3dB_rc_pfft3f_queuein(nbi,dcpl_mb(vpsi(1)))

            psi_shift   = psi_shift   + nshift1
            if (dohfx) Hpsir_shift = Hpsir_shift + nshift1
            ki         = ki + 1
            if (ki.gt.neall) then
               nbi = nbi + 1
               ki  = 1
               if (nbi.le.nbrillq) then
                  psi_shift   = cpsi_data_get_chnk(psi_r_tag, nbi)
                  if (dohfx) 
     >               Hpsir_shift = cpsi_data_get_chnk(Hpsi_r_tag,nbi)
               endif
            end if
         end if

         if ((C3dB_rc_pfft3_queue_filled()).or.(nbi.gt.nbrillq)) then
            call C3dB_rc_pfft3f_queueout(nbj,dcpl_mb(vpsi(1)))

*           **** add to Hpsi_k ****
            call Cram_cc_daxpy(nbj,(-scal1),
     >                         dcpl_mb(vpsi(1)),
     >                         dbl_mb(Hpsi_shift))

            Hpsi_shift = Hpsi_shift + nshift2
            kj         = kj + 1
            if (kj.gt.neall) then
               nbj = nbj + 1
               kj  = 1
               if (nbj.le.nbrillq) then
                  Hpsi_shift   = cpsi_data_get_chnk(Hpsi_k_tag,nbj)
               endif
            end if

         end if

         done = ((nbi.gt.nbrillq).and.(nbj.gt.nbrillq))
      end do


      value = .true.
      if (dohfx) call cpsi_data_dealloc(Hpsi_r_tag)

      value = value.and.BA_pop_stack(vpsi(2))
      value = value.and.BA_pop_stack(vall(2))
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
      subroutine cpsi_Horb(hfxon,nb,ispin,ms,
     >                     psi_r_tag,
     >                     vall,
     >                     orb_k,orb_r,
     >                     Horb_k)
      implicit none
      logical    hfxon
      integer    nb,ispin,ms
      integer    psi_r_tag
      complex*16 vall(*)
      complex*16 orb_k(*)
      complex*16 orb_r(*)
      complex*16 Horb_k(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      logical dohfx
      integer nx,ny,nz,npack1
      real*8  scal1

*     **** tempory space ****
      integer nfft3d,ne1,shiftk
      integer vpsi(2),Horb_r(2)
      logical value


*     **** external functions ****
      logical  band_HFX,band_HFX_relaxed
      integer  cpsi_ne
      real*8   lattice_omega
      logical  control_spin_orbit
      complex*16 zdotc
      external band_HFX,band_HFX_relaxed
      external lattice_omega
      external control_spin_orbit
      external zdotc,cpsi_ne


      dohfx = band_HFX().and.band_HFX_relaxed().and.hfxon

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

      value = BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (dohfx) then
         value = value.and.
     >          BA_push_get(mt_dcpl,nfft3d,'Horb_r',Horb_r(2),Horb_r(1))
         call dcopy(2*nfft3d,0.0d0,0,dbl_mb(Horb_r(1)),1)
      end if
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


*     **********************
*     **** HFX energy   ****
*     **********************
      if (dohfx) then
         call band_potential_HFX_orb(nb,ms,
     >                               psi_r_tag,
     >                               orb_r,
     >                               dcpl_mb(Horb_r(1)))
      end if



*     **************************************************
*     **** Hpsi_k = Hpsi_k + FFT[vall_r+xcp)*orb_r] ****
*     **************************************************

*     **** vpsi_r = (vall_r)*orb_r ****
      call C3dB_bb_Mul(1,vall(1+(ms-1)*nfft3d),
     >                   orb_r,
     >                   dcpl_mb(vpsi(1)))

      if (dohfx)
     >   call C3dB_bb_Sum2(1,dcpl_mb(Horb_r(1)),
     >                      dcpl_mb(vpsi(1)))


*     **** put vpsi_r in k-space ****
      call C3dB_rc_pfft3f(1,nb,dcpl_mb(vpsi(1)))
      call Cram_c_pack(nb,dcpl_mb(vpsi(1)))

*     **** add to Horb_k ****
      call Cram_cc_daxpy(nb,(-scal1),
     >                   dcpl_mb(vpsi(1)),
     >                   Horb_k)

      value = .true.
      if (dohfx) value = value.and.BA_pop_stack(Horb_r(2))
      value = value.and.BA_pop_stack(vpsi(2))
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

#include "bafdecls.fh"
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

      value = BA_push_get(mt_dcpl,2*nfft3d,'vpsi',vpsi(2),vpsi(1))
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
      call C3dB_bb_Mul(1,vall(1),
     >                   orb_r,
     >                   dcpl_mb(vpsi(1)))
      call C3dB_bb_Mul(1,vall(1+nfft3d),
     >                   orb_r(shiftr),
     >                   dcpl_mb(vpsi(1)+nfft3d))


*     **** put vpsi_r in k-space ****
      call C3dB_rc_pfft3f(1,nb,dcpl_mb(vpsi(1)))
      call Cram_c_pack(nb,dcpl_mb(vpsi(1)))

      call C3dB_rc_pfft3f(1,nb,dcpl_mb(vpsi(1)+nfft3d))
      call Cram_c_pack(nb,dcpl_mb(vpsi(1)+nfft3d))

*     **** add to Horb_k ****
      call Cram_cc_daxpy(nb,(-scal1),
     >                   dcpl_mb(vpsi(1)),
     >                   Horb_k)
      call Cram_cc_daxpy(nb,(-scal1),
     >                   dcpl_mb(vpsi(1)+nfft3d),
     >                   Horb_k(shiftk))


      value = BA_pop_stack(vpsi(2))
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
#include "bafdecls.fh"

*     **** local variables ****
      integer ms,n,n1(2),n2(2)
      integer nx,ny,nz
      integer psi_shift,Hpsi_shift,nshift1,nshift2
      real*8  scal1,scal2,fion(3)

*     **** tempory space ****
      integer nfft3d,nb,nbq,nbrillq
      integer vpsi(2)
      logical value

*     **** external functions ****
      integer  Pneb_nbrillq,Pneb_convert_nb
      integer  cpsi_data_get_chnk,cpsi_data_nsize
      real*8   lattice_omega
      external Pneb_nbrillq,Pneb_convert_nb
      external cpsi_data_get_chnk,cpsi_data_nsize
      external lattice_omega
      

      call C3dB_nx(1,nx)
      call C3dB_ny(1,ny)
      call C3dB_nz(1,nz)
      call C3dB_nfft3d(1,nfft3d)
c      call Cram_max_npack(npack1)

      value = BA_push_get(mt_dcpl,nfft3d,'vpsi',vpsi(2),vpsi(1))
      if (.not.value) call errquit('cpsi_H_vall:out of stack',0,MA_ERR)
      

      scal2 = 1.0d0/lattice_omega()
      scal1 = 1.0d0/dble(nx*ny*nz)
      n1(1) = 1
      n2(1) = ne(1)
      n1(2) = ne(1) + 1
      n2(2) = ne(1) + ne(2)
      nbrillq = Pneb_nbrillq()
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
      do nbq=1,nbrillq
cccc         nb = Pneb_convert_nb(nbq)
         psi_shift  = cpsi_data_get_chnk( psi_r_tag,nbq)
         Hpsi_shift = cpsi_data_get_chnk(Hpsi_k_tag,nbq)
    
         do ms=1,ispin
         do n=n1(ms),n2(ms)

*           *************************************
*           **** vpsi_r = (vall_r)*psi_r ****
*           *************************************
            call C3dB_bb_Mul(1,vall(1+(ms-1)*nfft3d),
     >                         dbl_mb(psi_shift),
     >                         dcpl_mb(vpsi(1)))

*           *******************************
*           **** get vpsi_r in k-space ****
*           *******************************
            !call C3dB_rc_fft3f(1,dcpl_mb(vpsi(1)))
            call C3dB_rc_pfft3f(1,nbq,dcpl_mb(vpsi(1)))
            call Cram_c_pack(nbq,dcpl_mb(vpsi(1)))

*           ***********************
*           **** add to Hpsi_k ****
*           ***********************
            call Cram_cc_daxpy(nbq,(-scal1),
     >                         dcpl_mb(vpsi(1)),
     >                         dbl_mb(Hpsi_shift))

            psi_shift  = psi_shift  + nshift1
            Hpsi_shift = Hpsi_shift + nshift2
         end do
         end do
      end do

      value = BA_pop_stack(vpsi(2))
      if(.not.value) call errquit('cpsi_H_vall:popping stack',0,MA_ERR)

      return
      end

