c
c     $Id: nominimize.F,v 1.1 2002-04-24 21:59:17 bylaska Exp $ 
c

      subroutine nominimize(E)
      implicit none
      real*8     E(*)

#include "mafdecls.fh"

*     **** local variables ****
      integer G1(2)
      real*8  E0

      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew


      logical value
      integer neall,npack1
      real*8 e_pol,e_vib,e_cav
      real*8  e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      

*     **** external functions ****
      integer  psi_ne,control_version
      real*8   dng_1ehartree
      real*8   psi_1ke
      real*8   psi_1vl,psi_1v_field
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
   
      external psi_ne,control_version
      external dng_1ehartree
      external psi_1ke
      external psi_1vl,psi_1v_field
      external psi_1vnl
      external rho_1exc
      external rho_1pxc
      external ewald_e,ion_ion_e
      external psi_1eorbit
      external linesearch

*     ***** QM/MM external functions ****
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_Energy_QMMM_Q
      real*8   pspw_qmmm_Energy_QMMM_LJ
      real*8   pspw_qmmm_Energy_MMMM_Q
      real*8   pspw_qmmm_Energy_MMMM_LJ
      real*8   pspw_qmmm_Energy_pol
      real*8   pspw_qmmm_Energy_vib
      real*8   pspw_qmmm_Energy_cav
      external pspw_qmmm_found
      external pspw_qmmm_Energy_QMMM_Q
      external pspw_qmmm_Energy_QMMM_LJ
      external pspw_qmmm_Energy_MMMM_Q
      external pspw_qmmm_Energy_MMMM_LJ
      external pspw_qmmm_Energy_pol
      external pspw_qmmm_Energy_vib
      external pspw_qmmm_Energy_cav

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

      call Pack_npack(1,npack1)
      neall = psi_ne(1)+psi_ne(2)


*     **** allocate G1 ****
      value = MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('nominimize:out of heap memory',0)


      call psi_1get_Tgradient(dcpl_mb(G1(1)),E0)
      Enew = E0

 30   call psi_2to1()
      call psi_1toelectron() 
c     call rho_2to1()
c     call dng_2to1()
      call psi_check()

      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
      if (control_version().eq.4) eion = ion_ion_e()

      eorbit   = psi_1eorbit()
      ehartree = dng_1ehartree()
      exc      = rho_1exc()
      pxc      = rho_1pxc()

      E(1)  = Enew + eion
      E(2)  = eorbit
      E(3)  = ehartree
      E(4)  = exc
      E(5)  = eion
      E(6)  = psi_1ke()
      E(7)  = psi_1vl()
      E(8)  = psi_1vnl()
      E(9)  = 2.0d0*ehartree
      E(10) = pxc
      if (pspw_qmmm_found()) then

*        **** get QM/MM  energies ****
         e_qmmm_e  = psi_1v_field()
         e_qmmm_q   = pspw_qmmm_Energy_QMMM_Q()
         e_qmmm_lj  = pspw_qmmm_Energy_QMMM_LJ()

         e_mmmm_q  = pspw_qmmm_Energy_MMMM_Q()
         e_mmmm_lj  = pspw_qmmm_Energy_MMMM_LJ()

         e_pol   = pspw_qmmm_Energy_pol()
         e_vib   = pspw_qmmm_Energy_vib()
         e_cav   = pspw_qmmm_Energy_cav()

         E(1) = E(1) + e_qmmm_q
     >               + e_qmmm_lj
     >               + e_mmmm_q
     >               + e_mmmm_lj
     >               + e_pol
     >               + e_vib
     >               + e_cav
         E(11) = e_qmmm_e
         E(12) = e_qmmm_q
         E(13) = e_qmmm_lj
         E(14) = e_mmmm_q
         E(15) = e_mmmm_lj
         E(16) = e_pol
         E(17) = e_vib
         E(18) = e_cav
      end if

*     **** get pspw_charge  energies ****
      if (pspw_charge_found()) then
         E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()
         E(21)  = pspw_charge_Energy_charge()
         E(1)   = E(1) + E(20) + E(21)
      end if

      value = MA_free_heap(G1(2))
      if (.not. value) 
     >  call errquit('nominimize:error freeing heap memory',0)

  
      return
      end
 

