c
c     $Id: bybminimize.F,v 1.2 2002-12-10 01:37:22 bylaska Exp $                       
c

*  ************************************************************
*  *                                                          *
*  *             Band by Band Kohn-Sham Minimizer             *
*  *                                                          *
*  *                                                          *
*  *                                                          *
*  ************************************************************
  
      subroutine bybminimize(E,deltae,deltac,current_iteration)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration

#include "mafdecls.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer H0(2),G0(2),G1(2)
      real*8  E0,dE0


      real*8     sum0,sum1,sum3,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      integer it,it_in,i
      real*8 tmin,deltat
      real*8 max_sigma

      logical value,precondition
      integer neall,npack1
      real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      real*8  e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      


*     **** external functions ****
      integer  control_it_in,psi_ne,control_version
      real*8   control_tole,control_tolc
      real*8   rho_error,psi_1energy
      real*8   dng_1ehartree
      real*8   psi_1ke
      real*8   psi_1vl,psi_1v_field
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
      real*8   linesearch
   
      external control_it_in,psi_ne,control_version
      external control_tole,control_tolc
      external rho_error,psi_1energy
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
      logical  pspw_charge_found,control_precondition
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found,control_precondition
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge


      it_in = control_it_in()
      Eold = E(1)
      precondition = control_precondition()

*     **********************
*     **** bybminimizer ****
*     **********************
      neall = psi_ne(1)+psi_ne(2)

      do it=1,it_in
         call psi_1toelectron()
         call electron_gen_vall() 
         do i=neall,1,-1
           call psi_1KS_update_orb(precondition,5,
     >                             control_tole(),
     >                             0.3d0,i)
         end do
      end do

      Enew = psi_1energy() 
      call psi_check()

      deltac=0.0d0

*     *****************************
*     **** set return energies **** - This is duplicate code
*     *****************************
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

      deltae = E(1)-Eold
      return
      end
 

