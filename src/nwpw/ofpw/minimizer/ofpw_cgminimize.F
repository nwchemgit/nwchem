c
c     $Id$                       
c

*  ************************************************************
*  *                MPI ofpw_cgminimize routine               *
*  *                                                          *
*  ************************************************************
  
      subroutine ofpw_cgminimize(E,deltae,deltac,
     >                           current_iteration,it_in)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      integer    it_in

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer H0(2),G0(2),G1(2)
      real*8  E0,dE0


      logical    precondition
      common / cgsd_block2 / precondition

      real*8     sum0,sum1,sum3,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      integer it
      real*8 tmin,deltat
      real*8 max_sigma

      logical value
      integer neall,npack1
      real*8 e_lj,e_q,e_spring
 


*     **** external functions ****
      integer  psi_neq,control_version
      real*8   control_tole,control_tolc
      real*8   psi_geodesic_energy
      real*8   psi_geodesic_denergy
      real*8   rho_error
      real*8   dng_1ehartree
      real*8   psi_1ke
      real*8   psi_1vl,psi_1v_field
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
      real*8   linesearch
   
      external psi_neq,control_version
      external control_tole,control_tolc
      external psi_geodesic_energy
      external psi_geodesic_denergy
      external rho_error
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
      real*8   pspw_qmmm_LJ_E
      real*8   pspw_qmmm_Q_E
      real*8   pspw_qmmm_spring_E
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E
      external pspw_qmmm_Q_E
      external pspw_qmmm_spring_E

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

      logical  control_precondition,pspw_HFX
      external control_precondition,pspw_HFX
       
  


      deltae = 0.0d0
      deltac = 0.0d0



      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
      if (control_version().eq.4) eion = ion_ion_e()

      eorbit   = 0.0d0
      ehartree = 0.0d0
      exc      = 0.0d0
      pxc      = 0.0d0

      E(1)  = Enew + eion
      E(2)  = eorbit
      E(3)  = ehartree
      E(4)  = exc
      E(5)  = eion
      E(6)  = 0.0d0
      E(7)  = 0.0d0
      E(8)  = 0.0d0
      E(9)  = 2.0d0*ehartree
      E(10) = pxc

*     **** set pspw_qmmm energies ****
      if (pspw_qmmm_found()) then
         e_lj     = pspw_qmmm_LJ_E()
         e_q      = pspw_qmmm_Q_E()
         e_spring = pspw_qmmm_spring_E()
         E(1)  = E(1) + e_lj + e_q + e_spring

         E(11) = e_lj
         E(12) = e_q
         E(13) = e_spring
      end if

*     **** get pspw_charge  energies ****
      if (pspw_charge_found()) then
         !E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()
         E(21)  = pspw_charge_Energy_charge()
         E(1)   = E(1) + E(20) + E(21)
      end if

      return
      end
 

