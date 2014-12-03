c
c     $Id$                       
c

*  ************************************************************
*  *                                                          *
*  *               Limited Memory BFGS routine                *
*  *                                                          *
*  ************************************************************
  
      subroutine paw_bfgsminimize(E,deltae,deltac,
     >                        current_iteration,
     >                        minimizer)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      integer    minimizer

#include "bafdecls.fh"
#include "errquit.fh"
#include "paw_energy_kin_atom.fh"
#include "paw_energy_vloc_atom.fh"
#include "paw_energy_ion_atom.fh"
#include "paw_energy_core_atom.fh"
#include "paw_energy_hartree_atom.fh"
#include "paw_xc.fh"


*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-2)
       
      integer G0(2),S0(2)
      real*8  E0,dE0

      logical    precondition
      common / paw_cgsd_block2 / precondition

      real*8     tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / paw_cgsd_block / Enew,Eold,Estart

      real*8 tmin,dte,sum0,sum1
      common / paw_bfgs_block / tmin,dte,sum0,sum1

      integer i,it,it_in
      real*8 deltat
      real*8 max_sigma,dt,kappa

      logical value
      integer neall,npack1

      double precision evloc_pw,evloc_atom,ehartr_pw,eke
      double precision ekin_atom
      double precision eion_atom
      double precision ecore_atom
      double precision ecore_ion_atom
      double precision ecore_self_atom
      double precision ehartree_atom
      double precision exc_atom,pxc_atom



*     **** external functions ****
      integer  paw_psi_ne,control_version
      real*8   control_tole,control_tolc
      real*8   paw_psi_geodesic_energy
      real*8   paw_psi_geodesic_denergy
      real*8   paw_rho_error
      real*8   paw_dng_1ehartree
      real*8   paw_psi_1ke
      real*8   paw_psi_1vl
      real*8   paw_rho_1exc
      real*8   paw_rho_1pxc
      real*8   paw_psi_1eorbit
      real*8   linesearch
      real*8   paw_mult_energy_atom_comp !**no header file for paw_mult**
      real*8   paw_mult_energy_atom_self
      real*8   paw_mult_energy_atom_mult

      external paw_psi_ne,control_version
      external control_tole,control_tolc
      external paw_psi_geodesic_energy
      external paw_psi_geodesic_denergy
      external paw_rho_error
      external paw_dng_1ehartree
      external paw_psi_1ke
      external paw_psi_1vl
      external paw_rho_1exc
      external paw_rho_1pxc
      external paw_psi_1eorbit
      external linesearch
      external paw_mult_energy_atom_comp
      external paw_mult_energy_atom_self
      external paw_mult_energy_atom_mult


      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      integer  control_lmbfgs_size,control_it_in
      external control_lmbfgs_size,control_it_in

      dt = control_time_step()
c     dte = dt/dsqrt(control_fake_mass())
c     dte = dsqrt(control_fake_mass())/dt
c     dte = 1.0d0

  
      call Pack_npack(1,npack1)
      neall = paw_psi_ne(1)+paw_psi_ne(2)


*     **** allocate G0, S0 ****
      value = BA_push_get(mt_dcpl,npack1*neall,
     >                     'S0',S0(2),S0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G0',G0(2),G0(1))
      if (.not. value) call errquit('cgminimize:out of heap memory',0,
     &       MA_ERR)



*     ***** get the initial gradient ****
      call paw_psi_1get_residual(dcpl_mb(S0(1)),dcpl_mb(G0(1)),E0)
      call dcopy(2*neall*npack1,dcpl_mb(S0(1)),1,dcpl_mb(G0(1)),1)
      call paw_psi_1gen_Tangent(dcpl_mb(G0(1)))

      if (current_iteration.eq.1) then

*        ***** use the initial gradient for the direction ****
         call paw_lmbfgs_init(control_lmbfgs_size(),dcpl_mb(G0(1)))
         
         call Grsm_gg_Copy(npack1,neall,
     >                      dcpl_mb(G0(1)),
     >                      dcpl_mb(S0(1)))
         tmin  = 10*deltat_min
      else
         call paw_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
      end if

*     ***** set to gradient if direction is assending *****
      call Grsm_gg_trace(npack1,neall,
     >                   dcpl_mb(S0(1)),
     >                   dcpl_mb(S0(1)),
     >                   kappa)
      if (kappa.le.1.0d-10) call Grsm_gg_Copy(npack1,neall,
     >                                      dcpl_mb(G0(1)),
     >                                      dcpl_mb(S0(1)))


     

*     ******************************************
*     ****                                  ****
*     **** Start of BFGS iteration loop     ****
*     ****                                  ****
*     ******************************************
      it_in = control_it_in()
      tole  = control_tole()
      tolc  = control_tolc()
      do it=2,it_in

*        **** initialize the geoedesic line data structure ****
         call paw_geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)

*        ******* line search *********
         if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        paw_psi_geodesic_energy,
     >                        paw_psi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
         call paw_psi_geodesic_final(tmin)
         deltac = paw_rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            go to 30
         end if
     
*        **** get the new gradient - also updates densities****
         call paw_psi_2get_residual(2,dcpl_mb(S0(1)),dcpl_mb(G0(1)),E0)
         call dcopy(2*neall*npack1,dcpl_mb(S0(1)),1,dcpl_mb(G0(1)),1)
         call paw_psi_2gen_Tangent(dcpl_mb(G0(1)))


         call paw_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))


*        ***** set to gradient if direction is assending *****
         call Grsm_gg_trace(npack1,neall,
     >                      dcpl_mb(S0(1)),
     >                      dcpl_mb(S0(1)),kappa)
         if (kappa.le.1.0d-10) call Grsm_gg_Copy(npack1,neall,
     >                                         dcpl_mb(G0(1)),
     >                                         dcpl_mb(S0(1)))

*        **** make psi1 <--- psi2(tmin) ****
         call paw_psi_2to1()


      end do

*     **** initialize the geoedesic line data structure ****
      call paw_geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)

*     ******* line search *********
      if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

      Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        paw_psi_geodesic_energy,
     >                        paw_psi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
      call paw_psi_geodesic_final(tmin)
      deltac = paw_rho_error()

*     **** free memory ****
 30   value = BA_pop_stack(G0(2))
      value = value.and.
     >        BA_pop_stack(S0(2))
      if (.not. value) 
     >  call errquit('paw_bfgsminimize:error freeing heap memory',
     >               0,MA_ERR)

      call paw_psi_2to1()
      call paw_psi_1toelectron() 
      call paw_psi_check()

*     **** plane-wave energies ****
      eorbit    = paw_psi_1eorbit()
      eke       = paw_psi_1ke()
      ehartr_pw = paw_dng_1ehartree()
      exc       = paw_rho_1exc()
      pxc       = paw_rho_1pxc()
      evloc_pw  = paw_psi_1vl()


*     **** atomic energies ***
      ekin_atom     = paw_energy_kin_atom()
      evloc_atom    = paw_energy_vloc_atom()

      ehartree_atom   = paw_energy_hartree_atom()
      eion_atom       = paw_energy_ion_atom()
      ecore_atom      = paw_energy_core_atom()
      ecore_ion_atom  = paw_energy_core_ion_atom()
      ecore_self_atom = paw_energy_core_self_atom()

      exc_atom        = paw_energy_xc_atom()


*     *** kinetic energy
      E(2) = eke
      E(3) = ekin_atom

*     *** coulomb contributions
      E(4) = ehartr_pw
      E(5) = eion_atom + ecore_atom + ehartree_atom +
     >       ecore_ion_atom + ecore_self_atom +
     >       paw_mult_energy_atom_self() +
     >       paw_mult_energy_atom_comp()

      E(6)=paw_mult_energy_atom_mult()

*     *** exch-correlation
      E(7) = exc
      E(8) = exc_atom

*     *** local pseudopot ***
      E(9)  = evloc_pw
      E(10) = evloc_atom


*     *** total energy ***
      Eold=E(1)
      E(1) = 0.0d0
      do i=2,10
       E(1) = E(1) + E(i)
      end do


      return
      end
 

