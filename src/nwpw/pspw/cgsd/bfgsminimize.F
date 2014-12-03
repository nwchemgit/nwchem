c
c     $Id$                       
c

*  ************************************************************
*  *                                                          *
*  *               Limited Memory BFGS routine                *
*  *                                                          *
*  ************************************************************
  
      subroutine bfgsminimize(E,deltae,deltac,
     >                        current_iteration,
     >                        minimizer)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      integer    minimizer

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-2)
       
      integer G0(2),S0(2)
      real*8  E0,dE0

      logical    precondition
      common / cgsd_block2 / precondition

      real*8     tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      real*8 tmin,dte,sum0,sum1
      common / bfgs_block / tmin,dte,sum0,sum1

      integer it,it_in
      real*8 deltat
      real*8 max_sigma,dt,kappa

      logical value,pawexist
      integer neall,npack1
      !real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      !real*8 e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      real*8 e_lj,e_q,e_spring
      real*8 ehfx,phfx
      real*8 exc_atom,pxc_atom,ehartree_atom,ecmp_cmp,ecmp_pw


*     **** external functions ****
      integer  control_it_in,psi_neq,control_version
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
   
      external control_it_in,psi_neq,control_version
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
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

*     ***** PAW functions ****
      logical  psp_pawexist
      real*8   psp_kinetic_core,psi_1ke_atom
      real*8   psi_1ncmp_vloc,psi_1vloc_atom
      real*8   psp_ion_core,psi_1valence_core_atom,psi_1hartree_atom
      real*8   psi_1hartree_cmp_cmp,psi_1hartree_cmp_pw
      external psp_pawexist
      external psp_kinetic_core,psi_1ke_atom
      external psi_1ncmp_vloc,psi_1vloc_atom
      external psp_ion_core,psi_1valence_core_atom,psi_1hartree_atom
      external psi_1hartree_cmp_cmp,psi_1hartree_cmp_pw

      real*8   control_time_step,control_fake_mass,psi_1meta_gga_pxc
      external control_time_step,control_fake_mass,psi_1meta_gga_pxc
      integer  control_lmbfgs_size
      external control_lmbfgs_size
      logical  control_precondition,pspw_HFX,psp_U_psputerm,meta_found
      external control_precondition,pspw_HFX,psp_U_psputerm,meta_found
      logical  nwpw_meta_gga_on,ion_disp_on
      external nwpw_meta_gga_on,ion_disp_on
      real*8   ion_disp_energy
      external ion_disp_energy


      dt = control_time_step()
c     dte = dt/dsqrt(control_fake_mass())
c     dte = dsqrt(control_fake_mass())/dt
c     dte = 1.0d0

      if (current_iteration.eq.1) then
         precondition = control_precondition()
      end if

  
      call Pack_npack(1,npack1)
      neall = psi_neq(1)+psi_neq(2)
      pawexist = psp_pawexist()


*     **** allocate G0, S0 ****
      value = BA_push_get(mt_dcpl,npack1*neall,
     >                     'S0',S0(2),S0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G0',G0(2),G0(1))
      if (.not. value) call errquit('bfgsminimize:out of heap memory',0,
     &       MA_ERR)
      call dcopy(2*npack1*neall,0.0d0,0,dcpl_mb(G0(1)),1)



*     ***** get the initial direction ****
      if (pawexist) then
         call psi_1get_STgradient(dcpl_mb(S0(1)),dcpl_mb(G0(1)),E0)
      else
         if (precondition) then
            call psi_1get_TMgradient(dcpl_mb(G0(1)),E0)
         else
            call psi_1get_Tgradient(dcpl_mb(G0(1)),E0)
         end if
      end if

      if (current_iteration.eq.1) then

*        ***** use the initial gradient for the direction ****
         if (minimizer.eq.3) then
           call pspw_diis_init(control_lmbfgs_size(),dcpl_mb(G0(1)))
         else
           call pspw_lmbfgs_init(control_lmbfgs_size(),dcpl_mb(G0(1)))
         end if
         call Grsm_gg_Copy(npack1,neall,
     >                      dcpl_mb(G0(1)),
     >                      dcpl_mb(S0(1)))
         tmin  = 10*deltat_min
      else
         if (minimizer.eq.3) then
           call pspw_diis(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
         else
           call pspw_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
         end if
      end if

*     ***** set to gradient if direction is assending *****
c      call Grsm_gg_trace(npack1,neall,
c     >                   dcpl_mb(S0(1)),
c     >                   dcpl_mb(S0(1)),
c     >                   kappa)
c      call D1dB_SumAll(kappa)
c      if (kappa.le.1.0d-10) call Grsm_gg_Copy(npack1,neall,
c     >                                      dcpl_mb(G0(1)),
c     >                                      dcpl_mb(S0(1)))


     

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
         call geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)

*        ******* line search *********
         if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
c         deltat = 1.0d0  !** debug **
         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        psi_geodesic_energy,
     >                        psi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
         call psi_geodesic_final(tmin)
         deltac = rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            if (.not.precondition) go to 30
            precondition = .false.
         end if
     
*        **** get the new gradient - also updates densities****
         if (pawexist) then
              call psi_2get_STgradient(2,dcpl_mb(S0(1)),
     >                                   dcpl_mb(G0(1)),E0)
         else
            if (precondition) then
              call psi_2get_TMgradient(2,dcpl_mb(G0(1)),E0)
            else
              call psi_2get_Tgradient(2,dcpl_mb(G0(1)),E0)
            end if
         end if

         if (minimizer.eq.3) then
           call pspw_diis(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
         else
           call pspw_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
         end if

*        ***** set to gradient if direction is assending *****
c         call Grsm_gg_trace(npack1,neall,
c     >                      dcpl_mb(S0(1)),
c     >                      dcpl_mb(S0(1)),kappa)
c      !write(*,*) "KAPPA=",kappa
c         if (kappa.le.1.0d-10) call Grsm_gg_Copy(npack1,neall,
c     >                                         dcpl_mb(G0(1)),
c     >                                         dcpl_mb(S0(1)))



*        **** make psi1 <--- psi2(tmin) ****
         call psi_2to1()

         !****DEBUG ****  check to see if S0 is orthogonal to psi1
         !call psi_1check_Tangent(dcpl_mb(S0(1)))
         !****DEBUG ****  

c        call rho_2to1()
c        call dng_2to1()

      end do

*     **** initialize the geoedesic line data structure ****
      call geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)

*     ******* line search *********
      if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
         deltat = tmin
      else
         deltat = deltat_min
      end if
c      deltat = 1.0d0  !** debug **
c      Eold = psi_geodesic_energy(0.0d0)
      Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        psi_geodesic_energy,
     >                        psi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
      call psi_geodesic_final(tmin)
      deltac = rho_error()

*     **** free memory ****
 30   value = BA_pop_stack(G0(2))
      value = value.and.
     >        BA_pop_stack(S0(2))
      if (.not. value) 
     >  call errquit('bfgsminimize:error freeing heap memory',0, MA_ERR)

      call psi_2to1()
      call psi_1toelectron() 
c     call rho_2to1()
c     call dng_2to1()
c     call psi_check()

      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
      if (control_version().eq.4) eion = ion_ion_e()

      eorbit   = psi_1eorbit()
      ehartree = dng_1ehartree()
      exc      = rho_1exc()
      pxc      = rho_1pxc()

*     **** Metadynamics GGA Tau potential term ****
      if (nwpw_meta_gga_on()) then
         pxc = pxc + psi_1meta_gga_pxc()
      end if

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
         e_lj     = pspw_qmmm_LJ_E()
         e_q      = pspw_qmmm_Q_E()
         e_spring = pspw_qmmm_spring_E()
         E(1)  = E(1) + e_lj + e_q + e_spring

         E(11) = e_lj
         E(12) = e_q
         E(13) = e_spring
      end if

*     **** get pspw_charge energies ****
      if (pspw_charge_found()) then
         E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()
         E(21)  = pspw_charge_Energy_charge()
         E(1)   = E(1) + E(20) + E(21)
      end if

*     **** paw energies ****
      if (pawexist) then
         E(35) = psi_1ke_atom()            !*** kinetic atom ***
         E(36) = psp_kinetic_core()        !*** kinetic core  - independent of psi ***

         E(37) = psi_1valence_core_atom()  !*** valence core energy ***

         E(38) = psi_1ncmp_vloc() !*** (ncmp)*Vlocal
         E(39) = psi_1vloc_atom() !*** Vlocal (local basis)-nv*Vsrloc+nloc*Vloc-ncmp*Vloc ***
         ehartree_atom = psi_1hartree_atom()
         ecmp_cmp      = psi_1hartree_cmp_cmp()
         ecmp_pw       = psi_1hartree_cmp_pw()
         E(40) = ehartree_atom             !*** vcoulomb atom  ***
         E(41) = ecmp_cmp                  !*** ncmp-ncmp coulomb energy ***
         E(42) = ecmp_pw                   !*** ncmp-pw coulomb energy ***

         call psi_1xc_atom(exc_atom,pxc_atom)
         E(43) = exc_atom                  !*** exc atom  ***
         E(44) = pxc_atom                  !*** pxc atom  ***

         E(45) = psp_ion_core()            !*** ion core energy - independent of psi ***
      end if

*     **** HFX terms ****
      if (pspw_HFX()) then
         call electron_HFX_energies(ehfx,phfx)
         E(26) = ehfx
         E(27) = phfx
      end if

*     **** DFT+U terms ****
      if (psp_U_psputerm()) then
         call electron_U_energies(ehfx,phfx)
         E(29) =  ehfx
         E(30) =  phfx
      end if

*     **** Metadynamics potential terms ****
      if (meta_found()) then
         call electron_meta_energies(ehfx,phfx)
         E(31) =  ehfx
         E(32) =  phfx
      end if

*     **** Dispersion energy ****
      if (ion_disp_on()) then
         E(33) = ion_disp_energy()
         E(1)  = E(1) + E(33)
      end if

      return
      end
