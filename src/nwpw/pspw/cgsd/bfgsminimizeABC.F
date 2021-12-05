c
c     $Id$                       
c

*  ************************************************************
*  *                                                          *
*  *           Limited Memory BFGS ABC routine                *
*  *                                                          *
*  ************************************************************
      subroutine bfgsminimizeABC(E,deltae,deltac,
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
       
      integer G0(2),S0(2),G1(2),G2(2)
      real*8  E0,dE0,Enew0,Ermdr

      logical    lin_on,lin_on0
      common / cgsd_block9 / lin_on,lin_on0

      logical    precondition
      common / cgsd_block2 / precondition

      real*8     tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      real*8 tmin,dte,sum0,sum1
      common / bfgs_block / tmin,dte,sum0,sum1

      integer it,it_in
      real*8 deltat,tmin0,deltae0,deltac0
      real*8 max_sigma,dt,kappa

      logical value,pawexist
      integer neall,npack1
      !real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      !real*8 e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      real*8 e_lj,e_q,e_spring
      real*8 ehfx,phfx
      real*8 exc_atom,pxc_atom,ehartree_atom,ecmp_cmp,ecmp_pw

      logical cosmo_on,cosmo1_on,V_APC_on,field_exist
      real*8  eapc,papc

*     **** define pointers ****
      integer Y,U,HY,HU
      integer ispin,ne(2)


*     **** external functions ****
      integer  control_it_in,psi_neq,control_version,psi_ispin,psi_ne
      real*8   control_tole,control_tolc
      real*8   psi_geodesic_energy
      real*8   psi_geodesic_denergy
      real*8   rho_error
      real*8   dng_1ehartree
      real*8   psi_1ke
      real*8   psi_1vl,psi_1v_field,dng_1vl_mm
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
      real*8   linesearch
   
      external control_it_in,psi_neq,control_version,psi_ispin,psi_ne
      external control_tole,control_tolc
      external psi_geodesic_energy
      external psi_geodesic_denergy
      external rho_error
      external dng_1ehartree
      external psi_1ke
      external psi_1vl,psi_1v_field,dng_1vl_mm
      external psi_1vnl
      external rho_1exc
      external rho_1pxc
      external ewald_e,ion_ion_e
      external psi_1eorbit
      external linesearch

*     ***** QM/MM external functions ****
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      real*8   pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge
      logical  pspw_Efield_found
      external pspw_Efield_found
      integer  pspw_Efield_type
      external pspw_Efield_type
      real*8   pspw_Efield_Energy_ion
      external pspw_Efield_Energy_ion
      real*8   dipole_Efield_e,dipole_Efield_p
      external dipole_Efield_e,dipole_Efield_p

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
      logical  pspw_Lin_HFX_on,nwpw_cosmo_on,nwpw_cosmo1_on
      external pspw_Lin_HFX_on,nwpw_cosmo_on,nwpw_cosmo1_on
      logical  pspw_V_APC_on
      external pspw_V_APC_on
      real*8   psi_1energy,psi_1vl_cosmo
      external psi_1energy,psi_1vl_cosmo
      real*8   nwpw_cosmo_EQionq,nwpw_cosmo_Eqq
      external nwpw_cosmo_EQionq,nwpw_cosmo_Eqq

      real*8   geodesic_ABC_energy,geodesic_ABC_denergy
      external geodesic_ABC_energy,geodesic_ABC_denergy

      integer  psi_get_psi_ptr,geodesic_get_U_ptr
      external psi_get_psi_ptr,geodesic_get_U_ptr



      dt = control_time_step()
c     dte = dt/dsqrt(control_fake_mass())
c     dte = dsqrt(control_fake_mass())/dt
c     dte = 1.0d0

      if (current_iteration.eq.1) then
         precondition = control_precondition()
         lin_on = pspw_Lin_HFX_on()
         lin_on0= lin_on
         !if (lin_on) call linesearch_maxiter_set(4)
      end if

  
      call Pack_npack(1,npack1)
      ispin = psi_ispin()
      ne(1) = psi_ne(1)
      ne(2) = psi_ne(2)
      neall = psi_neq(1)+psi_neq(2)
      pawexist = psp_pawexist()
      cosmo_on  = nwpw_cosmo_on()
      cosmo1_on = nwpw_cosmo1_on()
      V_APC_on  = pspw_V_APC_on()
      field_exist = pspw_charge_found().or.pspw_Efield_found()

      call geodesic_ABC_init()

*     **** allocate G0, S0, G1, and G2 ****
      value = BA_push_get(mt_dcpl,npack1*neall,
     >                     'S0',S0(2),S0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G0',G0(2),G0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G2',G2(2),G2(1))
      if (.not.value) call errquit('bfgsminimizeABC:out of heap',0,
     &       MA_ERR)
      !call dcopy(2*npack1*neall,0.0d0,0,dcpl_mb(G0(1)),1)

*     **** set ptrs ****
      Y = psi_get_psi_ptr(1)
      U = geodesic_get_U_ptr()
      HY = G1(1)
      HU = G2(1)

!$OMP PARALLEL private(it,deltat,Enew0,tmin0,deltae0,deltac0)
      call Parallel_shared_vector_zero(.true.,2*npack1*neall,
     >                                 dcpl_mb(G0(1)))

*     ***** get the initial direction ****
      if (lin_on) call pspw_Lin_HFX_reset()
      if (pawexist) then
         call psi_1get_STgradient(dcpl_mb(S0(1)),dcpl_mb(G0(1)),Enew0)
      else
         if (precondition) then
            call psi_1get_TMgradient(dcpl_mb(G0(1)),Enew0)
         else
            call psi_1get_Tgradient(dcpl_mb(G0(1)),Enew0)
         end if
      end if
      call psi_1get_remainder0(Ermdr)
!$OMP MASTER
      E0 = Enew0
!$OMP END MASTER

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
!$OMP MASTER
         tmin  = 10*deltat_min
!$OMP END MASTER
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
!$OMP MASTER
      it_in = control_it_in()
      tole  = control_tole()
      tolc  = control_tolc()
!$OMP END MASTER
!$OMP BARRIER
      do it=2,it_in

*        **** initialize the geoedesic line data structure ****
         call geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)

*        **** Copy Hpsi_k to HY then ****
*        **** generate Hpsi using psi_k=U and copy Hpsi_k to HU ****
*        **** compute A,B,C, <S0|H|S0>, <U|H|S0> and <U|H|U> matrices ****
         call geodesic_ABC_start(ispin,ne,npack1,Ermdr,
     >                        dcpl_mb(Y), dcpl_mb(U),
     >                        dcpl_mb(HY),dcpl_mb(HU))

*        ******* line search *********
         if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
 20      continue
         tmin0 = tmin
         deltae0 = deltae
         call nwpw_cosmo_fix_q()
         Enew0 = linesearch(0.0d0,E0,dE0,deltat,
     >                        geodesic_ABC_energy,
     >                        geodesic_ABC_denergy,
     >                        0.50d0,tmin0,deltae0,2)
         call nwpw_cosmo_relax_q()
!$OMP MASTER
         tmin = tmin0
         deltae = deltae0
         Enew = Enew0
!$OMP END MASTER
!$OMP BARRIER

         !*** turn off Lin HFX ****
         if (lin_on .and.
     >       ((dabs(deltae).lt.0.1d0*tole).or.(deltae.gt.0.0d0))) then
!$OMP MASTER
            lin_on = .false.
!$OMP END MASTER
            call pspw_Lin_HFX_off(.true.)
            go to 20
         end if
         call pspw_Lin_HFX_off(.false.)

         call psi_geodesic_final(tmin)
         call psi_2toelectron()
         deltac0 = rho_error()
!$OMP MASTER
         lin_on = lin_on0
         deltac = deltac0
!$OMP END MASTER
!$OMP BARRIER


*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            if (.not.precondition) go to 30
            precondition = .false.
         end if

     
*        **** get the new gradient - also updates densities****
         if (lin_on) call pspw_Lin_HFX_reset()
         if (pawexist) then
              call psi_2get_STgradient(2,dcpl_mb(S0(1)),
     >                                   dcpl_mb(G0(1)),Enew0)
         else
            if (precondition) then
              call psi_2get_TMgradient(2,dcpl_mb(G0(1)),Enew0)
            else
              if (pspw_Lin_HFX_on()) then
                 call psi_2get_Tgradient(1,dcpl_mb(G0(1)),Enew0)
              else
                 call psi_2get_Tgradient(2,dcpl_mb(G0(1)),Enew0)
              end if
            end if
         end if
         call psi_2get_remainder0(Ermdr)
!$OMP MASTER
      E0 = Enew0
!$OMP END MASTER
!$OMP BARRIER


         if (minimizer.eq.3) then
           call pspw_diis(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
         else
           call pspw_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
         end if

*        **** make psi1 <--- psi2(tmin) ****
         call psi_2to1()

      end do

*     **** initialize the geoedesic line data structure ****
      call geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)

*     **** Copy Hpsi_k to HY then ****
*     **** generate Hpsi using psi_k=U and copy Hpsi_k to HU ****
*     **** compute A,B,C, <S0|H|S0>, <U|H|S0> and <U|H|U> matrices ****
      call geodesic_ABC_start(ispin,ne,npack1,Ermdr,
     >                        dcpl_mb(Y), dcpl_mb(U),
     >                        dcpl_mb(HY),dcpl_mb(HU))

*     ******* line search *********
      if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
         deltat = tmin
      else
         deltat = deltat_min
      end if
c      deltat = 1.0d0  !** debug **
c      Eold = psi_geodesic_energy(0.0d0)
 25   continue
      tmin0 = tmin
      deltae0 = deltae
      call nwpw_cosmo_fix_q()
      Enew0 = linesearch(0.0d0,E0,dE0,deltat,
     >                        geodesic_ABC_energy,
     >                        geodesic_ABC_denergy,
     >                        0.50d0,tmin0,deltae0,2)
      call nwpw_cosmo_relax_q()
!$OMP MASTER
         tmin = tmin0
         deltae = deltae0
         Enew = Enew0
!$OMP END MASTER
!$OMP BARRIER

      !*** turn off Lin HFX ****
      if (lin_on .and.
     >       ((dabs(deltae).lt.0.1d0*tole).or.(deltae.gt.0.0d0))) then
!$OMP MASTER
         lin_on = .false.
!$OMP END MASTER
         call pspw_Lin_HFX_off(.true.)
         go to 25
      end if
      call pspw_Lin_HFX_off(.false.)

      call psi_geodesic_final(tmin)
      call psi_2toelectron()
      deltac0 = rho_error()
!$OMP MASTER
      lin_on = lin_on0
      deltac = deltac0
!$OMP END MASTER
!$OMP BARRIER

 30   CONTINUE
         

*     **** free memory ****
      value = BA_pop_stack(G2(2))
      value = value.and.BA_pop_stack(G1(2))
      value = value.and.BA_pop_stack(G0(2))
      value = value.and.BA_pop_stack(S0(2))
      if (.not. value) 
     >  call errquit('bfgsminimizeABC:freeing heap',0,MA_ERR)

      call geodesic_ABC_finalize()

      call psi_2to1()
c      call psi_1toelectron() 
      call rho_2to1()
      call dng_2to1()
c     call psi_check()

      Enew0 = 0.0d0
      if (control_version().eq.3) Enew0 = ewald_e()
      if (control_version().eq.4) Enew0 = ion_ion_e()
!$OMP MASTER
      eion = Enew0
!$OMP END MASTER
     

      Enew0    = psi_1eorbit()
!$OMP MASTER
      eorbit = Enew0
!$OMP END MASTER

      Enew0  = dng_1ehartree()
!$OMP MASTER
      ehartree = Enew0
!$OMP END MASTER

      Enew0    = rho_1exc()
!$OMP MASTER
      exc = Enew0
!$OMP END MASTER

      Enew0    = rho_1pxc()
!$OMP MASTER
      pxc = Enew0
!$OMP END MASTER

      Enew0 = psi_1ke()
!$OMP MASTER
      E(6)  = Enew0
!$OMP END MASTER

      Enew0 = psi_1vl()
!$OMP MASTER
      E(7) = Enew0
!$OMP END MASTER

      Enew0 = psi_1vnl()
!$OMP MASTER
      E(8) = Enew0
!$OMP END MASTER
!$OMP END PARALLEL

*     **** Metadynamics GGA Tau potential term ****
      if (nwpw_meta_gga_on()) then
         pxc = pxc + psi_1meta_gga_pxc()
      end if

      E(1)  = Enew + eion
      E(2)  = eorbit
      E(3)  = ehartree
      E(4)  = exc
      E(5)  = eion
      !E(6)  = psi_1ke()
      !E(7)  = psi_1vl()
      !E(8)  = psi_1vnl()
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

*        **** Eqm-mm energy ***
         E(14) = pspw_qmmm_LJ_Emix()
         E(14) = E(14) + pspw_qmmm_Q_Emix()
         E(14) = E(14) + dng_1vl_mm()
      end if

*     **** get pspw_charge and pspw_Efield energies ****
      if (field_exist) then
         if (pspw_Efield_type().eq.0) then
            E(49) = dipole_Efield_e()
            E(50) = dipole_Efield_p()
            E(51) = 0.0d0
         else
            E(49)  = psi_1v_field()
            E(50)  = pspw_charge_Energy_ion()
     >             + pspw_Efield_Energy_ion()
            E(51)  = pspw_charge_Energy_charge()
            E(1)   = E(1) + E(50) + E(51)
         end if
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


*     **** COSMO terms ****
      if (cosmo_on) then

         !*** cosmo1 ****
         if (cosmo1_on) then
            E(46) = psi_1vl_cosmo()
            E(47) = nwpw_cosmo_EQionq()
            E(48) = nwpw_cosmo_Eqq()

         !*** cosmo2 ****
         else
            call electron_apc_energies(eapc,papc)
            E(22) = eapc
            E(23) = papc
            E(46) = eapc
            E(47) = nwpw_cosmo_EQionq() !** E(Qion|q)
            E(48) = nwpw_cosmo_Eqq()    !** E(q|q)

            !E(1)  = E(1) + E(22) - E(23) + E(47) + E(48)
            !E(1)  = E(1) + E(22) - E(23)
         end if

      else if (V_APC_on) then
         call electron_apc_energies(eapc,papc)
         E(22) = eapc
         E(23) = papc
         !E(1)  = E(1) + eapc - papc
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
