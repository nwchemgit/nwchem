c
c     $Id$ 
c
*     ******************************************
*     *                                        *
*     *             nominimize                 *
*     *                                        *
*     ******************************************
      subroutine nominimize(E)
      implicit none
      real*8     E(*)

#include "bafdecls.fh"
#include "errquit.fh"

*     **** local variables ****
      integer G1(2)
      real*8  E0

      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew


      logical value
      integer neall,npack1
      !real*8 e_pol,e_vib,e_cav
      !real*8  e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      real*8  e_lj,e_q,e_spring
      real*8  ehsic,phsic,exsic,pxsic,ehfx,phfx
      real*8  exc_atom,pxc_atom,ehartree_atom,ecmp_cmp,ecmp_pw

      logical cosmo_on,cosmo1_on,V_APC_on
      real*8  eapc,papc
      

*     **** external functions ****
      integer  psi_neq,control_version
      real*8   dng_1ehartree
      real*8   psi_1ke
      real*8   psi_1vl,psi_1v_field,dng_1vl_mm
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
   
      external psi_neq,control_version
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

*     ***** pspw_SIC and pspw_HFX external functions ****
      logical  pspw_SIC,pspw_HFX,psp_U_psputerm,meta_found
      external pspw_SIC,pspw_HFX,psp_U_psputerm,meta_found

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

      logical  nwpw_meta_gga_on,ion_disp_on
      external nwpw_meta_gga_on,ion_disp_on
      real*8   psi_1meta_gga_pxc,ion_disp_energy
      external psi_1meta_gga_pxc,ion_disp_energy

      logical  nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      external nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      real*8   psi_1vl_cosmo,nwpw_cosmo_EQionq,nwpw_cosmo_Eqq
      external psi_1vl_cosmo,nwpw_cosmo_EQionq,nwpw_cosmo_Eqq


      call Pack_npack(1,npack1)
      neall = psi_neq(1)+psi_neq(2)
      cosmo_on  = nwpw_cosmo_on()
      cosmo1_on = nwpw_cosmo1_on()
      V_APC_on  = pspw_V_APC_on()

*     **** allocate G1 ****
      value = BA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('nominimize:out of heap memory',0,
     &       MA_ERR)


      call psi_1get_Tgradient(dcpl_mb(G1(1)),E0)
      Enew = E0

c      call psi_1toelectron() 

      eion = 0.0d0
      if (control_version().eq.3) eion = ewald_e()
      if (control_version().eq.4) eion = ion_ion_e()

      eorbit   = psi_1eorbit()
      ehartree = dng_1ehartree()
      exc      = rho_1exc()
      pxc      = rho_1pxc()


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

*     ***** qmmm energies *****
      if (pspw_qmmm_found()) then
         e_lj     = pspw_qmmm_LJ_E()     !*** independent of psi ***
         e_q      = pspw_qmmm_Q_E()      !*** independent of psi ***
         e_spring = pspw_qmmm_spring_E() !*** independent of psi ***
         E(1)  = E(1) + e_lj + e_q + e_spring

         E(11) = e_lj
         E(12) = e_q
         E(13) = e_spring

*        **** Eqm-mm energy ***
         E(14) = pspw_qmmm_LJ_Emix()
         E(14) = E(14) + pspw_qmmm_Q_Emix()
         E(14) = E(14) + dng_1vl_mm()

      end if

*     **** paw energies ****
      if (psp_pawexist()) then
         E(35) = psi_1ke_atom()            !*** kinetic atom ***
         E(36) = psp_kinetic_core()        !*** kinetic core  - independent of psi ***

         E(37) = psi_1valence_core_atom()  !*** valence core energy ***


         E(38) = psi_1ncmp_vloc()          !*** (ncmp)*Vlocal
         E(39) = psi_1vloc_atom()          !*** Vlocal (local basis) -  nv*Vsrloc + nloc*Vloc - ncmp*Vloc ***

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

         write(*,*) "E1          =",E(1)
c         write(*,*) "total energy=",E(3)+E(4)+E(5)+E(6)+E(7) 
c     >                            + E(35)+E(36)+E(37)+E(38)+E(39)
c     >                            + E(40)+E(41)+E(42)+E(43)+E(45)
         write(*,*) "total energy=",E(3)+E(4)+E(5)+E(6)+E(7) 
     >                            + E(35)+E(37)+E(38)+E(39)
     >                            + E(40)+E(41)+E(42)+E(43)
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



*     **** get pspw_charge  energies ****
      if (pspw_charge_found()) then
         E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()    !*** independent of psi ***
         E(21)  = pspw_charge_Energy_charge() !*** independent of psi ***
         E(1)   = E(1) + E(20) + E(21)
      end if


*     **** SIC corrections ****
      if (pspw_SIC()) then
         call electron_SIC_energies(ehsic,phsic,exsic,pxsic)
         E(22) = ehsic
         E(23) = exsic
         E(24) = phsic
         E(25) = pxsic
      end if

*     **** HFX terms ****
      if (pspw_HFX()) then
         call electron_HFX_energies(ehfx,phfx)
         E(26) = ehfx
         E(27) = phfx
      end if

*     **** HFX terms ****
      if (psp_U_psputerm()) then
         call electron_U_energies(ehfx,phfx)
         E(29) = ehfx
         E(30) = phfx
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



      value = BA_free_heap(G1(2))
      if (.not. value) 
     >  call errquit('nominimize:error freeing heap memory',0, MA_ERR)

  
      return
      end
 

