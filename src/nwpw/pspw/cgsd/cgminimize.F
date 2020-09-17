c
c     $Id$                       
c

*  ************************************************************
*  *                MPI cgminimize routine                    *
*  *               (Fletcher-Reeves' steps)                   *
*  *                                                          *
*  *  This is a developing cgsdv3 parallel code wrtten for    *
*  *  NWChem                                                  *
*  *                                                          *
*  ************************************************************
  
      subroutine cgminimize(E,deltae,deltac,current_iteration,it_in)
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
       
c      complex*16 H0(nfft3d,nemax)
c      complex*16 G1(nfft3d,nemax)
      integer H0(2),R1(2),G1(2)
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

      logical value,pawexist
      integer neall,npack1
      real*8 e_lj,e_q,e_spring
      real*8 ehfx,phfx
      real*8 exc_atom,pxc_atom,ehartree_atom,ecmp_cmp,ecmp_pw
      real*8 Enew0,tmin0,deltae0,deltac0

      logical cosmo_on,cosmo1_on,V_APC_on
      real*8  eapc,papc
      

*     **** external functions ****
      integer  psi_neq,control_version
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
   
      external psi_neq,control_version
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
      real*8   pspw_qmmm_LJ_E
      real*8   pspw_qmmm_Q_E
      real*8   pspw_qmmm_spring_E
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E
      external pspw_qmmm_Q_E
      external pspw_qmmm_spring_E
      real*8   pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix
      external pspw_qmmm_LJ_Emix,pspw_qmmm_Q_Emix

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

      logical  control_precondition,pspw_HFX,psp_U_psputerm,meta_found
      external control_precondition,pspw_HFX,psp_U_psputerm,meta_found
      logical  nwpw_meta_gga_on,ion_disp_on
      external nwpw_meta_gga_on,ion_disp_on
      real*8   psi_1meta_gga_pxc,ion_disp_energy
      external psi_1meta_gga_pxc,ion_disp_energy
      real*8   psi_1energy
      external psi_1energy

      logical  nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      external nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      real*8   psi_1vl_cosmo,nwpw_cosmo_EQionq,nwpw_cosmo_Eqq
      external psi_1vl_cosmo,nwpw_cosmo_EQionq,nwpw_cosmo_Eqq

  
      call Pack_npack(1,npack1)
      neall = psi_neq(1)+psi_neq(2)
      pawexist = psp_pawexist()
      cosmo_on  = nwpw_cosmo_on()
      cosmo1_on = nwpw_cosmo1_on()
      V_APC_on  = pspw_V_APC_on()

      Estart = Enew
      if (current_iteration.eq.1) then
         precondition = control_precondition()
      end if

*     **** allocate H0, G1, and R1 ****
      value = BA_alloc_get(mt_dcpl,npack1*neall,
     >                     'H0',H0(2),H0(1))
c      value = value.and.
c     >        BA_alloc_get(mt_dcpl,npack1*neall,
c     >                     'G0',G0(2),G0(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      if (pawexist) then
         value = value.and.
     >           BA_alloc_get(mt_dcpl,npack1*neall,
     >                        'R1',R1(2),R1(1))
      end if
      if (.not. value) call errquit('cgminimize:out of heap memory',0,
     &       MA_ERR)
      !call dcopy(2*npack1*neall,0.0d0,0,dcpl_mb(G1(1)),1)

!$OMP PARALLEL private(it,deltat,Enew0,tmin0,deltae0,deltac0)
      call Parallel_shared_vector_zero(.true.,2*npack1*neall,
     >                                 dcpl_mb(G1(1)))

*     ***** get the initial gradient and direction ****
      call pspw_Lin_HFX_reset()
      if (pawexist) then
         call psi_1get_STgradient(dcpl_mb(R1(1)),dcpl_mb(G1(1)),Enew0)
         call Grsm_gg_trace(npack1,neall,
     >                      dcpl_mb(G1(1)),
     >                      dcpl_mb(R1(1)),sum1)
      else
         if (precondition) then
            call psi_1get_TMgradient(dcpl_mb(G1(1)),Enew0)
         else
            call psi_1get_Tgradient(dcpl_mb(G1(1)),Enew0)
         end if
         call Grsm_gg_trace(npack1,neall,
     >                      dcpl_mb(G1(1)),
     >                      dcpl_mb(G1(1)),sum1)
      end if
!$OMP MASTER
      E0 = Enew0
!$OMP END MASTER
      call D1dB_SumAll(sum1)

      call Grsm_gg_Copy(npack1,neall,
     >                   dcpl_mb(G1(1)),
     >                   dcpl_mb(H0(1)))


*     ******************************************
*     ****                                  ****
*     **** Start of conjugate gradient loop ****
*     ****                                  ****
*     ******************************************
      !it_in = control_it_in()
!$OMP MASTER
      tole  = control_tole()
      tolc  = control_tolc()
      tmin  = deltat_min
!$OMP END MASTER
      do it=2,it_in

*        **** initialize the geoedesic line data structure ****
         call geodesic_start(dcpl_mb(H0(1)),max_sigma,dE0)

*        ******* line search *********
         if (tmin.gt.deltat_min) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
c        Eold = psi_geodesic_energy(0.0d0)
         tmin0   = tmin
         deltae0 = deltae
         Enew0 = linesearch(0.0d0,E0,dE0,deltat,
     >                        psi_geodesic_energy,
     >                        psi_geodesic_denergy,
     >                        tole,tmin0,deltae0,1)
c        Enew = linesearch(0.0d0,E0,dE0,deltat,
c    >                        psi_geodesic_energy,
c    >                        psi_geodesic_denergy,
c    >                        0.5d0,tmin,deltae,2)
!$OMP MASTER
         tmin = tmin0
         deltae = deltae0
         Enew = Enew0
!$OMP END MASTER
!$OMP BARRIER
         call psi_geodesic_final(tmin)
         deltac0 = rho_error()
!$OMP MASTER
         deltac = deltac0
!$OMP END MASTER
!$OMP BARRIER


*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            if (.not.precondition)  go to 30
            precondition = .false.
         end if
   

*        **** transport the previous search directions ****
c         call Grsm_gg_Copy(npack1,neall,
c     >                     dcpl_mb(G1(1)),
c     >                     dcpl_mb(G0(1)))

         call psi_1geodesic_transport(tmin,dcpl_mb(H0(1)))
c         call psi_1geodesic_Gtransport(tmin,dcpl_mb(G0(1)))

*        **** make psi1 <--- psi2(tmin) ****
         call psi_2to1()
c        call rho_2to1()
c        call dng_2to1()

*        **** get the new gradient - also updates densities****
         call pspw_Lin_HFX_reset()
         sum0  = sum1
         if (pawexist) then
            call psi_1get_STgradient(dcpl_mb(R1(1)),dcpl_mb(G1(1)),
     >                               Enew0)
            call Grsm_gg_trace(npack1,neall,
     >                         dcpl_mb(G1(1)),
     >                         dcpl_mb(R1(1)),sum1)
         else
            if (precondition) then
              call psi_1get_TMgradient(dcpl_mb(G1(1)),Enew0)
            else
              call psi_1get_Tgradient(dcpl_mb(G1(1)),Enew0)
            end if
            call Grsm_gg_trace(npack1,neall,
     >                           dcpl_mb(G1(1)),
     >                           dcpl_mb(G1(1)),
     >                           sum1)
         end if
!$OMP MASTER
        E0 = Enew0
!$OMP END MASTER
         call D1dB_SumAll(sum1)

c         call Grsm_ggg_Sub(npack1,neall,
c     >                       dcpl_mb(G1(1)),
c     >                       dcpl_mb(G0(1)),
c     >                       dcpl_mb(G0(1)))
c         call Grsm_gg_trace(npack1,neall,
c     >                        dcpl_mb(G0(1)),
c     >                        dcpl_mb(G1(1)),
c     >                        sum3)
c         call D1dB_SumAll(sum3)

c        write(*,*) "sum1:",sum1,sum0,sum3,max_sigma
c        write(*,*) "Enew:",Enew

*        **** the new direction using Fletcher-Reeves ****
         if ( (dabs(deltae).le.(1.0d-2)).and.
     >        (tmin.gt.deltat_min)) then

!$OMP MASTER
           if (sum0.gt.1.0d-9) then
c              scale = sum3/sum0
              scale = sum1/sum0
           else
              scale = 0.0d0
           end if
!$OMP END MASTER
!$OMP BARRIER

c           call Grsm_gg_dScale(npack1,neall,scale,
c     >                         dcpl_mb(H0(1)),
c     >                         dcpl_mb(H0(1)))
c           call Grsm_ggg_Sum(npack1,neall,
c     >                       dcpl_mb(G1(1)),
c     >                       dcpl_mb(H0(1)),
c     >                       dcpl_mb(H0(1)))
           call Grsm_gg_dScale1(npack1,neall,scale,
     >                         dcpl_mb(H0(1)))
           call Grsm_ggg_Sum2(npack1,neall,
     >                       dcpl_mb(G1(1)),
     >                       dcpl_mb(H0(1)))


*          ***** set to gradient if direction is assending *****
*          call Grsm_gg_trace(neall,H0,G1,kappa)
*          if (kappa.le.0.0d0) call Grsm_gg_Copy(neall,G1,H0)

*        **** the new direction using steepest-descent ****
         else
              call Grsm_gg_Copy(npack1,neall,
     >                          dcpl_mb(G1(1)),
     >                          dcpl_mb(H0(1)))
         end if


      end do

*     **** initialize the geoedesic line data structure ****
      call geodesic_start(dcpl_mb(H0(1)),max_sigma,dE0)


*     ******* line search *********
      if (tmin.gt.deltat_min) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

c      Eold = psi_geodesic_energy(0.0d0)
      tmin0   = tmin
      deltae0 = deltae
      Enew0 = linesearch(0.0d0,E0,dE0,deltat,
     >                        psi_geodesic_energy,
     >                        psi_geodesic_denergy,
     >                        tole,tmin0,deltae0,1)
!$OMP MASTER
         tmin = tmin0
         deltae = deltae0
         Enew = Enew0
!$OMP END MASTER
!$OMP BARRIER
c     Enew = linesearch(0.0d0,E0,dE0,deltat,
c    >                        psi_geodesic_energy,
c    >                        psi_geodesic_denergy,
c    >                        0.5d0,tmin,deltae,2)
      call psi_geodesic_final(tmin)
      deltac0 = rho_error()
!$OMP MASTER
         deltac = deltac0
!$OMP END MASTER

 30   call psi_2to1()
      call psi_1toelectron() 
      !Enew = psi_1energy() 
c     call rho_2to1()
c     call dng_2to1()
c     call psi_check()
      !write(*,*) "Elast=",Enew

      Enew0 = 0.0d0
      if (control_version().eq.3) Enew0 = ewald_e()
      if (control_version().eq.4) Enew0 = ion_ion_e()
!OMP MASTER
      eion = Enew0
!OMP END MASTER

      Enew0    = psi_1eorbit()
!OMP MASTER
      eorbit   = Enew0
!OMP END MASTER

      Enew0    = dng_1ehartree()
!OMP MASTER
      ehartree = Enew0
!OMP END MASTER

      Enew0    = rho_1exc()
!OMP MASTER
      exc      = Enew0
!OMP END MASTER

      Enew0    = rho_1pxc()
!OMP MASTER
      pxc      = Enew0
!OMP END MASTER

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

*     **** set pspw_qmmm energies ****
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
            E(1)  = E(1) + E(22) - E(23)
         end if

      else if (V_APC_on) then
         call electron_apc_energies(eapc,papc)
         E(22) = eapc
         E(23) = papc
         E(1)  = E(1) + eapc - papc
      end if

*     **** get pspw_charge  energies ****
      if (pspw_charge_found()) then
         E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()
         E(21)  = pspw_charge_Energy_charge()
         E(1)   = E(1) + E(20) + E(21)
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


      value = BA_free_heap(G1(2))
      value = value.and.BA_free_heap(H0(2))
      if (pawexist) 
     >   value = value.and.BA_free_heap(R1(2))
      if (.not. value) 
     >  call errquit('cgminimize:error freeing heap memory',0, MA_ERR)

  
      return
      end
 

