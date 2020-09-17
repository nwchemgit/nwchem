c
c     $Id$                       
c

*  ************************************************************
*  *                                                          *
*  *               Limited Memory BFGS routine                *
*  *                                                          *
*  ************************************************************
  
      subroutine bfgsminimize3(E,deltae,deltac,
     >                        current_iteration,
     >                        minimizer)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      integer    minimizer

#include "bafdecls.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-4)
       
      integer G0(2),S0(2),Y0(2)
      real*8  E0,dE0


      real*8     tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      real*8 tmin,dte,sum0,sum1
      common / bfgs_block / tmin,dte,sum0,sum1

      integer it,it_in
      real*8 deltat
      real*8 max_sigma,dt

      logical value
      integer neall,npack1
      !real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      !real*8 e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      real*8 e_lj,e_q,e_spring
      real*8 ehfx,phfx

      logical cosmo_on,cosmo1_on,V_APC_on
      real*8  eapc,papc


*     **** external functions ****
      integer  control_it_in,psi_neq,control_version
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
   
      external control_it_in,psi_neq,control_version
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


      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      integer  control_lmbfgs_size
      external control_lmbfgs_size
   
      logical  pspw_HFX,psp_U_psputerm,meta_found
      external pspw_HFX,psp_U_psputerm,meta_found
      logical  nwpw_meta_gga_on,ion_disp_on
      external nwpw_meta_gga_on,ion_disp_on
      real*8   psi_1meta_gga_pxc,ion_disp_energy
      external psi_1meta_gga_pxc,ion_disp_energy

      logical  nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      external nwpw_cosmo_on,nwpw_cosmo1_on,pspw_V_APC_on
      real*8   psi_1vl_cosmo,nwpw_cosmo_EQionq,nwpw_cosmo_Eqq
      external psi_1vl_cosmo,nwpw_cosmo_EQionq,nwpw_cosmo_Eqq

  

      dt = control_time_step()
c     dte = dt/dsqrt(control_fake_mass())
c     dte = dsqrt(control_fake_mass())/dt
c     dte = 1.0d0

      
      call Pack_npack(1,npack1)
      neall = psi_neq(1)+psi_neq(2)
      cosmo_on  = nwpw_cosmo_on()
      cosmo1_on = nwpw_cosmo1_on()
      V_APC_on  = pspw_V_APC_on()


*     **** allocate G0, S0 ****
      value = BA_push_get(mt_dcpl,npack1*neall,
     >                     'S0',S0(2),S0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'G0',G0(2),G0(1))
      value = value.and.
     >        BA_push_get(mt_dcpl,npack1*neall,
     >                     'Y0',Y0(2),Y0(1))
      if (.not. value) 
     >   call errquit('bfgsminimize3:out of heap memory',0,0)
      call dcopy(2*npack1*neall,0.0d0,0,dcpl_mb(G0(1)),1)



*     ***** get the initial direction ****
      call psi_1get_psi(dcpl_mb(Y0(1)))
      call psi_1get_Gradient(dcpl_mb(G0(1)),E0)

      
      if (current_iteration.eq.1) then

*        ***** use the initial gradient for the direction ****
         call pspw_lmbfgs3_init(control_lmbfgs_size(),
     >                          dcpl_mb(Y0(1)),
     >                          dcpl_mb(G0(1)))
         
         call Grsm_gg_Copy(npack1,neall,
     >                      dcpl_mb(G0(1)),
     >                      dcpl_mb(S0(1)))
         call psi_1gen_Tangent(dcpl_mb(S0(1)))
         tmin  = 10*deltat_min
      else
         call pspw_lmbfgs3(tmin,dcpl_mb(Y0(1)),
     >                          dcpl_mb(G0(1)),
     >                          dcpl_mb(S0(1)))
         call psi_1gen_Tangent(dcpl_mb(S0(1)))
      end if

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
         !deltat = 1.0d0  !** debug **
         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        psi_geodesic_energy,
     >                        psi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
         call psi_geodesic_final(tmin)
         deltac = rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            go to 30
         end if
     
*        **** get the new gradient - also updates densities****
         call psi_2get_psi(dcpl_mb(Y0(1)))
         call psi_2get_Gradient(2,dcpl_mb(G0(1)),E0)
         call pspw_lmbfgs3(tmin,dcpl_mb(Y0(1)),
     >                          dcpl_mb(G0(1)),
     >                          dcpl_mb(S0(1)))
         call psi_2gen_Tangent(dcpl_mb(S0(1)))


*        **** make psi1 <--- psi2(tmin) ****
         call psi_2to1()

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
 30   value = BA_pop_stack(Y0(2))
      value = value.and.
     >        BA_pop_stack(G0(2))
      value = value.and.
     >        BA_pop_stack(S0(2))
      if (.not. value) 
     >  call errquit('bfgsminimize3:error freeing heap memory',0,0)

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

*        **** Eqm-mm energy ***
         E(14) = pspw_qmmm_LJ_Emix()
         E(14) = E(14) + pspw_qmmm_Q_Emix()
         E(14) = E(14) + dng_1vl_mm()

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


*     **** get pspw_charge energies ****
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

*     **** Metadynamics GGA Tau potential term ****
      if (nwpw_meta_gga_on()) then
         E(10) = E(10) - psi_1meta_gga_pxc()
      end if

*     **** Dispersion energy ****
      if (ion_disp_on()) then
         E(33) = ion_disp_energy()
         E(1)  = E(1) + E(33)
      end if



      return
      end
 

