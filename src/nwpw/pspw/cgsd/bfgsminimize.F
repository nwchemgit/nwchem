c
c     $Id: bfgsminimize.F,v 1.2 2001-12-21 22:54:35 bylaska Exp $                       
c

*  ************************************************************
*  *                                                          *
*  *               Limited Memory BFGS routine                *
*  *                                                          *
*  ************************************************************
  
      subroutine bfgsminimize(E,deltae,deltac,current_iteration)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration

#include "mafdecls.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer G0(2),S0(2)

      real*8     tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      real*8 tmin
      common / bfgs_block / tmin

      integer it,it_in
      real*8 deltat
      real*8 max_sigma,dt,dte

      logical value
      integer neall,npack1
      real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav




*     **** external functions ****
      integer  control_it_in,psi_ne,control_version
      real*8   control_tole,control_tolc
      real*8   psi_geodesic_energy
      real*8   psi_geodesic_denergy
      real*8   rho_error
      real*8   dng_1ehartree
      real*8   psi_1ke
      real*8   psi_1vl
      real*8   psi_1vnl
      real*8   rho_1exc
      real*8   rho_1pxc
      real*8   ewald_e,ion_ion_e
      real*8   psi_1eorbit
      real*8   linesearch
   
      external control_it_in,psi_ne,control_version
      external control_tole,control_tolc
      external psi_geodesic_energy
      external psi_geodesic_denergy
      external rho_error
      external dng_1ehartree
      external psi_1ke
      external psi_1vl
      external psi_1vnl
      external rho_1exc
      external rho_1pxc
      external ewald_e,ion_ion_e
      external psi_1eorbit
      external linesearch

*     ***** QM/MM external functions ****
      logical  pspw_qmmm_found
      real*8   pspw_qmmm_Energy_ionMM
      real*8   pspw_qmmm_Energy_QMMM
      real*8   pspw_qmmm_Energy_MMMM
      real*8   pspw_qmmm_Energy_pol
      real*8   pspw_qmmm_Energy_vib
      real*8   pspw_qmmm_Energy_cav
      external pspw_qmmm_found
      external pspw_qmmm_Energy_ionMM
      external pspw_qmmm_Energy_QMMM
      external pspw_qmmm_Energy_MMMM
      external pspw_qmmm_Energy_pol
      external pspw_qmmm_Energy_vib
      external pspw_qmmm_Energy_cav

      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass

      dt = control_time_step()
c     dte = dt/dsqrt(control_fake_mass())
      dte = dsqrt(control_fake_mass())/dt


      write(*,*) "start,tmin:",tmin
       
  
      call Pack_npack(1,npack1)
      neall = psi_ne(1)+psi_ne(2)


*     **** allocate G0, S0 ****
      value = MA_push_get(mt_dcpl,npack1*neall,
     >                     'S0',S0(2),S0(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1*neall,
     >                     'G0',G0(2),G0(1))
      if (.not. value) call errquit('cgminimize:out of heap memory',0)



*        ***** get the initial direction ****
      if (current_iteration.eq.1) then

*        ***** use the initial gradient for the direction ****
         call psi_1get_Tgradient(dcpl_mb(G0(1)))
         call Grsm_gg_dScale(npack1,neall,dte,
     >                       dcpl_mb(G0(1)),dcpl_mb(G0(1)))
         call Grsm_gg_Copy(npack1,neall,
     >                      dcpl_mb(G0(1)),
     >                      dcpl_mb(S0(1)))
         call pspw_lmbfgs_init(5,dcpl_mb(G0(1)))
      else
         call pspw_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
      end if

*     ******************************************
*     ****                                  ****
*     **** Start of BFGS iteration loop     ****
*     ****                                  ****
*     ******************************************
      it_in = control_it_in()
      tole  = control_tole()
      tolc  = control_tolc()
      tmin  = deltat_min
      do it=2,it_in

*        **** initialize the geoedesic line data structure ****
         call geodesic_start(dcpl_mb(S0(1)),max_sigma)

*        ******* line search *********
         if (tmin.gt.0.0d0) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
         Enew = linesearch(0.0d0,deltat,
     >                        psi_geodesic_energy,
     >                        psi_geodesic_denergy,
     >                        tole,tmin,deltae)
         call psi_geodesic_final(tmin)
      write(*,*) "tmin:",tmin
         deltac = rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) 
     >      go to 30
     
*        **** get the new gradient - also updates densities****
         call psi_2get_Tgradient(dcpl_mb(G0(1)))
         call Grsm_gg_dScale(npack1,neall,dte,
     >                       dcpl_mb(G0(1)),dcpl_mb(G0(1)))
         call pspw_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))


*        **** make psi1 <--- psi2(tmin) ****
         call psi_2to1()
c        call rho_2to1()
c        call dng_2to1()

      end do

*     **** initialize the geoedesic line data structure ****
      call geodesic_start(dcpl_mb(S0(1)),max_sigma)

*     ******* line search *********
      if (tmin.gt.0.0d0) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

c      Eold = psi_geodesic_energy(0.0d0)
      Enew = linesearch(0.0d0,deltat,
     >                        psi_geodesic_energy,
     >                        psi_geodesic_denergy,
     >                        tole,tmin,deltae)
      call psi_geodesic_final(tmin)
      deltac = rho_error()
      write(*,*) "END,tmin:",tmin

*     **** free memory ****
 30   value = MA_pop_stack(G0(2))
      value = value.and.
     >        MA_pop_stack(S0(2))
      if (.not. value) 
     >  call errquit('bfgsminimize:error freeing heap memory',0)

      call psi_2to1()
      call psi_1toelectron() 
c     call rho_2to1()
c     call dng_2to1()
      call psi_check()

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
         e_ionmm = pspw_qmmm_Energy_ionMM()
         e_qmmm  = pspw_qmmm_Energy_QMMM()
         e_mmmm  = pspw_qmmm_Energy_MMMM()
         e_pol   = pspw_qmmm_Energy_pol()
         e_vib   = pspw_qmmm_Energy_vib()
         e_cav   = pspw_qmmm_Energy_cav()

         E(1) = E(1) + e_ionmm
     >               + e_qmmm
     >               + e_mmmm
     >               + e_pol
     >               + e_vib
     >               + e_cav
         E(11) = e_qmmm
         E(12) = e_mmmm
         E(13) = e_pol
         E(14) = e_vib
         E(15) = e_cav
      end if


  
      return
      end
 

