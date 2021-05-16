c
c     $Id$                       
c

*  ************************************************************
*  *                                                          *
*  *               Limited Memory BFGS routine                *
*  *                                                          *
*  ************************************************************
  
      subroutine c_bfgsminimize(E,deltae,deltac,
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
       
      integer G0_tag,S0_tag,G0_shift,S0_shift
      real*8  E0,dE0

      logical    precondition
      common / cgsd_block2 / precondition

      real*8     tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / c_cgsd_block / Enew,Eold,Estart

      real*8 tmin,dte,sum0,sum1
      common / c_bfgs_block / tmin,dte,sum0,sum1

      integer it,it_in
      real*8 deltat
      real*8 max_sigma,dt
      real*8  ehfx,phfx

      integer ispin,neq(2),neall,npack1,nbrillq


*     **** external functions ****
      integer  control_it_in,cpsi_neq,cpsi_ispin,cpsi_nbrillq
      integer  cpsi_data_push_stack,cpsi_data_get_allptr
      real*8   control_tole,control_tolc
      real*8   cpsi_geodesic_energy
      real*8   cpsi_geodesic_denergy
      real*8   c_rho_error
      real*8   c_dng_1ehartree
      real*8   cpsi_1ke
      real*8   cpsi_1vl
      real*8   cpsi_1vnl
      real*8   cpsi_1vnlso
      real*8   c_rho_1exc
      real*8   c_rho_1pxc
      real*8   ewald_e
      real*8   cpsi_1eorbit
      real*8   linesearch
   
      external control_it_in,cpsi_neq,cpsi_ispin,cpsi_nbrillq
      external cpsi_data_push_stack,cpsi_data_get_allptr
      external control_tole,control_tolc
      external cpsi_geodesic_energy
      external cpsi_geodesic_denergy
      external c_rho_error
      external c_dng_1ehartree
      external cpsi_1ke
      external cpsi_1vl
      external cpsi_1vnl
      external cpsi_1vnlso
      external c_rho_1exc
      external c_rho_1pxc
      external ewald_e
      external cpsi_1eorbit
      external linesearch


      real*8   control_time_step,control_fake_mass,ion_disp_energy
      external control_time_step,control_fake_mass,ion_disp_energy
      integer  control_lmbfgs_size
      external control_lmbfgs_size
      logical  control_precondition,band_hfx,ion_disp_on
      external control_precondition,band_hfx,ion_disp_on

      dt = control_time_step()
c     dte = dt/dsqrt(control_fake_mass())
c     dte = dsqrt(control_fake_mass())/dt
c     dte = 1.0d0


      if (current_iteration.eq.1) then
         precondition = control_precondition()
      end if

  
      call Cram_max_npack(npack1)
      ispin = cpsi_ispin()
      neq(1) = cpsi_neq(1)
      neq(2) = cpsi_neq(2)
      neall = neq(1)+neq(2)
      nbrillq = cpsi_nbrillq()


*     **** allocate G0, S0 ****
      S0_tag = cpsi_data_push_stack(nbrillq,neall,2*npack1)
      G0_tag = cpsi_data_push_stack(nbrillq,neall,2*npack1)
      S0_shift = cpsi_data_get_allptr(S0_tag)
      G0_shift = cpsi_data_get_allptr(G0_tag)


*     ***** get the initial direction ****
      call cpsi_1get_Tgradient(G0_tag,E0)

      if (current_iteration.eq.1) then

*        ***** use the initial gradient for the direction ****
         call band_lmbfgs_init(control_lmbfgs_size(),dbl_mb(G0_shift))
         call BGrsm_ff_Copy_tag(G0_tag,S0_tag)
         tmin  = 10*deltat_min
      else
         call band_lmbfgs(tmin,dbl_mb(G0_shift),dbl_mb(S0_shift))
         call Pneb_ff_traceall_tag(0,0,G0_tag,S0_tag,sum1)
         if (dabs(sum1).gt.1.0d3) call BGrsm_ff_Copy_tag(G0_tag,S0_tag)
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
         call c_geodesic_start(S0_tag,max_sigma,dE0)


*        ******* line search *********
         if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
c         deltat = 1.0d0  !** debug **
         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        cpsi_geodesic_energy,
     >                        cpsi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
         call cpsi_geodesic_final(tmin)
         deltac = c_rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            if (.not.precondition) go to 30
            precondition = .false.
         end if
     
*        **** get the new gradient - also updates densities****
         call cpsi_2get_Tgradient(2,G0_tag,E0)
         call band_lmbfgs(tmin,dbl_mb(G0_shift),dbl_mb(S0_shift))

         call Pneb_ff_traceall_tag(0,0,G0_tag,S0_tag,sum1)
         if (dabs(sum1).gt.1.0d3) call BGrsm_ff_Copy_tag(G0_tag,S0_tag)

*        **** make psi1 <--- psi2(tmin) ****
         call cpsi_2to1()
c        call c_rho_2to1()
c        call c_dng_2to1()

      end do

*     **** initialize the geoedesic line data structure ****
      call c_geodesic_start(S0_tag,max_sigma,dE0)

*     ******* line search *********
      if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
         deltat = tmin
      else
         deltat = deltat_min
      end if
      Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        cpsi_geodesic_energy,
     >                        cpsi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
      call cpsi_geodesic_final(tmin)
      deltac = c_rho_error()

*     **** free memory ****
 30   call cpsi_data_pop_stack(G0_tag)
      call cpsi_data_pop_stack(S0_tag)
 
      call cpsi_2to1()
      call cpsi_1toelectron() 
c     call c_rho_2to1()
c     call c_dng_2to1()
c      call cpsi_check()

      eion = ewald_e()

      eorbit   = cpsi_1eorbit()
      ehartree = c_dng_1ehartree()
      exc      = c_rho_1exc()
      pxc      = c_rho_1pxc()

      E(1)  = Enew + eion
      E(2)  = eorbit
      E(3)  = ehartree
      E(4)  = exc
      E(5)  = eion
      E(6)  = cpsi_1ke()
      E(7)  = cpsi_1vl()
      E(8)  = cpsi_1vnl()
      E(9)  = 2.0d0*ehartree
      E(10) = pxc
      E(11) = cpsi_1vnlso()


*     **** HFX terms ****
      if (band_HFX()) then
         call c_electron_HFX_energies(ehfx,phfx)
         E(26) = ehfx
         E(27) = phfx
      end if

*     **** Dispersion energy ****
      if (ion_disp_on()) then
         E(33) = ion_disp_energy()
         E(1)  = E(1) + E(33)
      end if


      return
      end
 

