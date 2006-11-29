c
c     $Id: c_bfgsminimize.F,v 1.2 2006-11-29 02:34:54 d3p708 Exp $                       
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
#include "errquit.fh"
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      integer    minimizer

#include "mafdecls.fh"

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
      common / c_cgsd_block / Enew,Eold,Estart

      real*8 tmin,dte,sum0,sum1
      common / c_bfgs_block / tmin,dte,sum0,sum1

      integer it,it_in
      real*8 deltat
      real*8 max_sigma,dt

      logical value
      integer ispin,ne(2),neall,npack1,nbrill




*     **** external functions ****
      integer  brillioun_nbrillioun
      integer  control_it_in,cpsi_ne,cpsi_ispin
      real*8   control_tole,control_tolc
      real*8   cpsi_geodesic_energy
      real*8   cpsi_geodesic_denergy
      real*8   c_rho_error
      real*8   c_dng_1ehartree
      real*8   cpsi_1ke
      real*8   cpsi_1vl
      real*8   cpsi_1vnl
      real*8   c_rho_1exc
      real*8   c_rho_1pxc
      real*8   c_ewald_e
      real*8   cpsi_1eorbit
      real*8   linesearch
   
      external brillioun_nbrillioun
      external control_it_in,cpsi_ne,cpsi_ispin
      external control_tole,control_tolc
      external cpsi_geodesic_energy
      external cpsi_geodesic_denergy
      external c_rho_error
      external c_dng_1ehartree
      external cpsi_1ke
      external cpsi_1vl
      external cpsi_1vnl
      external c_rho_1exc
      external c_rho_1pxc
      external c_ewald_e
      external cpsi_1eorbit
      external linesearch


      real*8   control_time_step,control_fake_mass
      external control_time_step,control_fake_mass
      integer  control_lmbfgs_size
      external control_lmbfgs_size
      logical  control_precondition
      external control_precondition

      dt = control_time_step()
c     dte = dt/dsqrt(control_fake_mass())
c     dte = dsqrt(control_fake_mass())/dt
c     dte = 1.0d0

      write(0,*)"WA"
      if (current_iteration.eq.1) then
         precondition = control_precondition()
      end if

  
      call Cram_max_npack(npack1)
      ispin = cpsi_ispin()
      ne(1) = cpsi_ne(1)
      ne(2) = cpsi_ne(2)
      neall = ne(1)+ne(2)
      nbrill = brillioun_nbrillioun()


*     **** allocate G0, S0 ****
      value = MA_push_get(mt_dcpl,npack1*neall*nbrill,
     >                     'S0',S0(2),S0(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1*neall*nbrill,
     >                     'G0',G0(2),G0(1))
      if (.not. value) 
     > call errquit('c_bfgsminimize:out of heap memory',0,MA_ERR)


      write(0,*)"ND"

*     ***** get the initial direction ****
      call cpsi_1get_Tgradient(dcpl_mb(G0(1)),E0)

      write(0,*)"tx"
      if (current_iteration.eq.1) then

*        ***** use the initial gradient for the direction ****
         write(0,*)"lmbfgs_size= ",control_lmbfgs_size() 
         call band_lmbfgs_init(control_lmbfgs_size(),dcpl_mb(G0(1)))
         
         write(0,*)"ar"
         call Pneb_ff_copy(0,0,npack1,
     >                      dcpl_mb(G0(1)),
     >                      dcpl_mb(S0(1)))
         tmin  = 10*deltat_min
      else
         call band_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
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

         write(0,*)"it = ",it
         write(0,*)"looper1"

*        **** initialize the geoedesic line data structure ****
         call c_geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)
         write(0,*)"ut"

*        ******* line search *********
         if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
c         deltat = 1.0d0  !** debug **
         write(0,*)"looper2"
         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        cpsi_geodesic_energy,
     >                        cpsi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
         write(0,*)"looper3"
         call cpsi_geodesic_final(tmin)
         write(0,*)"looper4"
         deltac = c_rho_error()
         write(0,*)"ca"

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            if (.not.precondition) go to 30
            precondition = .false.
         end if
     
*        **** get the new gradient - also updates densities****
         write(0,*)"looper5"
         call cpsi_2get_Tgradient(2,dcpl_mb(G0(1)),E0)
         write(0,*)"looper6"
cccccccccc crashes here !!!!!!!!!!!!!!
         call band_lmbfgs(tmin,dcpl_mb(G0(1)),dcpl_mb(S0(1)))
ccccccccccccccccc
         write(0,*)"looper7"


      write(0,*)"ny"

*        **** make psi1 <--- psi2(tmin) ****
         call cpsi_2to1()
         write(0,*)"looper8"
c        call c_rho_2to1()
c        call c_dng_2to1()

      end do

*     **** initialize the geoedesic line data structure ****
         write(0,*)"looper9"
      call c_geodesic_start(dcpl_mb(S0(1)),max_sigma,dE0)
         write(0,*)"looper10"

*     ******* line search *********
      if ((tmin.gt.deltat_min).and.(tmin.lt.1.0d4)) then
         deltat = tmin
      else
         deltat = deltat_min
      end if
         write(0,*)"looper11"
      Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        cpsi_geodesic_energy,
     >                        cpsi_geodesic_denergy,
     >                        0.50d0,tmin,deltae,2)
         write(0,*)"looper12"
      call cpsi_geodesic_final(tmin)
         write(0,*)"looper13"
      deltac = c_rho_error()
         write(0,*)"looper14"

*     **** free memory ****
 30   value = MA_pop_stack(G0(2))
      value = value.and.
     >        MA_pop_stack(S0(2))
      if (.not. value) 
     > call errquit('c_bfgsminimize:error freeing heap',0,MA_ERR)

      call cpsi_2to1()
      call cpsi_1toelectron() 
c     call c_rho_2to1()
c     call c_dng_2to1()
      call cpsi_check()

      write(0,*)"ma"
      
      eion = c_ewald_e()

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

      return
      end
 

