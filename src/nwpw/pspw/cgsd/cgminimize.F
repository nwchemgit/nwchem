c     $Id: cgminimize.F,v 1.7 2001-12-27 20:54:54 bylaska Exp $                       
*  ************************************************************
*  *                MPI cgminimize routine                    *
*  *               (Fletcher-Reeves' steps)                   *
*  *                                                          *
*  *  This is a developing cgsdv3 parallel code wrtten for    *
*  *  NWChem                                                  *
*  *                                                          *
*  ************************************************************
  
      subroutine cgminimize(E,deltae,deltac,current_iteration)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration

#include "mafdecls.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
c      complex*16 H0(nfft3d,nemax)
c      complex*16 G1(nfft3d,nemax)
      integer H0(2),G0(2),G1(2)

      logical    precondition
      common / cgsd_block2 / precondition

      real*8     sum0,sum1,sum3,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      integer it,it_in
      real*8 tmin,deltat
      real*8 max_sigma

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

      logical  control_precondition
      external control_precondition
       
  
      call Pack_npack(1,npack1)
      neall = psi_ne(1)+psi_ne(2)


*     **** allocate H0, G0, and G1 ****
      value = MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'H0',H0(2),H0(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G0',G0(2),G0(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('cgminimize:out of heap memory',0)

      Estart = Enew
      if (current_iteration.eq.1) then
         precondition = control_precondition()
      end if

*     ***** get the initial gradient and direction ****
      if (precondition) then
         call psi_1get_TMgradient(dcpl_mb(G1(1)))
      else
         call psi_1get_Tgradient(dcpl_mb(G1(1)))
      end if

      call Grsm_gg_trace(npack1,neall,
     >                   dcpl_mb(G1(1)),
     >                   dcpl_mb(G1(1)),sum1)

c     write(*,*) "sum1a,Estart,:",sum1,Estart
      call Grsm_gg_Copy(npack1,neall,
     >                   dcpl_mb(G1(1)),
     >                   dcpl_mb(H0(1)))

*     ******************************************
*     ****                                  ****
*     **** Start of conjugate gradient loop ****
*     ****                                  ****
*     ******************************************
      it_in = control_it_in()
      tole  = control_tole()
      tolc  = control_tolc()
      tmin  = deltat_min
      do it=2,it_in

*        **** initialize the geoedesic line data structure ****
         call geodesic_start(dcpl_mb(H0(1)),max_sigma)
c        Eold= psi_geodesic_energy(0.0d0)
c        write(*,*) "Eold:", Eold,max_sigma

*        ******* line search *********
         if (tmin.gt.deltat_min) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
c        Eold = psi_geodesic_energy(0.0d0)
         Enew = linesearch(0.0d0,deltat,
     >                        psi_geodesic_energy,
     >                        psi_geodesic_denergy,
     >                        tole,tmin,deltae)
         call psi_geodesic_final(tmin)
         deltac = rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            if (.not.precondition)  go to 30
            precondition = .false.
         end if
   

*        **** transport the previous search directions ****
         call Grsm_gg_Copy(npack1,neall,
     >                     dcpl_mb(G1(1)),
     >                     dcpl_mb(G0(1)))

         call psi_1geodesic_transport(tmin,dcpl_mb(H0(1)))
         call psi_1geodesic_Gtransport(tmin,dcpl_mb(G0(1)))

*        **** make psi1 <--- psi2(tmin) ****
         call psi_2to1()
c        call rho_2to1()
c        call dng_2to1()

*        **** get the new gradient - also updates densities****
         if (precondition) then
           call psi_1get_TMgradient(dcpl_mb(G1(1)))
         else
           call psi_1get_Tgradient(dcpl_mb(G1(1)))
         end if

         call Grsm_ggg_Sub(npack1,neall,
     >                       dcpl_mb(G1(1)),
     >                       dcpl_mb(G0(1)),
     >                       dcpl_mb(G0(1)))

         sum0  = sum1
         call Grsm_gg_trace(npack1,neall,
     >                        dcpl_mb(G1(1)),
     >                        dcpl_mb(G1(1)),
     >                        sum1)
         call Grsm_gg_trace(npack1,neall,
     >                        dcpl_mb(G0(1)),
     >                        dcpl_mb(G1(1)),
     >                        sum3)

c        write(*,*) "sum1:",sum1,sum0,sum3,max_sigma
c        write(*,*) "Enew:",Enew

*        **** the new direction using Fletcher-Reeves ****
         if ( (dabs(deltae).le.(1.0d-2)).and.
     >        (tmin.gt.deltat_min)) then

           if (sum0.gt.1.0d-15) then
c              scale = sum3/sum0
              scale = sum1/sum0
           else
              scale = 0.0d0
           end if

           call Grsm_gg_dScale(npack1,neall,scale,
     >                         dcpl_mb(H0(1)),
     >                         dcpl_mb(H0(1)))
           call Grsm_ggg_Sum(npack1,neall,
     >                       dcpl_mb(G1(1)),
     >                       dcpl_mb(H0(1)),
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
      call geodesic_start(dcpl_mb(H0(1)),max_sigma)

*     ******* line search *********
      if (tmin.gt.deltat_min) then
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

 30   call psi_2to1()
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

      value = MA_free_heap(G1(2))
      value = value.and.MA_free_heap(G0(2))
      value = value.and.MA_free_heap(H0(2))
      if (.not. value) 
     >  call errquit('cgminimize:error freeing heap memory',0)

  
      return
      end
 

