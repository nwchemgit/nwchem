c
c     $Id: cgminimize2.F,v 1.6 2003-10-17 22:57:10 carlfahl Exp $                       
c

*  **************************************************************
*  *                MPI cgminimize2 routine                     *
*  *               (Fletcher-Reeves' steps)                     *
*  *                                                            *
*  *  This is a developing Stiefel conjugate gradient minimizer *
*  *  written for NWChem.                                       *
*  *                                                            *
*  *                                                            *
*  **************************************************************
  
      subroutine cgminimize2(E,deltae,deltac,current_iteration)
      implicit none
#include "errquit.fh"
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration

#include "mafdecls.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer H0(2),G1(2)
      real*8  E0,dE0

      real*8     sum0,sum1,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      integer it,it_in
      real*8 tmin,deltat
      real*8 max_sigma

      logical value
      integer neall,npack1
      real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      real*8 e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      real*8 ehsic,exsic,phsic,pxsic      


*     **** external functions ****
      integer  control_it_in,psi_ne,control_version
      real*8   control_tole,control_tolc
      real*8   psi_geodesic2_energy
      real*8   psi_geodesic2_denergy
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
   
      external control_it_in,psi_ne,control_version
      external control_tole,control_tolc
      external psi_geodesic2_energy
      external psi_geodesic2_denergy
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
      real*8   pspw_qmmm_Energy_QMMM_Q
      real*8   pspw_qmmm_Energy_QMMM_LJ
      real*8   pspw_qmmm_Energy_MMMM_Q
      real*8   pspw_qmmm_Energy_MMMM_LJ
      real*8   pspw_qmmm_Energy_pol
      real*8   pspw_qmmm_Energy_vib
      real*8   pspw_qmmm_Energy_cav
      external pspw_qmmm_found
      external pspw_qmmm_Energy_QMMM_Q
      external pspw_qmmm_Energy_QMMM_LJ
      external pspw_qmmm_Energy_MMMM_Q
      external pspw_qmmm_Energy_MMMM_LJ
      external pspw_qmmm_Energy_pol
      external pspw_qmmm_Energy_vib
      external pspw_qmmm_Energy_cav

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

      logical  control_precondition,pspw_SIC
      external control_precondition,pspw_SIC
     
  
      call Pack_npack(1,npack1)
      neall = psi_ne(1)+psi_ne(2)


*     **** allocate H0 and G1 ****
      value = MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'H0',H0(2),H0(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('cgminimize2:out of heap memory',0,
     &       MA_ERR)

      Estart = Enew

*     ***** get the initial gradient and direction ****
      call psi_1get_TSgradient(dcpl_mb(G1(1)),E0)


      call Grsm_gg_trace(npack1,neall,
     >                   dcpl_mb(G1(1)),
     >                   dcpl_mb(G1(1)),sum1)

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
         call psi_1geodesic2_start(dcpl_mb(H0(1)),max_sigma,dE0)


*        ******* line search *********
         if (tmin.gt.deltat_min) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        psi_geodesic2_energy,
     >                        psi_geodesic2_denergy,
     >                        tole,tmin,deltae,1)
         call psi_geodesic2_final(tmin)
         deltac = rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            go to 30
         end if
   

*        **** transport the previous search directions ****
         call psi_1geodesic2_transport(tmin,dcpl_mb(H0(1)))

*        **** make psi1 <--- psi2(tmin) ****
         call psi_2to1()

*        **** get the new gradient - also updates densities****
         call psi_1get_TSgradient(dcpl_mb(G1(1)),E0)


         sum0  = sum1
         call Grsm_gg_trace(npack1,neall,
     >                        dcpl_mb(G1(1)),
     >                        dcpl_mb(G1(1)),
     >                        sum1)

*        **** the new direction using Fletcher-Reeves ****
         if ( (dabs(deltae).le.(1.0d-2)).and.
     >        (tmin.gt.deltat_min)) then

           if (sum0.gt.1.0d-15) then
              scale = sum1/sum0  ! Fletcher-Reeves
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


*        **** the new direction using steepest-descent ****
         else
              call Grsm_gg_Copy(npack1,neall,
     >                          dcpl_mb(G1(1)),
     >                          dcpl_mb(H0(1)))
         end if

      end do


*     **** initialize the geoedesic line data structure ****
      call psi_1geodesic2_start(dcpl_mb(H0(1)),max_sigma,dE0)


*     ******* line search *********
      if (tmin.gt.deltat_min) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

      Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        psi_geodesic2_energy,
     >                        psi_geodesic2_denergy,
     >                        tole,tmin,deltae,1)

      call psi_geodesic2_final(tmin)
      deltac = rho_error()

 30   call psi_2to1()
      call psi_1toelectron() 
      call psi_check()

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

*        **** get QM/MM  energies ****
         e_qmmm_e  = psi_1v_field()
         e_qmmm_q   = pspw_qmmm_Energy_QMMM_Q()
         e_qmmm_lj  = pspw_qmmm_Energy_QMMM_LJ()

         e_mmmm_q  = pspw_qmmm_Energy_MMMM_Q()
         e_mmmm_lj  = pspw_qmmm_Energy_MMMM_LJ()

         e_pol   = pspw_qmmm_Energy_pol()
         e_vib   = pspw_qmmm_Energy_vib()
         e_cav   = pspw_qmmm_Energy_cav()

         E(1) = E(1) + e_qmmm_q
     >               + e_qmmm_lj
     >               + e_mmmm_q
     >               + e_mmmm_lj
     >               + e_pol
     >               + e_vib
     >               + e_cav
         E(11) = e_qmmm_e
         E(12) = e_qmmm_q
         E(13) = e_qmmm_lj
         E(14) = e_mmmm_q
         E(15) = e_mmmm_lj
         E(16) = e_pol
         E(17) = e_vib
         E(18) = e_cav
      end if



*     **** get pspw_charge  energies ****
      if (pspw_charge_found()) then
         E(19)  = psi_1v_field()
         E(20)  = pspw_charge_Energy_ion()
         E(21)  = pspw_charge_Energy_charge()
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

      value = MA_free_heap(G1(2))
      value = value.and.MA_free_heap(H0(2))
      if (.not. value) 
     >  call errquit('cgminimize2:error freeing heap memory',2, MA_ERR)

  
      return
      end
 

