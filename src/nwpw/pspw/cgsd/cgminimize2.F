c
c     $Id: cgminimize2.F,v 1.17 2009-02-07 03:50:32 bylaska Exp $                       
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
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
      integer H0(2),G1(2)
      real*8  E0,dE0

      real*8     sum0(2),sum1(2),scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      integer it,it_in,ms,ispin,ne(2)
      real*8 tmin,deltat
      real*8 max_sigma

      logical value
      integer neall,npack1
      !real*8 e_ionmm,e_qmmm,e_mmmm,e_pol,e_vib,e_cav
      !real*8 e_qmmm_e,e_qmmm_q,e_qmmm_lj,e_mmmm_q,e_mmmm_lj
      real*8 e_lj,e_q,e_spring
      real*8 ehsic,exsic,phsic,pxsic,ehfx,phfx


*     **** external functions ****
      integer  control_it_in,psi_neq,control_version,control_ispin
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
   
      external control_it_in,psi_neq,control_version,control_ispin
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
      real*8   pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E
      external pspw_qmmm_found
      external pspw_qmmm_LJ_E,pspw_qmmm_Q_E,pspw_qmmm_spring_E

*     ***** pspw_charge external functions ****
      logical  pspw_charge_found
      real*8   pspw_charge_Energy_ion,pspw_charge_Energy_charge
      external pspw_charge_found
      external pspw_charge_Energy_ion,pspw_charge_Energy_charge

      logical  control_precondition,pspw_SIC,pspw_HFX
      external control_precondition,pspw_SIC,pspw_HFX
     
  
      call Pack_npack(1,npack1)
      ispin = control_ispin()
      ne(1) = psi_neq(1)
      ne(2) = psi_neq(2)
      neall = ne(1)+ne(2)


*     **** allocate H0 and G1 ****
      value = MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'H0',H0(2),H0(1))
      value = value.and.
     >        MA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('cgminimize2:out of heap memory',0,
     &       MA_ERR)
      call dcopy(2*npack1*neall,0.0d0,0,dcpl_mb(G1(1)),1)

      Estart = Enew

*     ***** get the initial gradient and direction ****
      call psi_1get_TSgradient(dcpl_mb(G1(1)),E0)


      do ms=1,ispin
        call Grsm_gg_trace(npack1,ne(ms),
     >                     dcpl_mb(G1(1)+(ms-1)*npack1*ne(1)),
     >                     dcpl_mb(G1(1)+(ms-1)*npack1*ne(1)),
     >                     sum1(ms))
        call D1dB_SumAll(sum1(ms))
      end do

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
     >                        0.5d0,tmin,deltae,2)
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


         do ms=1,ispin
           sum0(ms)  = sum1(ms)
           call Grsm_gg_trace(npack1,ne(ms),
     >                        dcpl_mb(G1(1)+(ms-1)*npack1*ne(1)),
     >                        dcpl_mb(G1(1)+(ms-1)*npack1*ne(1)),
     >                        sum1(ms))
           call D1dB_SumAll(sum1(ms))
         end do

*        **** the new direction using Fletcher-Reeves ****
         if ( (dabs(deltae).le.(1.0d-2)).and.
     >        (tmin.gt.deltat_min)) then

           do ms=1,ispin
             if (sum0(ms).gt.1.0d-15) then
                scale = sum1(ms)/sum0(ms)  ! Fletcher-Reeves
             else
                scale = 0.0d0
             end if

c             call Grsm_gg_dScale(npack1,ne(ms),scale,
c     >                         dcpl_mb(H0(1)+(ms-1)*npack1*ne(1)),
c     >                         dcpl_mb(H0(1)+(ms-1)*npack1*ne(1)))
             call Grsm_gg_dScale1(npack1,ne(ms),scale,
     >                         dcpl_mb(H0(1)+(ms-1)*npack1*ne(1)))
           end do
c           call Grsm_ggg_Sum(npack1,neall,
c     >                       dcpl_mb(G1(1)),
c     >                       dcpl_mb(H0(1)),
c     >                       dcpl_mb(H0(1)))
           call Grsm_ggg_Sum2(npack1,neall,
     >                       dcpl_mb(G1(1)),
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
     >                        0.5d0,tmin,deltae,2)

      call psi_geodesic2_final(tmin)
      deltac = rho_error()

 30   call psi_2to1()
      call psi_1toelectron() 
c      call psi_check()

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

*     **** HFX terms ****
      if (pspw_HFX()) then
         call electron_HFX_energies(ehfx,phfx)        
         E(26) = ehfx
         E(27) = phfx
      end if

      value = MA_free_heap(G1(2))
      value = value.and.MA_free_heap(H0(2))
      if (.not. value) 
     >  call errquit('cgminimize2:error freeing heap memory',2, MA_ERR)

  
      return
      end
 

