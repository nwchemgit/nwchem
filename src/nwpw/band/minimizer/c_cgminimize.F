*
*     $Id: c_cgminimize.F,v 1.16 2006-12-16 01:54:00 bylaska Exp $                       
*


*  ************************************************************
*  *                MPI cgminimize routine                    *
*  *               (Fletcher-Reeves' steps)                   *
*  *                                                          *
*  *  This is a developing cgsdv3 parallel code wrtten for    *
*  *  NWChem                                                  *
*  *                                                          *
*  ************************************************************
  
      subroutine c_cgminimize(E,deltae,deltac,current_iteration,it_in)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      integer    it_in

#include "mafdecls.fh"
#include "errquit.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
c      complex*16 H0(nfft3d,nemax)
c      complex*16 G1(nfft3d,nemax)
      integer H0(2),G0(2),G1(2)
      real*8  E0,dE0

      real*8     sum0,sum1,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / c_cgsd_block / Enew,Eold,Estart

      integer it
      real*8 tmin,deltat
      real*8 max_sigma

      logical value
      integer ispin,ne(2)
      integer neall,nfft3d,npack1,nbrill


*     **** external functions ****
      integer  brillioun_nbrillioun
      integer  cpsi_neq,control_version,cpsi_ispinq,cpsi_nbrillq
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
      logical  cpsi_spin_orbit
      external brillioun_nbrillioun
      external cpsi_neq,control_version,cpsi_ispinq,cpsi_nbrillq
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
      external cpsi_spin_orbit  
       
  
      call C3dB_nfft3d(1,nfft3d)
      call Cram_max_npack(npack1)
      ispin  = cpsi_ispinq()
      ne(1)  = cpsi_neq(1)
      ne(2)  = cpsi_neq(2)
      neall  = ne(1)+ne(2)
      nbrill = cpsi_nbrillq()


*     **** allocate H0, G0, and G1 ****
      value = MA_push_get(mt_dcpl,npack1*neall*nbrill,
     >                     'H0',H0(2),H0(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1*neall*nbrill,
     >                     'G0',G0(2),G0(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,npack1*neall*nbrill,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('c_cgminimize:out of heap memory',0,
     &       MA_ERR)


      Estart = Enew
      Eold= Enew

*     ***** get the initial gradient and direction ****
      call cpsi_1get_Tgradient(dcpl_mb(G1(1)),E0)


      call Pneb_ff_traceall(0,0,
     >                      dcpl_mb(G1(1)),
     >                      dcpl_mb(G1(1)),
     >                      npack1,sum1)

      call BGrsm_ff_Copy(dcpl_mb(G1(1)),
     >                   dcpl_mb(H0(1)))

 

*     ******************************************
*     ****                                  ****
*     **** Start of conjugate gradient loop ****
*     ****                                  ****
*     ******************************************

      tole  = control_tole()
      tolc  = control_tolc()
      tmin  = deltat_min

      do it=2,it_in

*        **** initialize the geoedesic line data structure ****
         call c_geodesic_start(dcpl_mb(H0(1)),max_sigma,dE0)


*        ******* line search *********
         if (tmin.gt.0.0d0) then
            deltat = tmin
         else
            deltat = deltat_min
         end if

         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        cpsi_geodesic_energy,
     >                        cpsi_geodesic_denergy,
     >                        tole,tmin,deltae,1)

c         write(*,*)'e0= ',e0
c         write(*,*)'de0 = ',de0
c         write(*,*)'tmin = ',tmin
c         call flush(6)
         call cpsi_geodesic_final(tmin)

         deltac = c_rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) 
     >      go to 30
     
   

*        **** transport the previous search directions ****
         call BGrsm_ff_Copy(dcpl_mb(G1(1)),
     >                      dcpl_mb(G0(1)))


         call cpsi_1geodesic_transport(tmin,dcpl_mb(H0(1)))
         call cpsi_1geodesic_Gtransport(tmin,dcpl_mb(G0(1)))

*        **** make psi1 <--- psi2(tmin) ****
         call cpsi_2to1()

*        **** get the new gradient - also updates densities****
         call cpsi_1get_Tgradient(dcpl_mb(G1(1)),E0)

         call BGrsm_fff_Sub(dcpl_mb(G1(1)),
     >                      dcpl_mb(G0(1)),
     >                      dcpl_mb(G0(1)))

         sum0  = sum1

         call Pneb_ff_traceall(0,0,
     >                        dcpl_mb(G1(1)),
     >                        dcpl_mb(G1(1)),npack1,
     >                        sum1)

*        **** the new direction using Fletcher-Reeves ****
         if (dabs(deltae).le.(1.0d-2).and.(tmin.gt.0.0d0)) then

           if (sum0.gt.1.0d-15) then
              scale = sum1/sum0
           else
              scale = 0.0d0
           end if

           call BGrsm_ff_dscale(scale,
     >                          dcpl_mb(H0(1)),
     >                          dcpl_mb(H0(1)))

           call BGrsm_fff_sum(dcpl_mb(G1(1)),
     >                        dcpl_mb(H0(1)),
     >                        dcpl_mb(H0(1)))

*        **** the new direction using steepest-descent ****
         else
              call BGrsm_ff_Copy(dcpl_mb(G1(1)),
     >                           dcpl_mb(H0(1)))

         end if

      end do

*     **** initialize the geoedesic line data structure ****

      call c_geodesic_start(dcpl_mb(H0(1)),max_sigma,dE0)

*     ******* line search *********
      if (tmin.gt.0.0d0) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

      Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        cpsi_geodesic_energy,
     >                        cpsi_geodesic_denergy,
     >                        tole,tmin,deltae,1)

      call cpsi_geodesic_final(tmin)

      deltac = c_rho_error()

*     **** free memory ****
 30   value = MA_pop_stack(G1(2))
      value = value.and.MA_pop_stack(G0(2))
      value = value.and.MA_pop_stack(H0(2))
      if (.not. value) 
     >  call errquit('c_cgminimize:error freeing stack memory',0,
     &       MA_ERR)
      
      

      call cpsi_2to1()


      call cpsi_1toelectron() 


      call cpsi_check()

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
 

