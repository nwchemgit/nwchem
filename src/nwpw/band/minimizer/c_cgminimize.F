*
*     $Id: c_cgminimize.F,v 1.3 2001-12-21 16:50:29 bylaska Exp $                       
*


*  ************************************************************
*  *                MPI cgminimize routine                    *
*  *               (Fletcher-Reeves' steps)                   *
*  *                                                          *
*  *  This is a developing cgsdv3 parallel code wrtten for    *
*  *  NWChem                                                  *
*  *                                                          *
*  ************************************************************
  
      subroutine c_cgminimize(E,deltae,deltac)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac

#include "mafdecls.fh"

*     **** local variables ****

      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
c      complex*16 H0(nfft3d,nemax)
c      complex*16 G1(nfft3d,nemax)
      integer H0(2),G0(2),G1(2)

      real*8     sum0,sum1,sum3,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / cgsd_block / Enew,Eold,Estart

      integer it,it_in
      real*8 tmin,deltat
      real*8 max_sigma

      logical value
      integer ispin,ne(2)
      integer neall,nfft3d,nbrill


*     **** external functions ****
      integer  brillioun_nbrillioun
      integer  control_it_in,cpsi_ne,control_version,cpsi_ispin
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
      external control_it_in,cpsi_ne,control_version,cpsi_ispin
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

       
  
      call C3dB_nfft3d(1,nfft3d)
      ispin  = cpsi_ispin()
      ne(1)  = cpsi_ne(1)
      ne(2)  = cpsi_ne(2)
      neall  = ne(1)+ne(2)
      nbrill = brillioun_nbrillioun()


*     **** allocate H0, G0, and G1 ****
      value = MA_push_get(mt_dcpl,nfft3d*neall*nbrill,
     >                     'H0',H0(2),H0(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d*neall*nbrill,
     >                     'G0',G0(2),G0(1))
      value = value.and.
     >        MA_push_get(mt_dcpl,nfft3d*neall*nbrill,
     >                     'G1',G1(2),G1(1))
      if (.not. value) call errquit('c_cgminimize:out of heap memory',0)

      Estart = Enew



*     ***** get the initial gradient and direction ****
      call cpsi_1get_Tgradient(dcpl_mb(G1(1)))
      call Stfl_gg_traceall(ispin,ne,nfft3d,neall,nbrill,
     >                   dcpl_mb(G1(1)),
     >                   dcpl_mb(G1(1)),sum1)

c     write(*,*) "sum1a,Estart,:",sum1,Estart
      call Stfl_gg_Copy(ispin,ne,nfft3d,neall,nbrill,
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
         call c_geodesic_start(dcpl_mb(H0(1)),max_sigma)
c        Eold= cpsi_geodesic_energy(0.0d0)
c        write(*,*) "Eold:", Eold,max_sigma

*        ******* line search *********
         if (tmin.gt.0.0d0) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
c        Eold = cpsi_geodesic_energy(0.0d0)
         Enew = linesearch(0.0d0,deltat,
     >                        cpsi_geodesic_energy,
     >                        cpsi_geodesic_denergy,
     >                        tole,tmin,deltae)
         call cpsi_geodesic_final(tmin)
         deltac = c_rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) 
     >      go to 30
     
   

*        **** transport the previous search directions ****
         call Stfl_gg_Copy(ispin,ne,nfft3d,neall,nbrill,
     >                     dcpl_mb(G1(1)),
     >                     dcpl_mb(G0(1)))

         call cpsi_1geodesic_transport(tmin,dcpl_mb(H0(1)))
         call cpsi_1geodesic_Gtransport(tmin,dcpl_mb(G0(1)))

*        **** make psi1 <--- psi2(tmin) ****
         call cpsi_2to1()
c        call c_rho_2to1()
c        call c_dng_2to1()

*        **** get the new gradient - also updates densities****
         call cpsi_1get_Tgradient(dcpl_mb(G1(1)))

         call Stfl_ggg_Sub(ispin,ne,nfft3d,neall,nbrill,
     >                       dcpl_mb(G1(1)),
     >                       dcpl_mb(G0(1)),
     >                       dcpl_mb(G0(1)))

         sum0  = sum1
         call Stfl_gg_traceall(ispin,ne,nfft3d,neall,nbrill,
     >                        dcpl_mb(G1(1)),
     >                        dcpl_mb(G1(1)),
     >                        sum1)
         call Stfl_gg_traceall(ispin,ne,nfft3d,neall,nbrill,
     >                        dcpl_mb(G0(1)),
     >                        dcpl_mb(G1(1)),
     >                        sum3)

c        write(*,*) "sum1:",sum1,sum0,sum3,max_sigma
c        write(*,*) "Enew:",Enew

*        **** the new direction using Fletcher-Reeves ****
         if (dabs(deltae).le.(1.0d-2).and.(tmin.gt.0.0d0)) then

           if (sum0.gt.1.0d-15) then
c              scale = sum3/sum0
              scale = sum1/sum0
           else
              scale = 0.0d0
           end if

           call Stfl_gg_dScale(ispin,ne,nfft3d,neall,nbrill,
     >                         scale,
     >                         dcpl_mb(H0(1)),
     >                         dcpl_mb(H0(1)))
           call Stfl_ggg_Sum(ispin,ne,nfft3d,neall,nbrill,
     >                       dcpl_mb(G1(1)),
     >                       dcpl_mb(H0(1)),
     >                       dcpl_mb(H0(1)))


*          ***** set to gradient if direction is assending *****
*          call Grsm_gg_traceall(neall,H0,G1,kappa)
*          if (kappa.le.0.0d0) call Grsm_gg_Copy(neall,G1,H0)

*        **** the new direction using steepest-descent ****
         else
              call Stfl_gg_Copy(ispin,ne,nfft3d,neall,nbrill,
     >                          dcpl_mb(G1(1)),
     >                          dcpl_mb(H0(1)))
         end if

      end do

*     **** initialize the geoedesic line data structure ****
      call c_geodesic_start(dcpl_mb(H0(1)),max_sigma)

*     ******* line search *********
      if (tmin.gt.0.0d0) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

c      Eold = cpsi_geodesic_energy(0.0d0)
c     Eold = cpsi_geodesic_energy(0.0d0)
c     Enew = Eold
c     tmin = 0.0d0
      Enew = linesearch(0.0d0,deltat,
     >                        cpsi_geodesic_energy,
     >                        cpsi_geodesic_denergy,
     >                        tole,tmin,deltae)
      call cpsi_geodesic_final(tmin)
      deltac = c_rho_error()

*     **** free memory ****
      value = MA_pop_stack(G1(2))
      value = value.and.MA_pop_stack(G0(2))
      value = value.and.MA_pop_stack(H0(2))
      if (.not. value) 
     >  call errquit('c_cgminimize:error freeing stack memory',0)

 30   call cpsi_2to1()
      call cpsi_1toelectron() 
c     call c_rho_2to1()
c     call c_dng_2to1()
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
 

