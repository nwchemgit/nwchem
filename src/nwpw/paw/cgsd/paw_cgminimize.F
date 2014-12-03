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
  
      subroutine paw_cgminimize(E,deltae,deltac,current_iteration,it_in)
      implicit none
      real*8     E(*)
      real*8     deltae,deltac
      integer    current_iteration
      integer    it_in


#include "bafdecls.fh"
#include "errquit.fh"
#include "paw_energy_kin_atom.fh"
#include "paw_energy_vloc_atom.fh"
#include "paw_energy_ion_atom.fh"
#include "paw_energy_core_atom.fh"
#include "paw_energy_hartree_atom.fh"
#include "paw_xc.fh"


*     **** local variables ****
      real*8  deltat_min
      parameter (deltat_min=1.0d-3)
       
c      complex*16 H0(nfft3d,nemax)
c      complex*16 G1(nfft3d,nemax)
      integer H0(2),G1(2),R1(2)
      real*8  E0,dE0

      logical    precondition
      common / paw_cgsd_block2 / precondition

      real*8     sum0,sum1,sum3,scale,tole,tolc
      real*8     ehartree,eorbit,exc,pxc,eion
      real*8     Enew,Eold,Estart
      common / paw_cgsd_block / Enew,Eold,Estart

      integer it
      real*8 tmin,deltat
      real*8 max_sigma,w

      logical value
      integer i,neall,npack1,k

      double precision evloc_pw,evloc_atom,ehartr_pw,eke
      double precision ekin_atom
      double precision eion_atom
      double precision ecore_atom
      double precision ecore_ion_atom
      double precision ecore_self_atom
      double precision ehartree_atom
      double precision exc_atom,pxc_atom



*     **** external functions ****
      integer  paw_psi_ne,control_version
      real*8   control_tole,control_tolc
      real*8   paw_psi_geodesic_energy
      real*8   paw_psi_geodesic_denergy
      real*8   paw_rho_error
      real*8   paw_dng_1ehartree
      real*8   paw_psi_1ke
      real*8   paw_psi_1vl
      real*8   paw_rho_1exc
      real*8   paw_rho_1pxc
      real*8   paw_psi_1eorbit
      real*8   linesearch
      real*8   paw_mult_energy_atom_comp !**no header file for paw_mult**
      real*8   paw_mult_energy_atom_self
      real*8   paw_mult_energy_atom_mult
   
      external paw_psi_ne,control_version
      external control_tole,control_tolc
      external paw_psi_geodesic_energy
      external paw_psi_geodesic_denergy
      external paw_rho_error
      external paw_dng_1ehartree
      external paw_psi_1ke
      external paw_psi_1vl
      external paw_rho_1exc
      external paw_rho_1pxc
      external paw_psi_1eorbit
      external linesearch
      external paw_mult_energy_atom_comp
      external paw_mult_energy_atom_self
      external paw_mult_energy_atom_mult


      call Pack_npack(1,npack1)
      neall = paw_psi_ne(1)+paw_psi_ne(2)

*     **** allocate H0, G1, and R1 ****
      value = BA_alloc_get(mt_dcpl,npack1*neall,
     >                     'H0',H0(2),H0(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*neall,
     >                     'G1',G1(2),G1(1))
      value = value.and.
     >        BA_alloc_get(mt_dcpl,npack1*neall,
     >                     'R1',R1(2),R1(1))
      if (.not. value) 
     >  call errquit('paw_cgminimize:out of heap memory',0,MA_ERR)

      Estart = Enew


*     ***** get the initial gradient and direction ****
      call paw_psi_1get_residual(dcpl_mb(R1(1)),dcpl_mb(G1(1)),E0)
c      call ke_Precondition2(npack1,neall,dcpl_mb(R1(1)),dcpl_mb(G1(1)))
      call dcopy(2*neall*npack1,dcpl_mb(R1(1)),1,dcpl_mb(G1(1)),1)
      call paw_psi_1gen_Tangent(dcpl_mb(G1(1)))


      call Grsm_gg_trace(npack1,neall,
     >                   dcpl_mb(G1(1)),
     >                   dcpl_mb(R1(1)),sum1)

      call Grsm_gg_Copy(npack1,neall,
     >                   dcpl_mb(G1(1)),
     >                   dcpl_mb(H0(1)))

*     ******************************************
*     ****                                  ****
*     **** Start of conjugate gradient loop ****
*     ****                                  ****
*     ******************************************
      !it_in = control_it_in()
      tole  = control_tole()
      tolc  = control_tolc()
      tmin  = deltat_min
      do it=2,it_in

*        **** initialize the geoedesic line data structure ****
         call paw_geodesic_start(dcpl_mb(H0(1)),max_sigma,dE0)


*        ******* line search *********
         if (tmin.gt.deltat_min) then
            deltat = tmin
         else
            deltat = deltat_min
         end if
         Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        paw_psi_geodesic_energy,
     >                        paw_psi_geodesic_denergy,
     >                        tole,tmin,deltae,1)

         call paw_psi_geodesic_final(tmin)
         deltac = paw_rho_error()

*        **** exit loop early ****
         if ((dabs(deltae).lt.tole).and.(deltac.lt.tolc)) then
            go to 30
         end if
   

*        **** transport the previous search directions ****
         call paw_psi_1geodesic_transport(tmin,dcpl_mb(H0(1)))
         

*        **** make psi1 <--- psi2(tmin) ****
         call paw_psi_2to1()

*        **** get the new gradient - also updates densities****
         call paw_psi_1get_residual(dcpl_mb(R1(1)),dcpl_mb(G1(1)),E0)
c         call ke_Precondition2(npack1,neall,
c     >                         dcpl_mb(R1(1)),
c     >                         dcpl_mb(G1(1)))
         call dcopy(2*neall*npack1,dcpl_mb(R1(1)),1,dcpl_mb(G1(1)),1)
         call paw_psi_1gen_Tangent(dcpl_mb(G1(1)))
      
      

         sum0  = sum1
         call Grsm_gg_trace(npack1,neall,
     >                        dcpl_mb(G1(1)),
     >                        dcpl_mb(R1(1)),
     >                        sum1)


*        **** the new direction using Fletcher-Reeves ****
         if ( (dabs(deltae).le.(1.0d-2)).and.
     >        (tmin.gt.deltat_min)) then

           if (sum0.gt.1.0d-15) then
c              scale = sum3/sum0
              scale = sum1/sum0
           else
              scale = 0.0d0
           end if

c           call Grsm_gg_dScale(npack1,neall,scale,
c     >                         dcpl_mb(H0(1)),
c     >                         dcpl_mb(H0(1)))
c           call Grsm_ggg_Sum(npack1,neall,
c     >                       dcpl_mb(G1(1)),
c     >                       dcpl_mb(H0(1)),
c     >                       dcpl_mb(H0(1)))
           call Grsm_gg_dScale1(npack1,neall,scale,dcpl_mb(H0(1)))
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
      call paw_geodesic_start(dcpl_mb(H0(1)),max_sigma,dE0)


*     ******* line search *********
      if (tmin.gt.deltat_min) then
         deltat = tmin
      else
         deltat = deltat_min
      end if

      Enew = linesearch(0.0d0,E0,dE0,deltat,
     >                        paw_psi_geodesic_energy,
     >                        paw_psi_geodesic_denergy,
     >                        tole,tmin,deltae,1)

      call paw_psi_geodesic_final(tmin)
      deltac = paw_rho_error()


 30   call paw_psi_2to1()
      call paw_psi_1toelectron() 
c      call paw_psi_check()


*     **** plane-wave energies ****
      eorbit    = paw_psi_1eorbit()
      eke       = paw_psi_1ke()
      ehartr_pw = paw_dng_1ehartree()
      exc       = paw_rho_1exc()
      pxc       = paw_rho_1pxc()
      evloc_pw  = paw_psi_1vl()


*     **** atomic energies ***
      ekin_atom     = paw_energy_kin_atom()
      evloc_atom    = paw_energy_vloc_atom()

      ehartree_atom   = paw_energy_hartree_atom()
      eion_atom       = paw_energy_ion_atom()
      ecore_atom      = paw_energy_core_atom()
      ecore_ion_atom  = paw_energy_core_ion_atom()
      ecore_self_atom = paw_energy_core_self_atom()

      exc_atom        = paw_energy_xc_atom()



*     *** kinetic energy
      E(2) = eke
      E(3) = ekin_atom

*     *** coulomb contributions
      E(4) = ehartr_pw
      E(5) = eion_atom + ecore_atom + ehartree_atom +
     >       ecore_ion_atom + ecore_self_atom +
     >       paw_mult_energy_atom_self() +
     >       paw_mult_energy_atom_comp()

      E(6)=paw_mult_energy_atom_mult()

*     *** exch-correlation
      E(7) = exc
      E(8) = exc_atom

*     *** local pseudopot ***
      E(9)  = evloc_pw
      E(10) = evloc_atom


*     *** total energy ***
      Eold=E(1)
      E(1) = 0.0d0
      do i=2,10
       E(1) = E(1) + E(i)
      end do

      E(11) = eorbit


      value = BA_free_heap(R1(2))
      value = value.and.BA_free_heap(G1(2))
      value = value.and.BA_free_heap(H0(2))
      if (.not. value) 
     >  call errquit('cgminimize:error freeing heap memory',0, MA_ERR)

  
      return
      end
 

